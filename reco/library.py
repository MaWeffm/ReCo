# pylint: disable=line-too-long

"""Library class for ReCo."""

from collections import Counter
from pathlib import Path
from subprocess import check_output

import os
import uuid

import pandas as pd

from .attribute_descriptors import LoggerDesc
from .tools import alphabet
from .tools import read_df
from .tools import remove_duplicates
from .tools import reverse_complement

from .reco_config import Config


class Library:
    """Class representing a gRNA library.

    Parameters
    ----------
    logger : logging.Logger
        A logging.Logger object.
    library_file : str
        Path to the library file in one of the following formats:
        "*.xlsx", "*.csv", "*.tsv", "*.txt".

    Methods
    -------
    from_file
        Directly instantiate by reading a library file.
    read_library_file
        Read a library file in one of the following formats:
        "*.xlsx", "*.csv", "*.tsv", "*.txt".
    add_guide
        Add a guide to the library.
    reverse_complement
        Reverse-complement all guide sequences in the library.
    to_fastq
        Write the library to a fasta file, required to build a bowtie2 index.
    to_bowtie2_index
        Build a bowtie2 index from the library.
    """
    logger = LoggerDesc()

    def __init__(self, logger=None, library_file=None):
        self.logger = logger
        self.library_file = library_file
        self.fasta_file = None
        self.bowtie2_index = None

        self.sequences = []

        self._iter_current = 0

    @classmethod
    def from_file(cls, logger=None, library_file=None):
        """Build Library instance directly from library file."""
        return cls(logger=logger, library_file=library_file)

    @property
    def library_file(self):
        """Get the library file."""
        return self._library_file

    @library_file.setter
    def library_file(self, new_library_file):
        """Set the library file. Loads and sets up the library dataframe when called.

        Parameters
        ----------
        new_library_file : str
            The path of a library file.
        """
        if not os.path.exists(new_library_file):
            raise FileNotFoundError(
                f"Library file '{new_library_file}' does not exist!"
            )
        if os.path.splitext(new_library_file)[-1] not in Config.LIBRARY_FILE_FORMATS:
            raise ValueError(
                f"Invalid library file format! Try one of the following types: "
                f"{', '.join(Config.LIBRARY_FILE_FORMATS)}"
            )
        self._library_file = new_library_file
        self.lib_df, self.sequence_length = self.read_library_file(
            library_file=self.library_file
        )
        self.sequences = list(self.lib_df["sequence"].values)
        self.name = Path(self._library_file).stem

        self._iter_max = self.lib_df.shape[0]

    def read_library_file(self, library_file=None):
        """Read the library file in one of the following formats: "*.xlsx", "*.csv", "*.tsv", "*.txt".

        Parameters
        ----------
        library_file : str
            Path to the library file in one of the following formats: "*.xlsx", "*.csv", "*.tsv", "*.txt".
        Returns
        -------
        cleaned_library_df : DataFrame
            The library as a pandas dataframe.
        sequence_length : int
            Length of the library sequences.
        """

        library_df = read_df(
            self.logger,
            df_file=library_file,
            columns=[0, 1],
            header=None,
            names=["guide", "sequence"],
        )

        library_df["sequence"] = library_df["sequence"].str.upper()
        cleaned_library_df = remove_duplicates(
            self.logger, dataframe=library_df, columns=["guide", "sequence"]
        )
        cleaned_library_df = remove_invalid_sequences(cleaned_library_df, self.logger)
        cleaned_library_df, sequence_length = check_sequence_lengths(
            cleaned_library_df, self.logger, keep="most abundant"
        )
        cleaned_library_df.sort_values(by="guide", ascending=True, inplace=True)
        cleaned_library_df.reset_index(drop=True, inplace=True)

        return cleaned_library_df, sequence_length

    def add_guide(self, guide_name, sequence):
        """Add a gRNA to the library.

        Parameters
        ----------
        guide_name : str
            Unique gRNA sequence.
        sequence : str
            Unique gRNA sequence.
        Returns
        -------
        None
        """
        self.logger.info(f"Adding '{guide_name}': '{sequence}'")
        if alphabet(sequence) == "DNA":
            if self.lib_df is None:
                lib = {"guide": [guide_name], "sequence": [sequence.upper()]}
                self.lib_df = pd.DataFrame.from_dict(lib, orient="columns")
                self.sequence_length = len(sequence)
            else:
                if len(sequence) != self.sequence_length:
                    raise ValueError(
                        f"The sequence length ({len(sequence)}) does not match the sequence length of the library ({self.sequence_length})!"
                    )
                if guide_name in self.lib_df["guide"].values:
                    raise ValueError(
                        f"Guide name '{guide_name}' already exists in library!"
                    )
                if sequence in self.lib_df["sequence"].values:
                    raise ValueError(
                        f"Sequence '{sequence}' already exists in library!"
                    )

                self.lib_df = pd.concat(
                    [
                        self.lib_df,
                        pd.Series({"guide": guide_name, "sequence": sequence.upper()})
                        .to_frame()
                        .T,
                    ],
                    ignore_index=True,
                )
                self.lib_df.sort_values(by="guide", ascending=True, inplace=True)
                my_sequences = self.sequences
                my_sequences.append(sequence.upper())
                self.sequences = my_sequences
                self._iter_max += 1
                self.logger.info("  OK!")
        else:
            raise ValueError(
                f"The provided sequence '{sequence}' is not a DNA sequence!"
            )

    def reverse_complement(self):
        """Reverse complement all sequences in the library.
        Returns
        -------
        None
        """
        self.logger.info(f"Reverse-complementing the library '{self.name}'")
        if self.lib_df is not None:
            self.lib_df["sequence"] = self.lib_df["sequence"].apply(reverse_complement)
        else:
            raise ValueError("Library does not exist!")

    def to_fasta(self, output_dir=None, name_suffix=""):
        """Write the library to a fasta file, required to build a bowtie2 index.

        Parameters
        ----------
        output_dir : str
            Directory to which the fasta file is written.
        name_suffix : str
            Add a suffix to the name to enable usage of the same library for different cassettes,
            e.g., a multiplex sample.

        Returns
        -------
        None
        """
        if os.path.isdir(output_dir):
            self.fasta_file = os.path.join(
                output_dir, f"{self.name}_{name_suffix}.fasta"
            )
            self.logger.info(f"Writing library to: '{self.fasta_file}'")
            with open(self.fasta_file, "w", encoding="utf-8") as writefile:
                for guide in self:
                    writefile.write(f">{guide['guide_name']}\n")
                    writefile.write(f"{guide['sequence']}\n")
        else:
            raise ValueError("Please provide a valid path to write the fasta file!")

    def to_bowtie2_index(self, output_dir=None):
        """Build a bowtie2 index from the library.

        Parameters
        ----------
        output_dir : str
            Directory to which the index is written.

        Returns
        -------
        str
        """
        if os.path.isdir(output_dir):
            if self.fasta_file is None:
                self.to_fasta(
                    output_dir, name_suffix=str(uuid.uuid1()).replace("-", "_")
                )

            self.bowtie2_index = os.path.join(
                output_dir, Path(self.fasta_file).stem + "_index"
            )

            self.logger.info("Creating bowtie2 index: %s", self.bowtie2_index)

            command = ["bowtie2-build", "-q", "-f", self.fasta_file, self.bowtie2_index]
            bowtie2_output = check_output(command)
            self.logger.info(bowtie2_output)
            return self.bowtie2_index
        return None

    def __iter__(self):
        """Iterator."""
        self._iter_current = 0
        return self

    def __next__(self):
        """Iterator next."""
        if self._iter_current < self._iter_max:
            self._iter_current += 1
            entry = self.lib_df.loc[self._iter_current - 1]
            return {"guide_name": entry["guide"], "sequence": entry["sequence"]}
        raise StopIteration

    def __len__(self):
        """Size of library."""
        if self.lib_df is not None:
            return self.lib_df.shape[0]
        return 0

    def __repr__(self):
        """Representation."""
        return f"{self.__class__.__name__} '{self.name}'"

    def __str__(self):
        """String representation."""
        return f"<{self.__class__.__name__} '{self.name}'>"


def check_sequence_lengths(dataframe, logger, keep="most abundant"):
    """Determine the length of the sequenes in the library.

    Parameters
    ----------
    dataframe : pd.DataFrame
        A pandas dataframe with a column 'sequence'.
    logger : logging.Logger
        A logging.Logger object.
    keep : str
        A string indicating the rule to select the length of the gRNA sequence. Can be 'shortest', 'longest', or 'most abundant'.

    Returns
    -------
    dataframe : pd.DataFrame
        The dataframe.
    selection[keep] : int
        The length of the library sequences.
    """
    dataframe["seq_length"] = dataframe["sequence"].apply(len)
    length_counter = Counter(dataframe["seq_length"])

    selection = {
        "shortest": min(length_counter.keys()),
        "longest": max(length_counter.keys()),
        "most abundant": length_counter.most_common(1)[0][0],
    }

    if keep in selection:
        if len(length_counter) != 1:
            logger.warning(
                f"Library sequences are not equal in length: {length_counter}. Will keep the sequences with the {keep} length ({selection[keep]})!"
            )
            logger.warning(f"{dataframe[dataframe['seq_length'] != selection[keep]]}")
        dataframe = dataframe[dataframe["seq_length"] == selection[keep]]
    else:
        raise ValueError(
            f"Argument {keep} is not supported! Try on of these: longest, shortest, most abundant."
        )

    dataframe = dataframe.drop("seq_length", axis=1)
    return dataframe, selection[keep]


def remove_invalid_sequences(dataframe, logger):
    """Removes sequences from a 'sequence' column in a pandas dataframe.

    Parameters
    ----------
    dataframe : pd.DataFrame
        A pandas dataframe with a column 'sequence'.
    logger : logging.Logger
        A logging.Logger object.

    Returns
    -------
    dataframe : pd.DataFrame
        The cleaned dataframe
    """
    dataframe["dna_seq"] = dataframe["sequence"].apply(alphabet)
    invalid_sequences_df = dataframe[dataframe["dna_seq"] != "DNA"]
    if invalid_sequences_df.shape[0] != 0:
        logger.warning(f"Invalid sequences found: {invalid_sequences_df.shape[0]}.")
        logger.warning(f"{invalid_sequences_df}")

    dataframe = dataframe[dataframe["dna_seq"] == "DNA"]
    dataframe = dataframe.drop("dna_seq", axis=1)
    return dataframe
