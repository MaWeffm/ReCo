# pylint: disable=line-too-long

"""Vector class for ReCo."""

from pathlib import Path
import os

from Bio import SeqIO
from snapgene_reader import snapgene_file_to_dict

from .attribute_descriptors import LoggerDesc
from .attribute_descriptors import NameDesc
from .tools import alphabet
from .tools import enzymes_in_seq
from .tools import reverse_complement

from .reco_config import Config


class Vector:
    """Class representing a DNA vector.

    Parameters
    ----------
    name : str
        Name of the vector.
    logger : logging.Logger
        A logging.Logger object.
    vector_file : str
        Path to the library in one of the following formats:
        ".dna", ".fasta", ".fa", ".gb", ".gbk".
    Methods
    -------

    """
    logger = LoggerDesc()
    name = NameDesc()

    def __init__(self, name=None, logger=None, vector_file=None):
        self.name = name
        self.logger = logger
        self.vector_file = vector_file

        self.logger.info("Vector: %s (%s)", self.name, self.vector_file)

    @classmethod
    def from_file(cls, name=None, logger=None, vector_file=None):
        """Build vector instance from vector file."""
        return cls(name=name, logger=logger, vector_file=vector_file)

    @property
    def vector_file(self):
        """Get the vector file."""
        return self._vector_file

    @vector_file.setter
    def vector_file(self, new_vector_file):
        """Set the vector file."""

        if os.path.isfile(new_vector_file):
            if os.path.splitext(new_vector_file)[-1] not in Config.VECTOR_FILE_FORMATS:
                raise ValueError(
                    f"Invalid library file format! Try one of the following types: "
                    f"{', '.join(Config.VECTOR_FILE_FORMATS)}"
                )
            self._vector_file = new_vector_file

            vector_file_path_suffix = Path(self._vector_file).suffixes[-1]
            if vector_file_path_suffix == ".dna":
                vector = snapgene_file_to_dict(self._vector_file)
                if vector["isDNA"]:
                    self.sequence = vector["seq"].upper()
                else:
                    raise ValueError(f"The sequence in '{self._vector_file}' is not a DNA sequence!")
            elif vector_file_path_suffix in [".fasta", ".fa"]:
                self.sequence = str(SeqIO.read(self._vector_file, "fasta").seq).upper()
            elif vector_file_path_suffix in [".gb", ".gbk"]:
                self.sequence = str(SeqIO.read(self._vector_file, "genbank").seq).upper()
        elif alphabet(new_vector_file) == "DNA":
            self.sequence = new_vector_file
        else:
            raise ValueError(f"This is nor a file and neither a DNA sequence: '{new_vector_file}'")

    def find_template(self, template_length, prefix=None):
        """Find template sequence and potential restriction enzyme recognition sites for a given prefix."""
        if prefix and isinstance(prefix, str):
            prefix = prefix.upper()
            if alphabet(prefix) == "DNA":
                # template_sequence = None
                prefix_rc = reverse_complement(prefix)

                prefix_index = self.sequence.find(prefix)
                prefix_rc_index = self.sequence.find(prefix_rc)

                if prefix_index > -1:
                    template_sequence = self.sequence[prefix_index + len(prefix): prefix_index + len(prefix) + template_length]
                elif prefix_rc_index > -1:
                    template_sequence = reverse_complement(self.sequence[prefix_rc_index - template_length: prefix_rc_index])
                else:
                    raise ValueError(f"Could not find the prefix '{prefix}' or its reverse complement '{prefix_rc}'")
                template_name = enzymes_in_seq(template_sequence)
                return {"template_name": template_name, "template_sequence": template_sequence}

            raise ValueError("Prefix must be DNA string!")
        raise ValueError("Prefix is required as a string!")

    def __len__(self):
        """Length of vector sequence."""
        if self.sequence:
            return len(self.sequence)
        return 0

    def __repr__(self):
        """Representation."""
        return f"{self.__class__.__name__} '{self.name}'"

    def __str__(self):
        """String representation."""
        return f"<{self.__class__.__name__} '{self.name}'>"
