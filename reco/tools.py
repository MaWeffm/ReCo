# pylint: disable=line-too-long

"""Provide tools for ReCo, mainly for sequences."""

import os
import zipfile

from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

from .reco_config import Config


def alphabet(text):
    """Determine the alphabet of a string, DNA or text.

    Parameters
    ----------
    text : str
        A string.

    Returns
    -------
    str
        'DNA' if the string contains only the letter 'A', 'C', 'G', and 'T' in upper or lower case, 'text' otherwise.
    """
    dna = {"A", "C", "G", "T"}
    text_char_set = set(Counter(text.upper()).keys())
    if all(n in dna for n in text_char_set):
        return "DNA"
    return "text"


def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence, alphabet: {A, C, G, T, N}.

    Parameters
    ----------
    seq : str
        A DNA sequence.

    Returns
    -------
    str
        The reverse complement of a DNA sequence
    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

    bases = list(seq.upper())
    bases = reversed([complement.get(base, base) for base in bases])
    bases = "".join(bases)
    return bases


def enzymes_in_seq(seq):
    """Finds restriction enzyme recognition sites in a DNA sequence and returns a string containing their names separated with '/'.

    Parameters
    ----------
    seq : str
        A DNA sequence

    Returns
    -------
    str
        All matching enzymes separated by '/'.
    """
    seq = seq.upper()
    enzymes = "/".join(
        [
            e
            for e in Config.ENZYMES
            if Config.ENZYMES[e] in seq or reverse_complement(Config.ENZYMES[e]) in seq
        ]
    )
    if enzymes == "":
        enzymes = "unknown"

    return enzymes


def read_df(
    logger, df_file=None, columns=None, header=None, names=None, index_col=None
):
    """Read dataframe from file.

    Parameters
    ----------
    logger : logging.Logger
         A logging.Logger object.
    df_file : str
    columns : [str]
    header : [int]
    names : [str]
    index_col : [int]

    Returns
    -------
    dataframe : pd.DataFrame
    """
    logger.info(f"Reading '{df_file}'")
    if not df_file:
        raise ValueError("Please provide a valid sample sheet file!")
    if not os.path.exists(df_file):
        raise FileNotFoundError(f"File '{df_file}' does not exist!")

    file_extension = os.path.splitext(df_file)[-1]
    text_file = False
    sep = ","

    if file_extension == ".xlsx":
        text_file = False
    elif file_extension == ".csv":
        sep = ","
        text_file = True
    elif file_extension in [".tsv", ".txt", ""]:
        sep = "\t"
        text_file = True

    if text_file:
        dataframe = pd.read_csv(
            df_file,
            sep=sep,
            usecols=columns,
            index_col=index_col,
            names=names,
            header=header,
        )
    else:
        dataframe = pd.read_excel(
            df_file, usecols=columns, index_col=index_col, names=names, header=header
        )

    return dataframe


def remove_duplicates(logger, dataframe=None, columns=None):
    """Remove duplicated rows and rows with duplicated values in specified columns.

    Parameters
    ----------
    logger
    dataframe
    columns

    Returns
    -------

    """
    dataframe = dataframe.dropna(how="all")
    duplicated_entries = dataframe[dataframe.duplicated(keep=False)]
    if duplicated_entries.shape[0] != 0:
        duplicated_indices = list(duplicated_entries.index.values)
        logger.warning(
            f"Duplicate {'entry' if duplicated_entries.shape[0] == 1 else 'entries'} in dataframe: {duplicated_indices} ({duplicated_entries.shape[0]:n}). Will keep the first occurrence of all duplicates!"
        )

        dataframe = dataframe.drop_duplicates(keep="first")

    for col in columns:
        duplicates_in_col = dataframe[dataframe.duplicated(subset=col, keep=False)]
        if duplicates_in_col.shape[0] != 0:
            duplicated_indices = list(duplicates_in_col.index.values)
            logger.warning(
                f"Duplicate(s) found in column '{col}': {duplicated_indices} ({duplicates_in_col.shape[0]}). Will keep the first occurrence!"
            )

            dataframe = dataframe.drop_duplicates(subset=col, keep="first")

    return dataframe


def get_file_stem(filename):
    """Determine the file stem of a file by iteratively removing all file extensions

    Parameters
    ----------
    filename : str
        A filename.

    Returns
    -------
    str
        The stem of the file.

    Raises
    ------
    IOError
        If the filename does not exist.
    """
    if os.path.exists(filename):
        filename_path = Path(filename)
        filename_stem = filename_path.name

        for suffix in filename_path.suffixes:
            filename_stem = filename_stem.rsplit(suffix)[0]

        return filename_stem
    raise IOError(f"File '{filename}' does not exist!")


def guide_in_seq(guide, sequence, homology_5_length, homology_3_length):
    """Find a gRNA sequence and its reverse complement in a longer DNA sequence and determine the 5' and 3' homology sequences.

    Parameters
    ----------
    guide : str
        A gRNA sequence.
    sequence : str
        A DNA sequence.
    homology_5_length : int
        The length of the 5' homology.
    homology_3_length
        The length of the 3' homology.

    Returns
    -------
    tuple
        A tuple of two tuples containing the forward and reverse matching homologies, respectively. None if no match was found.
    """
    # print(guide, sequence)
    sequence = str(sequence)
    guide_match_start = sequence.find(
        guide
    )  # str.find() returns -1 if string is not found
    guide_rc_match_start = sequence.find(reverse_complement(guide))

    forward_match = None
    reverse_match = None

    if guide_match_start >= 0:
        homology_5 = sequence[guide_match_start - homology_5_length : guide_match_start]
        homology_3 = sequence[
            guide_match_start
            + len(guide) : guide_match_start
            + len(guide)
            + homology_3_length
        ]
        forward_match = (homology_5, homology_3)

    elif guide_rc_match_start >= 0:
        homology_5 = sequence[
            guide_rc_match_start - homology_5_length : guide_rc_match_start
        ]
        homology_3 = sequence[
            guide_rc_match_start
            + len(guide) : guide_rc_match_start
            + len(guide)
            + homology_3_length
        ]
        reverse_match = (homology_5, homology_3)

    return forward_match, reverse_match


def count_sam_reads(sam_file):
    """Count number of alignments in a SAM file.

    Parameters
    ----------
    sam_file : str
        The path to a SAM file.

    Returns
    -------
    int
        The number of alignments in a SAM file.
    """
    sam_reads = 0
    bw2sam = pysam.AlignmentFile(sam_file, "rb")
    for read in bw2sam.fetch():
        sam_reads += 1
    bw2sam.close()

    return sam_reads


def zip_file_list(file_list=None, output_file=None):
    """
    Zip a file list.

    Parameters
    ----------
    file_list : [str]
        A list of file paths.
    output_file : str
        The path to the zip file.

    Returns
    -------

    """
    with zipfile.ZipFile(output_file, "w") as zip_me:
        for file in file_list:
            zip_me.write(
                file, os.path.basename(file), compress_type=zipfile.ZIP_DEFLATED
            )
