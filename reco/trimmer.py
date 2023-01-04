# pylint: disable=line-too-long

"""Trimmer class for ReCo."""

import subprocess

from collections import Counter

import pyfastx

from .attribute_descriptors import LoggerDesc
from .reco_config import Config


class CutadaptTrimmer:
    """
    Class representing a cutadapt trimmer.

    Parameters
    ----------
    logger : logging.Logger
        A logging.Logger object.
    max_n : int
        The cutadapt max_n parameter, number of tolerated N bases per read.
    error_tolerance : float
        The cutadapt -e parameter, the maximum error rate in the homology.

    Methods
    ------
    run_single
        Run cutadapt trimming for SingleSample
    run_paired
        Run cutadapt trimming for PairedSample
    finalize_single
        Finalize trimming by counting unique potential gRNA sequences and writing them to a fastq file.
    finalize_paired
        Finalize trimming by counting unique potential gRNA sequences and writing them to a fastq file.
    """
    logger = LoggerDesc()

    def __init__(self, logger=None, max_n=Config.CUTADAPT_MAX_N, error_tolerance=Config.CUTADAPT_ERROR_TOLERANCE):
        self.logger = logger

        self.standard_command = [
            "cutadapt",
            "--discard-untrimmed",
            #"--max-n", str(max_n),
            #"--nextseq-trim=20",
            "-e", str(error_tolerance)
        ]

    def _run(self, command):
        """
        Run cutadapt with the provided command.

        Parameters
        ----------
        command : [str]
            The cutadapt command as a list of strings.
        Returns
        -------
        proc_output : str
            Returns the cutadapt output as a string.
        """
        proc_output = subprocess.check_output(command)
        proc_output = proc_output.decode("utf-8")
        self.logger.info(proc_output)
        return proc_output

    def run_single(self, sequence_length=0, homology=None, output_file=None, fastq_file=None, cores=1):
        """
        Run the cutadapt trimmer for a SingleSample.

        Parameters
        ----------
        sequence_length : int
            Length of the gRNA library sequences.
        homology : str
            The adapter sequence that cutadapt uses for trimming, this is the sequence directly upstream of the gRNA sequence.
        output_file : str
            The path to the output file.
        fastq_file : str
            The path to the fastq file.
        cores : int
            The numver of cores to utilize.

        Returns
        -------
        cutadapt_output : str
            The cutadapt output.
        """
        command_ext = [
            "--minimum-length", str(sequence_length),
            "-l", str(sequence_length),
            "-g", homology,
            "-o", output_file,
            fastq_file,
            "-j", str(cores)
        ]
        single_command = self.standard_command + command_ext
        cutadapt_output = self._run(single_command)

        return cutadapt_output

    def run_paired(self, pair_filter=None, sequence_length_1=0, sequence_length_2=0, homology_1=None, homology_2=None,
                   output_file_1=None, output_file_2=None, fastq_file_1=None, fastq_file_2=None, cores=1):
        """
        Run the cutadapt trimmer for a PairedSample.

        Parameters
        ----------
        pair_filter : str
            The pair-filter parameter for cutadapt.
        sequence_length_1 : int
            Length of the sequences in gRNA library 1.
        sequence_length_2 : int
            Length of the sequences in gRNA library 2.
        homology_1 : str
            The adapter sequence that cutadapt uses for trimming of read 1, this is the sequence directly upstream of the gRNA sequence.
         homology_2 : str
            The adapter sequence that cutadapt uses for trimming of read 2, this is the sequence directly upstream of the gRNA sequence.
        output_file_1 : str
            The path to the read 1 output file.
        output_file_2 : str
            The path to the read 2 output file.
        fastq_file_1 : str
            The path to the fastq file 1.
        fastq_file_2 : str
            The path to the fastq file 2.
        cores : int
            The numver of cores to utilize.

        Returns
        -------
        cutadapt_output : str
            The cutadapt output.
        """
        command_ext = [
            "--pair-filter", pair_filter,
            "--minimum-length", f"{sequence_length_1}:{sequence_length_2}",
            "-g", homology_1,
            "-G", homology_2,
            "-o", output_file_1,
            "-p", output_file_2,
            fastq_file_1,
            fastq_file_2,
            "-j", str(cores)
        ]
        paired_command = self.standard_command + command_ext
        cutadapt_output = self._run(paired_command)

        return cutadapt_output

    def finalize_single(self, cutadapt_output_file=None, finalized_output=None, sequence_length=0):
        """
        Finalize trimming of a SingleSample.

        Parameters
        ----------
        cutadapt_output_file : str
            The path to the cutadapt output file.
        finalized_output : str
            The path to the final trimmed file.
        sequence_length : int
            The length of the library sequences.

        Returns
        -------

        """
        trimmed_sequences = Counter()
        for name, seq, qual in pyfastx.Fastq(cutadapt_output_file, build_index=False):
            trimmed_sequences[seq] += 1

        with open(finalized_output, "w", encoding="utf-8") as writefile:
            for seq in trimmed_sequences.most_common():
                sequence = seq[0]
                count = seq[1]
                if len(sequence) == sequence_length:
                    write_sequence(writefile, sequence, str(count))

    def finalize_paired(self, cutadapt_output_file_1=None, cutadapt_output_file_2=None, finalized_output_1=None, finalized_output_2=None, sequence_length_1=0, sequence_length_2=0):
        """
        Finalize trimming of a PairedSample.

        Parameters
        ----------
        cutadapt_output_file_1 : str
            The path to the cutadapt output file 1.
        cutadapt_output_file_2 : str
            The path to the cutadapt output file 1.
        finalized_output_1 : str
            The path to the final trimmed file 1.
        finalized_output_2 : str
            The path to the final trimmed file 2.
        sequence_length_1 : int
            The length of library 1 sequences.
        sequence_length_2 : int
            The length of library 2 sequences.

        Returns
        -------

        """
        combinations = Counter()

        r1_set = Counter()
        r2_set = Counter()
        # count all read1/read2 combinations
        for read_1, read_2 in zip(
            pyfastx.Fastq(cutadapt_output_file_1, build_index=False),
            pyfastx.Fastq(cutadapt_output_file_2, build_index=False)
        ):
            seq_1 = read_1[1][: sequence_length_1]
            seq_2 = read_2[1][: sequence_length_2]
            combinations[(seq_1, seq_2)] += 1
            r1_set[seq_1] += 1
            r2_set[seq_2] += 1

        self.logger.info(f"{len(combinations):n} combinations after trimming.")

        for final_trimmed_file, read_set in zip(
            [
                finalized_output_1,
                finalized_output_2,
            ],
            [r1_set, r2_set],
        ):
            with open(final_trimmed_file, "w", encoding="utf-8") as writefile:
                for sequence in read_set:
                    write_sequence(writefile, sequence, read_set[sequence])

        return combinations

    def __repr__(self):
        """Representation."""
        return f"{self.__class__.__name__}"

    def __str__(self):
        """String representation."""
        return f"<{self.__class__.__name__}>"


def write_sequence(output_file_handle, sequence, id_suffix):
    """Write fastq read block to a handle, add id_suffix to sequence id.

    Parameters
    ----------
    output_file_handle
        A file handle.
    sequence : str
        A DNA sequence.
    id_suffix
        A suffix that is added to a read if it is not an empty string.

    Returns
    -------

    """
    if id_suffix != "":
        id_suffix = f"_{id_suffix}"
    else:
        id_suffix = ""

    output_file_handle.write(f"@{sequence}{id_suffix}\n")
    output_file_handle.write(f"{sequence}\n")
    output_file_handle.write("+\n")
    output_file_handle.write(f"{'H' * len(sequence)}\n")
