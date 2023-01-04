# pylint: disable=line-too-long

"""Cassette finder for ReCo."""

from collections import Counter

import random
import tqdm

import pyfastx

from .attribute_descriptors import LoggerDesc
from .reco_config import Config
from .tools import guide_in_seq


class CassetteFinder:
    """
    Class to find cassette information based on library sequences and vector information

    """
    logger = LoggerDesc()

    def __init__(self, logger=None, library=None, fastq_file=None, h5_length=Config.HOMOLOGY_5_LENGTH, h3_length=Config.HOMOLOGY_3_LENGTH):
        self.logger = logger
        self.library = library
        self.fastq_file = fastq_file

        self.h5_length = h5_length
        self.h3_length = h3_length

        self.warnings = []

    @classmethod
    def run(cls, logger=None, library=None, fastq_file=None, h5_length=Config.HOMOLOGY_5_LENGTH, h3_length=Config.HOMOLOGY_3_LENGTH,
            guides_to_test_percent=Config.GUIDES_TO_TEST_PERCENT, reads_to_test=Config.READS_TO_TEST):
        """
        Run cassette finder.

        Parameters
        ----------
        logger : logging.Logger
            A logging.Logger object.
        library : Library
            A library object.
        fastq_file : str
            The path to a fastq file.
        h5_length : int
            Length of the upstream homology.
        h3_length : int
            Length of the downstream homology.
        guides_to_test_percent: int
            The percentage of library sequences to test.
        reads_to_test : int
            The absolute number of reads to test

        Returns
        -------

        """
        c_finder = cls(logger=logger, library=library, fastq_file=fastq_file, h5_length=h5_length, h3_length=h3_length)

        return c_finder.find(guides_to_test_percent=guides_to_test_percent, reads_to_test=reads_to_test)

    def find(self, guides_to_test_percent=Config.GUIDES_TO_TEST_PERCENT, reads_to_test=Config.READS_TO_TEST):
        """
        Wrapper for find_guides_in_reads.

        Parameters
        ----------
        guides_to_test_percent: int
            The percentage of library sequences to test.
        reads_to_test : int
            The absolute number of reads to test

        Returns
        -------

        """
        rand_guides = self.get_rand_guides(guides_to_test_percent=guides_to_test_percent)
        test_reads = self.get_reads_to_test(reads_to_test=reads_to_test)

        return self.find_guides_in_reads(rand_guides, test_reads)

    def get_rand_guides(self, guides_to_test_percent):
        """
        Generate a random set of library sequences.

        Parameters
        ----------
        guides_to_test_percent : int
            The percentage of library sequences to test.

        Returns
        -------
        rand_guides:
            A list of sequences.
        """
        if guides_to_test_percent == 100:
            rand_guides = self.library.lib_df["sequence"].values
        else:
            rand_guides = random.sample(self.library.lib_df["sequence"].values, max(1, int((len(self.library) / 100) * guides_to_test_percent)))
        self.logger.info("Testing %s guides", len(rand_guides))
        return rand_guides

    def get_reads_to_test(self, reads_to_test):
        """
        Get reads to test.

        Parameters
        ----------
        reads_to_test : int
            The number of reads to test.

        Returns
        -------
        test_reads
            A list of reads.
        """
        test_reads = [""] * reads_to_test
        i = 0
        try:
            for name, seq, qual in pyfastx.Fastq(self.fastq_file, build_index=False):
                if i < Config.READS_TO_TEST:
                    test_reads[i] = seq
                else:
                    break
                i += 1
        except UnicodeDecodeError as exc:
            print("Unicode decode error:", exc)
        self.logger.info("Testing %s reads", len(test_reads))
        return test_reads

    def find_guides_in_reads(self, rand_guides, test_reads):
        """
        Determine the abundance of a set of gRNA sequences in a set of reads.

        Parameters
        ----------
        rand_guides : [str]
            A list of gRNA sequences.
        test_reads : [str]
            A list of read sequences.

        Returns
        -------
        dict
            {"lib_dir": "direction", "homology": reverse_homologies_5.most_common(1)[0][0]}
        """
        forward_matches = 0
        reverse_matches = 0

        forward_homologies_5 = Counter()
        forward_homologies_3 = Counter()

        reverse_homologies_5 = Counter()
        reverse_homologies_3 = Counter()

        for guide in tqdm.tqdm(
            rand_guides,
            desc=f"Determine library direction ({len(rand_guides):n} guides in {len(test_reads):n} reads)",
            total=len(rand_guides),
            unit=" guides"
        ):
            for read in test_reads:
                forward_match, reverse_match = guide_in_seq(
                    guide,
                    read,
                    homology_5_length=self.h5_length,
                    homology_3_length=self.h3_length,
                )
                if forward_match:
                    forward_matches += 1
                    forward_homologies_5[forward_match[0]] += 1
                    forward_homologies_3[forward_match[1]] += 1
                if reverse_match:
                    reverse_matches += 1
                    reverse_homologies_5[reverse_match[0]] += 1
                    reverse_homologies_3[reverse_match[1]] += 1

        self.logger.info(f"Forward matches: {forward_matches}")
        self.logger.info(
            f"Forward 5' homologies: {forward_homologies_5.most_common(10)}"
        )
        self.logger.info(
            f"Forward 3' homologies: {forward_homologies_3.most_common(10)}"
        )
        self.logger.info(f"Reverse matches: {reverse_matches}")
        self.logger.info(
            f"Reverse 5' homologies: {reverse_homologies_5.most_common(10)}"
        )
        self.logger.info(
            f"Reverse 3' homologies: {reverse_homologies_3.most_common(10)}"
        )

        if forward_matches == 0 and reverse_matches == 0:
            return {"lib_dir": None, "homology": None}
        elif forward_matches > reverse_matches:
            # direction is forward
            if (
                forward_homologies_5.most_common()[1][1]
                >= (forward_homologies_5.most_common()[0][1] / 100.0)
                * Config.HOMOLOGY_ABUNDANCE_DIFFERENCE_PERCENT
            ):
                self.logger.warning(
                    " => More than one forward homology found, sample mix-up or mislabeling?"
                )
                self.warnings.append(
                    "More than one forward homology found, sample mix-up or mislabeling?"
                )
                self.warnings.append(f" => {str(forward_homologies_5.most_common(5))}")
            return {"lib_dir": "forward", "homology": forward_homologies_5.most_common(1)[0][0]}
        else:
            # direction is reverse
            if (
                reverse_homologies_5.most_common()[1][1]
                >= (reverse_homologies_5.most_common()[0][1] / 100.0)
                * Config.HOMOLOGY_ABUNDANCE_DIFFERENCE_PERCENT
            ):
                self.logger.warning(
                    " => More than one reverse homology found, sample mix-up or mislabeling?"
                )
                self.warnings.append(
                    "More than one reverse homology found, sample mix-up or mislabeling?"
                )
                self.warnings.append(f" => {str(reverse_homologies_5.most_common(5))}")
            return {"lib_dir": "reverse", "homology": reverse_homologies_5.most_common(1)[0][0]}

    def __repr__(self):
        return f"<{self.__class__.__name__}>"

    def __str__(self):
        return f"{self.__class__.__name__}"
