# pylint: disable=line-too-long

"""Sample classes for ReCo."""

import datetime
import os.path
import shutil

from collections import Counter
from pathlib import Path

import numpy
import pandas as pd

from .aligner import Bowtie2Aligner
from .attribute_descriptors import FileDesc
from .attribute_descriptors import NameDesc
from .cassette_finder import CassetteFinder
from .library import Library
from .logger import setup_logger
from .plot import PlotPanel
from .readcount_writer import ReadcountWriter
from .reco_config import Config
from .tools import zip_file_list
from .trimmer import CutadaptTrimmer
from .vector import Vector

from ._version import get_versions


class Sample:
    """Sample parent class, should be treated as abstract.

    Parameters
    ----------
    name: str
        The name of the sample.
    vector_file : str
        Path to the library in one of the following formats:
        ".dna", ".fasta", ".fa", ".gb", ".gbk".
    expected_reads : int
        The number of expected reads for the sample.

    """
    name = NameDesc()
    vector_file = FileDesc()

    def __init__(self, name=None, vector_file=None, expected_reads=0, output_dir=None):
        self.name = name
        self.vector_file = vector_file

        self.expected_reads = expected_reads

        self.output_dir = None
        self.bowtie2_index_dir = None

        self.logfile = None

        self.output_dir = os.path.join(output_dir, self.name)
        self.bowtie2_index_dir = os.path.join(self.output_dir, "bowtie2_index")
        self.create_dirs()
        self.logfile = os.path.join(self.output_dir, f"ReCo_{self.name}.log")

        self.logger = setup_logger(self.name, self.logfile, stdout=True)
        self.logger.info("ReCo v%s: initializing sample '%s'", get_versions()['version'], self.name)
        self.logger.info("Created: %s", self.output_dir)
        self.logger.info("Created: %s", self.bowtie2_index_dir)

        self.vector = None
        if self.vector_file:
            self.vector = Vector.from_file(name=Path(self.vector_file).stem, logger=self.logger, vector_file=self.vector_file)

        self.observed_reads = 0
        self.trimmed_reads = 0

        self.rc_table = None

        self.failed_file = os.path.join(self.output_dir, f"{self.name}_failed.csv")
        self.top_failed_file = os.path.join(self.output_dir, f"{self.name}_top{Config.TOP_FAILED_SEQUENCES}_failed.csv")

        self.report_file = os.path.join(self.output_dir, "report.txt")

        self.zip_file = os.path.join(self.output_dir, f"{self.name}_results.zip")

        self.warnings = []

    @property
    def expected_reads(self):
        """Get number of expected reads for this sample."""
        return self._expected_reads

    @expected_reads.setter
    def expected_reads(self, new_expected_reads):
        """Set number of expected reads for this sample."""
        if new_expected_reads == "" or isinstance(new_expected_reads, (int, numpy.int64)):
            self._expected_reads = new_expected_reads
        else:
            raise ValueError("Expected reads must be an integer!")

    def create_dirs(self):
        """Create sample directories."""
        try:
            os.makedirs(self.bowtie2_index_dir, exist_ok=False)
        except FileExistsError as exc:
            raise FileExistsError(f"Output dir '{self.output_dir}' already exists!") from exc

    def get_trimming_rate(self, cutadapt_output=None, sample_type=None):
        """
        Extract the number of processed and written reads from cutadapt output

        Parameters
        ----------
        cutadapt_output : str
            Cutadapt output as a string
        sample_type : str
            The sample type.

        Returns
        -------

        """
        written_text = ""
        processed_text = ""

        if sample_type == "single":
            written_text = "Reads written (passing filters):"
            processed_text = "Total reads processed:"
        elif sample_type == "paired":
            written_text = "Pairs written (passing filters):"
            processed_text = "Total read pairs processed:"

        for line in cutadapt_output.split("\n"):
            if written_text in line:
                self.trimmed_reads = int(
                    line.split(":")[1].strip().split(" ")[0].replace(",", "")
                )
            if processed_text in line:
                self.observed_reads = int(
                    line.split(":")[1].strip().split(" ")[0].replace(",", "")
                )
        self.logger.info("Processed reads: %s", f"{self.observed_reads:n}")
        self.logger.info("Trimmed reads: %s", f"{self.trimmed_reads:n}")

    def __repr__(self):
        """Representation."""
        return f"{self.__class__.__name__} '{self.name}'"

    def __str__(self):
        """String representation."""
        return f"<{self.__class__.__name__} '{self.name}'>"


class SingleSample(Sample):
    """Class representing a single sample: a single read, a single gRNA."""
    def __init__(self, sample_dict=None, output_dir=None):
        super().__init__(name=sample_dict["Sample name"], vector_file=sample_dict["Vector"], expected_reads=sample_dict["Expected reads"], output_dir=output_dir)

        sd_fastq_1 = sample_dict["FastQ read 1"]
        sd_fastq_2 = sample_dict["FastQ read 2"]

        if sd_fastq_1:
            self.fastq_file = sd_fastq_1
        elif sd_fastq_2:
            self.fastq_file = sd_fastq_2
        self.cutadapt_output_file = os.path.join(self.output_dir, Path(self.fastq_file).stem + ".cutadapt_out")
        self.final_trimmed_file = os.path.join(Path(self.cutadapt_output_file).parent, Path(self.cutadapt_output_file).stem + ".final_cutadapt_out")
        self.sam_file = os.path.join(Path(self.final_trimmed_file).parent, Path(self.final_trimmed_file).stem + ".sam")
        self.logger.info("FastQ: %s", self.fastq_file)

        sd_lib_1 = sample_dict["Lib 1"]
        sd_lib_2 = sample_dict["Lib 2"]

        if sd_lib_1:  # if field is not empty use this
            self.lib_file = sd_lib_1
        elif sd_lib_2:  # otherwise use this
            self.lib_file = sd_lib_2
        self.lib = Library(logger=self.logger, library_file=self.lib_file)
        self.logger.info("Library size: %s", len(self.lib))

        self.cassette_information = None
        self.template_information = None

        self.bowtie2_index = None

        self.readcount_file = os.path.join(self.output_dir, f"{self.name}_final_guidecounts.csv")

    def info(self):
        """Return fastq and library information."""
        return {"fastq_file": self.fastq_file, "lib_file": self.lib_file}

    def run(self, remove_unused_files=True, cores=1):
        """
        Run single sample

        Parameters
        ----------
        remove_unused_files : bool
            Whether to remove the intermediate files.
        cores : int
            Number of cores to utilize

        Returns
        -------

        """
        print(f"Running {self.name}")
        cassette_information_time = datetime.datetime.now()
        yield "Determine cassette information ..."
        self.cassette_information = CassetteFinder.run(logger=self.logger, library=self.lib, fastq_file=self.fastq_file, h5_length=Config.HOMOLOGY_5_LENGTH, h3_length=Config.HOMOLOGY_3_LENGTH, guides_to_test_percent=Config.GUIDES_TO_TEST_PERCENT, reads_to_test=Config.READS_TO_TEST)
        if self.cassette_information["lib_dir"] is None and self.cassette_information["homology"] is None:
            self.logger.error("No gRNAs found! Are you using the correct library?")
            self.warnings.append("No gRNAs found! Are you using the correct library?")
            return "No gRNAs found! Are you using the correct library?"
        if self.cassette_information["lib_dir"] == "reverse":
            self.lib.reverse_complement()
        print("  Cassette determination in:", datetime.datetime.now() - cassette_information_time)

        template_information_time = datetime.datetime.now()
        yield "Determine 3Cs template sequence ..."
        if self.vector:
            self.template_information = self.vector.find_template(template_length=self.lib.sequence_length, prefix=self.cassette_information["homology"])
            self.lib.add_guide(guide_name=f"template-{self.template_information['template_name']}", sequence=self.template_information["template_sequence"])
            self.logger.info("Template: %s", self.template_information)
        start_time = datetime.datetime.now()
        yield "Build Bowtie2 index ..."
        self.bowtie2_index = self.lib.to_bowtie2_index(self.bowtie2_index_dir)
        print("  Template determination in:", datetime.datetime.now() - template_information_time)

        trimming_time = datetime.datetime.now()
        yield "Trim ..."
        sample_trimmer = CutadaptTrimmer(logger=self.logger, max_n=Config.CUTADAPT_MAX_N, error_tolerance=Config.CUTADAPT_ERROR_TOLERANCE)

        cutadapt_output = sample_trimmer.run_single(sequence_length=self.lib.sequence_length, homology=self.cassette_information["homology"],
                                                    output_file=self.cutadapt_output_file, fastq_file=self.fastq_file, cores=cores)
        self.get_trimming_rate(cutadapt_output=cutadapt_output, sample_type="single")
        sample_trimmer.finalize_single(cutadapt_output_file=self.cutadapt_output_file, finalized_output=self.final_trimmed_file, sequence_length=self.lib.sequence_length)
        print("  Trimming in:", datetime.datetime.now() - trimming_time)

        align_time = datetime.datetime.now()
        yield "Align ..."
        aligner = Bowtie2Aligner(logger=self.logger, seed_length=Config.BOWTIE2_SEED_LENGTH, seed_mismatch=Config.BOWTIE2_SEED_MISMATCH, interval_function=Config.BOWTIE2_INTERVAL_FUNCTION,
                                 bowtie2_index=self.bowtie2_index, final_trimmed_file=self.final_trimmed_file, sam_file=self.sam_file)
        aligner.run(cores=cores)
        print("  Alignment in:", datetime.datetime.now() - align_time)

        collect_time = datetime.datetime.now()
        yield "Collect read counts ..."
        rc_writer = ReadcountWriter(logger=self.logger)
        count_tables = rc_writer.counts_from_sam(sam_file=self.sam_file, library=self.lib)
        self.rc_table = rc_writer.write_counter(count_tables=count_tables, count_col_name=self.name, counts_file=self.readcount_file, failed_file=self.failed_file, top_failed_file=self.top_failed_file)
        # print(count_tables["guide_hitlist"].most_common(25))
        end_time = datetime.datetime.now()
        print("  Collected read counts in:", datetime.datetime.now() - collect_time)

        print(" ===========> Running time: ", end_time - start_time)

        yield "Generate plots ..."
        p_panel = PlotPanel(logger=self.logger, sample_name=self.name, sample_type="single", read_count_df=self.rc_table,
                            expected_reads=self.expected_reads, observed_reads=self.observed_reads, trimmed_reads=self.trimmed_reads,
                            aligned_reads=self.rc_table[self.name].sum(), template_marker="template", output_dir=self.output_dir)
        p_panel.generate()

        yield "Write report ..."
        self.write_report()

        yield "Build results package ..."
        zip_file_list(file_list=[self.readcount_file, self.top_failed_file, p_panel.pdf_file, p_panel.png_file, self.report_file],
                      output_file=self.zip_file)

        yield "Clean-up ..."
        if remove_unused_files:
            print("Clean up!")
            self.clean_up()
        else:
            print("No clean up, keep intermediate files!")

        self.get_trimming_rate(cutadapt_output=cutadapt_output, sample_type="single")
        self.logger.info("Aligned reads: %s", f"{self.rc_table[self.name].sum():n}")

        return "Done!"

    def clean_up(self):
        """
        Remove intermediate files.

        Returns
        -------

        """
        self.logger.info("Cleaning up!")
        os.remove(self.cutadapt_output_file)
        os.remove(self.final_trimmed_file)
        os.remove(self.sam_file)
        shutil.rmtree(self.bowtie2_index_dir)

    def write_report(self):
        """
        Write informative report.

        Returns
        -------

        """
        with open(self.report_file, "w", encoding="utf-8") as writefile:
            writefile.write(f"ReCo {get_versions()['version']}\n")
            writefile.write(f"{self}\n")
            writefile.write(f"Library: {self.lib_file} ({len(self.lib)} guides)\n")
            writefile.write(f"Vector: {self.vector_file}\n")
            writefile.write(f"FastQ: {self.fastq_file}\n")
            writefile.write(f"Homology: {self.cassette_information['homology']}\n")
            writefile.write(f"Library direction: {self.cassette_information['lib_dir']}\n")
            if self.template_information:
                writefile.write(f"Template name: {self.template_information['template_name']}\n")
                writefile.write(f"Template sequence: {self.template_information['template_sequence']}\n")
            else:
                writefile.write("No template vector provided.")
            writefile.write(f"Expected reads: {self.expected_reads:n}\n")
            writefile.write(f"Observed reads: {self.observed_reads:n} ({round(self.observed_reads / (self.expected_reads / 100), 2)}% of expected)\n")
            writefile.write(f"Trimmed reads: {self.trimmed_reads:n} ({round(self.trimmed_reads / (self.observed_reads / 100), 2)}% of observed)\n")
            writefile.write(f"Aligned reads: {self.rc_table[self.name].shape[0]:n} ({round(self.rc_table[self.name].shape[0] / (self.observed_reads / 100), 2)}% of observed)\n")


class PairedSample(Sample):
    """Class representing a paired sample: paired-end reads, one gRNA in each read."""

    def __init__(self, sample_dict=None, output_dir=None):
        super().__init__(name=sample_dict["Sample name"], vector_file=sample_dict["Vector"], expected_reads=sample_dict["Expected reads"], output_dir=output_dir)

        self.fastq_file_1 = sample_dict["FastQ read 1"]
        self.fastq_file_2 = sample_dict["FastQ read 2"]

        self.cutadapt_output_file_1 = os.path.join(self.output_dir, Path(self.fastq_file_1).stem + ".cutadapt_out")
        self.final_trimmed_file_1 = os.path.join(Path(self.cutadapt_output_file_1).parent, Path(self.cutadapt_output_file_1).stem + ".final_cutadapt_out")
        self.sam_file_1 = os.path.join(Path(self.final_trimmed_file_1).parent, Path(self.final_trimmed_file_1).stem + ".sam")

        self.cutadapt_output_file_2 = os.path.join(self.output_dir, Path(self.fastq_file_2).stem + ".cutadapt_out")
        self.final_trimmed_file_2 = os.path.join(Path(self.cutadapt_output_file_2).parent, Path(self.cutadapt_output_file_2).stem + ".final_cutadapt_out")
        self.sam_file_2 = os.path.join(Path(self.final_trimmed_file_2).parent, Path(self.final_trimmed_file_2).stem + ".sam")

        self.logger.info("FastQ 1: %s", self.fastq_file_1)
        self.logger.info("FastQ 2: %s", self.fastq_file_2)

        self.lib_file_1 = sample_dict["Lib 1"]
        self.lib_file_2 = sample_dict["Lib 2"]

        self.lib_1 = Library(logger=self.logger, library_file=self.lib_file_1)
        self.lib_2 = Library(logger=self.logger, library_file=self.lib_file_2)

        self.logger.info("Size library 1: %s", len(self.lib_1))
        self.logger.info("Size library 2: %s", len(self.lib_2))

        self.lib_1_bowtie2_index = None
        self.lib_2_bowtie2_index = None

        self.readcount_file = os.path.join(self.output_dir, f"{self.name}_combinations_final_guidecounts.csv")

        self.readcount_file_read_1 = os.path.join(self.output_dir, f"{self.name}_R1_final_guidecounts.csv")
        self.failed_file_read_1 = os.path.join(self.output_dir, f"{self.name}_R1_failed.csv")
        self.top_failed_file_read_1 = os.path.join(self.output_dir, f"{self.name}_top{Config.TOP_FAILED_SEQUENCES}_R1_failed.csv")

        self.readcount_file_read_2 = os.path.join(self.output_dir, f"{self.name}_R2_final_guidecounts.csv")
        self.failed_file_read_2 = os.path.join(self.output_dir, f"{self.name}_R2_failed.csv")
        self.top_failed_file_read_2 = os.path.join(self.output_dir, f"{self.name}_top{Config.TOP_FAILED_SEQUENCES}_R2_failed.csv")

        self.cassette_information_1 = None
        self.cassette_information_2 = None

        self.template_information_1 = None
        self.template_information_2 = None

        self.rc_table_1 = None
        self.rc_table_2 = None

        self.aligned_combis = 0

    def info(self):
        """Return fastq and library information."""
        return {"fastq_file_1": self.fastq_file_1, "fastq_file_2": self.fastq_file_2, "lib_file_1": self.lib_file_1, "lib_file_2": self.lib_file_2}

    def run(self, remove_unused_files=True, cores=1):
        """
        Run single sample

        Parameters
        ----------
        remove_unused_files : bool
            Whether to remove the intermediate files.
        cores : int
            Number of cores to utilize

        Returns
        -------

        """
        yield "Determine cassette information ..."
        self.cassette_information_1 = CassetteFinder.run(logger=self.logger, library=self.lib_1, fastq_file=self.fastq_file_1, h5_length=Config.HOMOLOGY_5_LENGTH, h3_length=Config.HOMOLOGY_3_LENGTH, guides_to_test_percent=Config.GUIDES_TO_TEST_PERCENT, reads_to_test=Config.READS_TO_TEST)
        self.cassette_information_2 = CassetteFinder.run(logger=self.logger, library=self.lib_2, fastq_file=self.fastq_file_2, h5_length=Config.HOMOLOGY_5_LENGTH, h3_length=Config.HOMOLOGY_3_LENGTH, guides_to_test_percent=Config.GUIDES_TO_TEST_PERCENT, reads_to_test=Config.READS_TO_TEST)
        if self.cassette_information_1["lib_dir"] is None and self.cassette_information_1["homology"] is None:
            self.logger.error("No gRNAs found in R1! Are you using the correct library?")
            self.warnings.append("No gRNAs found in R1! Are you using the correct library?")
            return "No gRNAs found in R1! Are you using the correct library?"
        if self.cassette_information_2["lib_dir"] is None and self.cassette_information_2["homology"] is None:
            self.logger.error("No gRNAs found in R2! Are you using the correct library?")
            self.warnings.append("No gRNAs found in R2! Are you using the correct library?")
            return "No gRNAs found in R2! Are you using the correct library?"
        if self.cassette_information_1["lib_dir"] == "reverse":
            self.lib_1.reverse_complement()
        if self.cassette_information_2["lib_dir"] == "reverse":
            self.lib_2.reverse_complement()

        yield "Determine 3Cs template sequences ..."
        if self.vector:
            self.template_information_1 = self.vector.find_template(template_length=self.lib_1.sequence_length, prefix=self.cassette_information_1["homology"])
            self.logger.info("Template 1: %s", self.template_information_1)
            self.lib_1.add_guide(guide_name=f"template-{self.template_information_1['template_name']}", sequence=self.template_information_1["template_sequence"])

            self.template_information_2 = self.vector.find_template(template_length=self.lib_2.sequence_length, prefix=self.cassette_information_2["homology"])
            self.logger.info("Template 2: %s", self.template_information_2)
            self.lib_2.add_guide(guide_name=f"template-{self.template_information_2['template_name']}", sequence=self.template_information_2["template_sequence"])

        start_time = datetime.datetime.now()
        yield "Build Bowtie2 indices ..."
        self.lib_1_bowtie2_index = self.lib_1.to_bowtie2_index(self.bowtie2_index_dir)
        self.lib_2_bowtie2_index = self.lib_2.to_bowtie2_index(self.bowtie2_index_dir)

        yield "Trim ..."
        sample_trimmer = CutadaptTrimmer(logger=self.logger, max_n=Config.CUTADAPT_MAX_N, error_tolerance=Config.CUTADAPT_ERROR_TOLERANCE)
        cutadapt_output_file_1 = os.path.join(self.output_dir, Path(self.fastq_file_1).stem + ".cutadapt_out")
        cutadapt_output_file_2 = os.path.join(self.output_dir, Path(self.fastq_file_2).stem + ".cutadapt_out")
        cutadapt_output = sample_trimmer.run_paired(pair_filter=Config.CUTADAPT_PAIR_FILTER, sequence_length_1=self.lib_1.sequence_length, sequence_length_2=self.lib_2.sequence_length,
                                                    homology_1=self.cassette_information_1["homology"], homology_2=self.cassette_information_2["homology"],
                                                    output_file_1=cutadapt_output_file_1, output_file_2=cutadapt_output_file_2,
                                                    fastq_file_1=self.fastq_file_1, fastq_file_2=self.fastq_file_2, cores=cores)
        self.get_trimming_rate(cutadapt_output=cutadapt_output, sample_type="paired")
        combinations = sample_trimmer.finalize_paired(cutadapt_output_file_1=self.cutadapt_output_file_1, cutadapt_output_file_2=self.cutadapt_output_file_2,
                                                      finalized_output_1=self.final_trimmed_file_1, finalized_output_2=self.final_trimmed_file_2,
                                                      sequence_length_1=self.lib_1.sequence_length, sequence_length_2=self.lib_2.sequence_length)

        yield "Align ..."
        aligner_1 = Bowtie2Aligner(logger=self.logger, seed_length=Config.BOWTIE2_SEED_LENGTH, seed_mismatch=Config.BOWTIE2_SEED_MISMATCH, interval_function=Config.BOWTIE2_INTERVAL_FUNCTION,
                                   bowtie2_index=self.lib_1_bowtie2_index, final_trimmed_file=self.final_trimmed_file_1, sam_file=self.sam_file_1)
        aligner_1.run(cores=cores)

        aligner_2 = Bowtie2Aligner(logger=self.logger, seed_length=Config.BOWTIE2_SEED_LENGTH, seed_mismatch=Config.BOWTIE2_SEED_MISMATCH, interval_function=Config.BOWTIE2_INTERVAL_FUNCTION,
                                   bowtie2_index=self.lib_2_bowtie2_index, final_trimmed_file=self.final_trimmed_file_2, sam_file=self.sam_file_2)
        aligner_2.run(cores=cores)

        yield "Collect read counts ..."
        rc_writer = ReadcountWriter(logger=self.logger)
        count_tables_1 = rc_writer.counts_from_sam(sam_file=self.sam_file_1, library=self.lib_1)
        count_tables_2 = rc_writer.counts_from_sam(sam_file=self.sam_file_2, library=self.lib_2)

        self.rc_table_1 = rc_writer.write_counter(count_tables=count_tables_1, count_col_name=self.name + "_R1", counts_file=self.readcount_file_read_1, failed_file=self.failed_file_read_1, top_failed_file=self.top_failed_file_read_1)
        self.rc_table_2 = rc_writer.write_counter(count_tables=count_tables_2, count_col_name=self.name + "_R2", counts_file=self.readcount_file_read_2, failed_file=self.failed_file_read_2, top_failed_file=self.top_failed_file_read_2)

        self.aligned_combis = self.write_combinations(combinations, count_tables_1["sequence_hitlist"], count_tables_2["sequence_hitlist"])

        read_count_df = pd.read_csv(self.readcount_file, index_col=["Guide_1", "Guide_2"])
        end_time = datetime.datetime.now()

        print(" ===========> Running time: ", end_time - start_time)
        yield "Generate plots ..."
        p_panel_combis = PlotPanel(logger=self.logger, sample_name=self.name, sample_type="paired", read_count_df=read_count_df,
                                   expected_reads=self.expected_reads, observed_reads=self.observed_reads, trimmed_reads=self.trimmed_reads,
                                   aligned_reads=read_count_df[self.name].sum(), template_marker="template", output_dir=self.output_dir)
        p_panel_combis.generate()

        p_panel_read_1 = PlotPanel(logger=self.logger, sample_name=self.name + "_R1", sample_type="single", read_count_df=self.rc_table_1,
                                   expected_reads=self.expected_reads, observed_reads=self.observed_reads, trimmed_reads=self.trimmed_reads,
                                   aligned_reads=self.rc_table_1[self.name + "_R1"].sum(), template_marker="template", output_dir=self.output_dir)
        p_panel_read_1.generate()

        p_panel_read_2 = PlotPanel(logger=self.logger, sample_name=self.name + "_R2", sample_type="single", read_count_df=self.rc_table_2,
                                   expected_reads=self.expected_reads, observed_reads=self.observed_reads, trimmed_reads=self.trimmed_reads,
                                   aligned_reads=self.rc_table_2[self.name + "_R2"].sum(), template_marker="template", output_dir=self.output_dir)
        p_panel_read_2.generate()

        self.get_trimming_rate(cutadapt_output=cutadapt_output, sample_type="paired")
        self.logger.info("Aligned reads: %s", f"{self.aligned_combis:n}")

        yield "Write report ..."
        self.write_report()

        yield "Build results package ..."
        zip_file_list(file_list=[
            self.readcount_file,
            self.readcount_file_read_1,
            self.readcount_file_read_2,
            self.top_failed_file,
            self.top_failed_file_read_1,
            self.top_failed_file_read_2,
            p_panel_combis.pdf_file, p_panel_combis.png_file,
            p_panel_read_1.pdf_file, p_panel_read_1.png_file,
            p_panel_read_2.pdf_file, p_panel_read_2.png_file,
            self.report_file],
            output_file=self.zip_file
        )

        yield "Clean-up ..."
        if remove_unused_files:
            print("Remove intermediate files!")
            self.clean_up()
        else:
            print("Keep intermediate files!")

        return "Done!"

    def write_combinations(self, combinations, count_table_1, count_table_2):
        """

        Parameters
        ----------
        combinations
        count_table_1
        count_table_2

        Returns
        -------

        """
        self.logger.info("Writing read counts to: %s", self.readcount_file)
        self.logger.info("Writing failed alignments to: %s", self.failed_file)

        lib_counts = Counter()
        failed_counts = Counter()

        for guide_1 in self.lib_1:
            for guide_2 in self.lib_2:
                lib_counts[(guide_1["guide_name"], guide_2["guide_name"])] = 0

        for combi, count in combinations.most_common():
            guide_1 = count_table_1.get(combi[0], "-")
            guide_2 = count_table_2.get(combi[1], "-")

            if (guide_1, guide_2) in lib_counts:
                lib_counts[(guide_1, guide_2)] = count
            else:
                failed_counts[(combi[0], combi[1], guide_1, guide_2)] += count

        aligned = 0
        with open(self.readcount_file, "w", encoding="utf-8") as writefile:
            writefile.write(f"Guide_1,Guide_2,{self.name}\n")
            for combi, count in lib_counts.most_common():
                writefile.write(f"{combi[0]},{combi[1]},{count}\n")
                aligned += count

        with open(self.failed_file, "w", encoding="utf-8") as failed_file_writer, open(self.top_failed_file, "w", encoding="utf-8") as top_failed_file_writer:
            failed_file_writer.write("Sequence_1,Sequence_2,Guide_1,Guide_2,count\n")
            top_failed_file_writer.write("Sequence_1,Sequence_2,Guide_1,Guide_2,count\n")
            combi_counter = 0
            for combi, count in failed_counts.most_common():
                failed_file_writer.write(f"{combi[0]},{combi[1]},{combi[2]},{combi[3]},{count}\n")

                if combi_counter < Config.TOP_FAILED_SEQUENCES:
                    top_failed_file_writer.write(f"{combi[0]},{combi[1]},{combi[2]},{combi[3]},{count}\n")
                combi_counter += 1

        return aligned

    def clean_up(self):
        """
        Remove intermediate files.

        Returns
        -------

        """
        self.logger.info("Cleaning up!")
        os.remove(self.cutadapt_output_file_1)
        os.remove(self.cutadapt_output_file_2)
        os.remove(self.final_trimmed_file_1)
        os.remove(self.final_trimmed_file_2)
        os.remove(self.sam_file_1)
        os.remove(self.sam_file_2)
        shutil.rmtree(self.bowtie2_index_dir)

    def write_report(self):
        """
        Write informative report.

        Returns
        -------

        """
        with open(self.report_file, "w", encoding="utf-8") as writefile:
            writefile.write(f"ReCo {get_versions()['version']}\n")
            writefile.write(f"{self}\n")
            writefile.write(f"Library 1: {self.lib_file_1} ({len(self.lib_1)} guides)\n")
            writefile.write(f"Library 2: {self.lib_file_2} ({len(self.lib_2)} guides)\n")
            writefile.write(f"Vector: {self.vector_file}\n")
            writefile.write(f"FastQ read 1: {self.fastq_file_1}\n")
            writefile.write(f"FastQ read 2: {self.fastq_file_2}\n")
            writefile.write(f"Homology read 1: {self.cassette_information_1['homology']}\n")
            writefile.write(f"Homology read 2: {self.cassette_information_2['homology']}\n")
            writefile.write(f"Library direction 1: {self.cassette_information_1['lib_dir']}\n")
            writefile.write(f"Library direction 2: {self.cassette_information_2['lib_dir']}\n")
            if self.template_information_1 and self.template_information_2:
                writefile.write(f"Template name 1: {self.template_information_1['template_name']}\n")
                writefile.write(f"Template name 2: {self.template_information_2['template_name']}\n")
                writefile.write(f"Template sequence read 1: {self.template_information_1['template_sequence']}\n")
                writefile.write(f"Template sequence read 2: {self.template_information_2['template_sequence']}\n")
            else:
                writefile.write("No template vector provided.")
            writefile.write(f"Expected combinations: {self.expected_reads:n}\n")
            writefile.write(f"Observed combinations: {self.observed_reads:n} ({round(self.observed_reads / (self.expected_reads / 100), 2)}% of expected)\n")
            writefile.write(f"Trimmed combinations: {self.trimmed_reads:n} ({round(self.trimmed_reads / (self.observed_reads / 100), 2)}% of observed)\n")
            writefile.write(f"Aligned combinations: {self.aligned_combis:n} ({round(self.aligned_combis / (self.observed_reads / 100), 2)}% of observed)\n")
