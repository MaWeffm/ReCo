"""Provide ReCo configuration parameters."""


class Config:  # pylint: disable=too-few-public-methods
    """The ReCo configuration class."""

    LIBRARY_FILE_FORMATS = [".xlsx", ".csv", ".tsv", ".txt", ""]
    SAMPLE_SHEET_FILE_FORMATS = [".xlsx", ".csv", ".tsv", ".txt", ""]
    VECTOR_FILE_FORMATS = [".dna", ".fasta", ".fa", ".gb", ".gbk"]
    FASTQ_FILE_FORMATS = [".fastq"]

    ENZYMES = {
        "PacI": "TTAATTAA",
        "I-SceI": "TAGGGATAACAGGGTAAT",
        "AgeI": "ACCGGT",
        "I-CeuI": "TAACGGTCCTAAGGTAGCGA",
        "EcoRI": "GAATTC",
        "BamHI": "GGATCC",
        "SnaBI": "TACGTA",
    }

    BOWTIE2_SEED_LENGTH = 12
    BOWTIE2_SEED_MISMATCH = 1
    BOWTIE2_INTERVAL_FUNCTION = "S,1.0,0.75"

    CUTADAPT_PAIR_FILTER = "any"
    CUTADAPT_MAX_N = 2
    CUTADAPT_ERROR_TOLERANCE = 0.1
    CUTADAPT_SEQUENCE_LENGTH = 20

    GUIDES_TO_TEST_PERCENT = 100
    READS_TO_TEST = 500

    HOMOLOGY_5_LENGTH = 12
    HOMOLOGY_3_LENGTH = 12

    HOMOLOGY_ABUNDANCE_DIFFERENCE_PERCENT = 20

    TOP_FAILED_SEQUENCES = 100
