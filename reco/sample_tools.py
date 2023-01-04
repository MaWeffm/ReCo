# pylint: disable=line-too-long

"""Sample validator class for ReCo."""

import os

import numpy


class SampleValidator:
    """
    Class representing a sample validator.

    Methods
    -------
    validate_single_sample
        Validate a sample of type 'single'.
    validate_paired_sample
        Validate a sample of type 'paired'.
    """
    def __init__(self):
        self.result = False
        self.valid = False

    @classmethod
    def from_dict(cls, sample_dict):
        """
        Directly run the validation.

        Parameters
        ----------
        sample_dict : dict
            A dictionary containing sample information

        Returns
        -------
        dict
            Dictionary: {"valid": bool, "message": str}
        """
        if sample_dict["Sample type"] == "single":
            return validate_single_sample(sample_dict)
        if sample_dict["Sample type"] == "paired":
            return validate_paired_sample(sample_dict)
        return {"valid": False, "message": "Unknown sample type!"}

    def __repr__(self):
        return f"<{self.__class__.__name__}>"

    def __str__(self):
        return f"{self.__class__.__name__}"


def validate_single_sample(sample_dict=None):
    """
    Validate a single sample.

    Parameters
    ----------
    sample_dict : dict
        Dictionary containing sample information

    Returns
    -------
    dict
        Dictionary: {"valid": bool, "message": str}
    """
    log_text = ""
    valid = True

    if sample_dict["Vector"]:
        if not os.path.exists(sample_dict["Vector"]):
            log_text += " Vector file does not exist!"
            valid = False
    if not (sample_dict["FastQ read 1"] or sample_dict["FastQ read 2"]):
        log_text += " One fastq file is required!"
        valid = False
    fastq_file = None
    if sample_dict["FastQ read 1"]:
        fastq_file = sample_dict["FastQ read 1"]
    elif sample_dict["FastQ read 2"]:
        fastq_file = sample_dict["FastQ read 2"]
    if not os.path.exists(fastq_file):
        log_text += f" FastQ file {fastq_file} does not exist!"
        valid = False
    if not (sample_dict["Lib 1"] or sample_dict["Lib 2"]):
        log_text += " One library file is required!"
        valid = False
    if not (sample_dict["Expected reads"] == 0 or isinstance(sample_dict["Expected reads"], (int, numpy.int64))):
        log_text += " Expected reads must be an integer!"
        valid = False

    return {"valid": valid, "message": log_text}


def validate_paired_sample(sample_dict=None):
    """
        Validate a paired sample.

        Parameters
        ----------
        sample_dict : dict
            Dictionary containing sample information

        Returns
        -------
        dict
            Dictionary: {"valid": bool, "message": str}
        """
    log_text = ""
    valid = True

    # for key in sample_dict:
    #     print(key, sample_dict[key])

    if sample_dict["Vector"]:
        if not os.path.exists(sample_dict["Vector"]):
            log_text += " Vector file does not exist!"
            valid = False

    if not sample_dict["FastQ read 1"]:
        log_text += " Fastq file 1 does not exist!"
        valid = False
    if not os.path.exists(sample_dict["FastQ read 1"]):
        log_text += " FastQ file 1 does not exist!"
        valid = False
    if not sample_dict["FastQ read 2"]:
        log_text += " Fastq file 2 does not exist!"
        valid = False
    if not os.path.exists(sample_dict["FastQ read 2"]):
        log_text += " FastQ file 2 does not exist!"
        valid = False
    if not sample_dict["Lib 1"]:
        log_text += " Library file 1 does not exist!"
        valid = False
    if not sample_dict["Lib 2"]:
        log_text += " Library file 2 does not exist!"
        valid = False

    if not (sample_dict["Expected reads"] == 0 or isinstance(sample_dict["Expected reads"], (int, numpy.int64))):
        log_text += " Expected reads must be an integer!"
        valid = False

    if valid and sample_dict["FastQ read 1"] == sample_dict["FastQ read 2"]:
        log_text += " Read 1 and read 2 are identical, exclude!"
        valid = False

    return {"valid": valid, "message": log_text}
