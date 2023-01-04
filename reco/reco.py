# pylint: disable=line-too-long

"""ReCo class."""

import datetime
import locale
import os

from .logger import setup_logger
from .sample_sheet import SampleSheet

from ._version import get_versions

locale.setlocale(locale.LC_ALL, "en_US.UTF-8")


class ReCo:
    """
    Class representing ReCo.

    Parameters
    ----------
    sample_sheet_file : str
        Path to a sample sheet file.
    output_dir : str
        Output folder to which all output is written.

    Methods
    -------
    run
        Run all samples.
    """

    reco_log = "ReCo.log"

    def __init__(self, sample_sheet_file=None, output_dir=None):
        self.date = datetime.datetime.now().replace(microsecond=0)
        self.output_dir = output_dir
        self.logfile = os.path.join(self.output_dir, self.reco_log)

        self.logger = setup_logger(
            f"ReCo {get_versions()['version']}", self.logfile, stdout=True
        )
        self.logger.info("Initializing %s", self)
        self.logger.info("Sample sheet: %s", sample_sheet_file)

        self.sample_sheet = SampleSheet.from_file(
            logger=self.logger,
            sample_sheet_file=sample_sheet_file,
            output_dir=self.output_dir,
        )

    @property
    def output_dir(self):
        """Get the output dir."""
        return self._output_dir

    @output_dir.setter
    def output_dir(self, new_output_dir):
        """Set output dir, create if it does not exist."""
        try:
            os.makedirs(new_output_dir, exist_ok=False)
        except FileExistsError as exc:
            raise FileExistsError(
                f"Output dir '{new_output_dir}' already exists!"
            ) from exc
        self._output_dir = new_output_dir

    def run(self, remove_unused_files=True, cores=1):
        """
        Run all samples.

        Parameters
        ----------
        remove_unused_files : bool, default=True
            Whether to remove intermediate files to save disk space. Intermediate files can be quite large, consider setting this option to True!
        cores : int, default=1
            Number of cores to utilize.

        Returns
        -------
        bool
            Whether ReCo successfully processed all samples.

        """
        self.logger.info("Running on %s cores", cores)
        self.logger.info("Remove intermediate files: %s", remove_unused_files)
        valid_sample_names = list(self.sample_sheet.samples.values())
        self.logger.info(
            "Found %s valid samples: %s", len(valid_sample_names), valid_sample_names
        )

        for sample in self.sample_sheet:
            for step in sample.run(
                remove_unused_files=remove_unused_files, cores=cores
            ):
                print(" ", step)

    def __repr__(self):
        return f"<{self.__class__.__name__} v{get_versions()['version']}>"

    def __str__(self):
        return f"{self.__class__.__name__} v{get_versions()['version']}, {self.date}"
