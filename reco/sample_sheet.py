# pylint: disable=line-too-long

"""Sample sheet class for ReCo."""

from .attribute_descriptors import FileDesc
from .attribute_descriptors import LoggerDesc
from .sample import PairedSample
from .sample import SingleSample
from .sample_tools import SampleValidator
from .tools import read_df
from .tools import remove_duplicates


class SampleSheet:
    """Class representing a sample sheet.

    Parameters
    ----------
    logger : logging.Logger
        A logging.Logger object
    sample_sheet_file : str
        Path to the sample sheet file in one of the following formats:
        "*.xlsx", "*.csv", "*.tsv", "*.txt".

    Methods
    -------
    """
    sample_sheet_file = FileDesc()
    logger = LoggerDesc()

    def __init__(self, logger=None, sample_sheet_file=None, output_dir=None):
        self.logger = logger
        self.sample_sheet_file = sample_sheet_file

        self.samples = {}

        self.sample_sheet_df = None

        self.output_dir = output_dir

        self._iter_current = 0
        self._iter_max = len(self.samples)

    @classmethod
    def from_file(cls, logger=None, sample_sheet_file=None, output_dir=None):
        """Build sample sheet instance directly from sample sheet file."""
        s_sheet = cls(logger=logger, sample_sheet_file=sample_sheet_file, output_dir=output_dir)
        s_sheet.read_sample_sheet_file()
        return s_sheet

    def read_sample_sheet_file(self):
        """
        Read the sample sheet in one of the following formats: "*.xlsx", "*.csv", "*.tsv", "*.txt".

        Returns
        -------

        """
        sample_sheet_df = read_df(
            self.logger, df_file=self.sample_sheet_file, columns=range(0, 10), header=0
        )
        cleaned_sample_sheet_df = remove_duplicates(
            self.logger, dataframe=sample_sheet_df, columns=["Sample name"]
        )
        cleaned_sample_sheet_df.reset_index(drop=True, inplace=True)

        cleaned_sample_sheet_df.fillna("", inplace=True)
        cleaned_sample_sheet_df.replace("", None, inplace=True)

        cleaned_sample_sheet_df["Sample name"] = cleaned_sample_sheet_df["Sample name"].str.translate(
            {ord(c): "_" for c in r"!@#$%^&*()[]{};:,./<>?\|`~-=+ "}
        )  # replace special characters with "_", otherwise file system/platform stuff will become ugly

        self.sample_sheet_df = cleaned_sample_sheet_df

        self.logger.info(f"Unique samples: {list(cleaned_sample_sheet_df['Sample name'].values)}")

        self.create_samples()

    def create_samples(self):
        """
        Validate the samples in the sample sheet.

        Returns
        -------

        """
        sample_counter = 0
        for index in self.sample_sheet_df.index:
            entry = self.sample_sheet_df.loc[index]
            sample_dict = {key: entry[key] for key in entry.keys()}
            valid = SampleValidator.from_dict(sample_dict)

            if valid["valid"]:
                if sample_dict["Sample type"] == "single":
                    self.samples[sample_counter] = SingleSample(sample_dict=sample_dict, output_dir=self.output_dir)
                    sample_counter += 1
                elif sample_dict["Sample type"] == "paired":
                    self.samples[sample_counter] = PairedSample(sample_dict=sample_dict, output_dir=self.output_dir)
                    sample_counter += 1
            else:
                self.logger.warning(valid)
        self._iter_max = len(self.samples)

    def __len__(self):
        """
        Number of samples.

        Returns
        -------
        int
            The number of samples in the sample sheet.
        """
        return len(self.samples)

    def __iter__(self):
        """Iterator."""
        self._iter_current = 0
        return self

    def __next__(self):
        """Iterator next."""
        if self._iter_current < self._iter_max:
            self._iter_current += 1
            return self.samples[self._iter_current - 1]
        raise StopIteration

    def __repr__(self):
        return f"<{self.__class__.__name__}"

    def __str__(self):
        return f"{self.__class__.__name__}"
