# pylint: disable=line-too-long

"""Read count writer class for ReCo."""

from collections import Counter

import datetime
import pysam
import pandas as pd
import tqdm

from .attribute_descriptors import LoggerDesc
from .reco_config import Config
from .tools import count_sam_reads


class ReadcountWriter:
    """
    Class representing a read count writer, writes table from SAM and library files.

    Methods
    -------
    counts_from_sam
    write_counter
    """
    logger = LoggerDesc()

    def __init__(self, logger=None):
        self.logger = logger

    def counts_from_sam(self, sam_file=None, library=None):
        """
        Determine gRNA counts from SAM and library files.

        Parameters
        ----------
        sam_file : str
            Path to a SAM file.
        library : library.Library
            A Library object.

        Returns
        -------
        dict
            {"guide_hitlist": guide_hitlist, "failed": failed}

        """
        count_start = datetime.datetime.now()
        guide_hitlist = Counter({guide["guide_name"]: 0 for guide in library})
        sequence_hitlist = {}
        failed = Counter()

        sam_reads = count_sam_reads(sam_file)
        sam_handle = pysam.AlignmentFile(sam_file, "rb")
        self.logger.info("Alignments: %s", f"{sam_reads:n}")
        for alignment in tqdm.tqdm(sam_handle.fetch(), desc=" Fetching Bowtie2 alignments", total=sam_reads, unit=" alignments"):
            alignment_info = test_alignment(alignment)

            reference_name = alignment_info["reference_name"]
            sequence = alignment_info["sequence"]
            count = alignment_info["count"]
            status = alignment_info["status"]
            score = alignment_info["score"]

            if status in ("primary", "tolerate"):
                guide_hitlist[reference_name] += count
                sequence_hitlist[sequence] = reference_name
            else:
                failed[f"{sequence},{status},{score},{reference_name}"] += count

        sam_handle.close()
        print("   Counted in", datetime.datetime.now() - count_start)
        return {"guide_hitlist": guide_hitlist, "sequence_hitlist": sequence_hitlist, "failed": failed}

    def write_counter(self, count_tables=None, count_col_name=None, counts_file=None, failed_file=None, top_failed_file=None):
        """
        Write count tables to files.

        Parameters
        ----------
        count_tables : dict
        count_col_name : str
            Name of the column that contains the count.
        counts_file : str
            Path to the final gRNA count file.
        failed_file : str
            Path to the file containing failed sequence alignments.
        top_failed_file : str
            Path to the file containing a short list of failed sequence alignments.

        Returns
        -------

        """
        write_start = datetime.datetime.now()
        self.logger.info("Writing read counts to: %s", counts_file)
        self.logger.info("Writing failed alignments to: %s", failed_file)
        counts = count_tables["guide_hitlist"]
        failed = count_tables["failed"]

        read_count_table = pd.DataFrame.from_dict(
            counts, orient="index", columns=[count_col_name]
        )

        read_count_table.sort_values(by=count_col_name, ascending=False, inplace=True)
        read_count_table.index.names = ["Guide"]
        read_count_table.to_csv(counts_file)

        with open(failed_file, "w", encoding="utf-8") as writefile, open(top_failed_file, "w", encoding="utf-8") as top_failed_file_writer:
            writefile.write("sequence,status,score,nearest gRNA,count\n")
            top_failed_file_writer.write("sequence,status,score,nearest gRNA,count\n")
            failed_seq_counter = 0
            for seq, count in sorted(failed.items(), key=lambda item: item[1], reverse=True):
                writefile.write(f"{seq},{failed[seq]}\n")
                if failed_seq_counter < Config.TOP_FAILED_SEQUENCES:
                    top_failed_file_writer.write(f"{seq},{failed[seq]}\n")
                failed_seq_counter += 1
        print("   Table written in", datetime.datetime.now() - write_start)
        return read_count_table

    def __repr__(self):
        """Representation."""
        return f"{self.__class__.__name__}"

    def __str__(self):
        """String representation."""
        return f"<{self.__class__.__name__}>"


def test_alignment(alignment):
    """
    Get alignments from an alignment.

    Parameters
    ----------
    alignment : pysam.AlignedSegment
        A pysam.AlignedSegment object.

    Returns
    -------
    dict
        {"reference_name": alignment.reference_name, "sequence": sequence, "count": count, "status": status, "score": as_tag}
    """
    my_alignment_split = str(alignment).split('\t', maxsplit=1)[0]
    count = int(my_alignment_split.split("_")[1])
    sequence = my_alignment_split.split("_")[0]
    as_tag = ""

    if alignment.has_tag("AS"):
        as_tag = alignment.get_tag("AS")
        if as_tag >= 2 * len(sequence):  # max possible alignment score
            if alignment.has_tag("XS"):
                xs_tag = alignment.get_tag("XS")
                if xs_tag <= as_tag - 2:
                    status = "tolerate"
                else:
                    status = "ambiguous"
            else:
                status = "primary"
        else:
            status = "insufficient"
    else:
        status = "failed"

    return {"reference_name": alignment.reference_name, "sequence": sequence, "count": count, "status": status, "score": as_tag}
