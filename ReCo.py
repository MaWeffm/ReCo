"""Main entry for ReCo."""

import argparse
import sys

from reco.reco import ReCo
from reco._version import get_versions


def gui(subparser):
    gui_command = subparser.add_parser("gui", help="Available soon: Run ReCo in a GUI")

    optional = gui_command.add_argument_group(
        title="Optional arguments", description=""
    )

    optional.add_argument(
        "-r", "--remove_unused_files", action="store_true", help="Remove unused files."
    )
    optional.add_argument("-j", "--cores", type=int, help="The number of cores to use.")


def cli(subparser):
    cli_command = subparser.add_parser("cli", help="Run ReCo from the command line")

    required = cli_command.add_argument_group(
        title="Required arguments", description=""
    )
    required.add_argument(
        "-s",
        "--sample_sheet",
        type=str,
        required=True,
        help="A sample sheet containing sample information. See documentation for file format.",
    )
    required.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="The path to an output dir. Will be created if it does not exist.",
    )

    optional = cli_command.add_argument_group(
        title="Optional arguments", description=""
    )

    optional.add_argument(
        "-r", "--remove_unused_files", action="store_true", help="Remove unused files."
    )
    optional.add_argument("-j", "--cores", type=int, help="The number of cores to use.")


def parse_args():
    """
    Parse command line arguments.

    Returns
    -------
    args : Namespace
        Namespace populated with attributes.
    """
    parent_parser = argparse.ArgumentParser(
        description=f"ReCo v{get_versions()['version']}: find gRNA read counts (ReCo) in fastq files."
    )
    parent_parser.add_argument(
        "-v", "--version", action="version", version=f"{get_versions()['version']}"
    )
    parent_parser.add_argument(
        "-r", "--remove_unused_files", action="store_true", help="Remove unused files."
    )
    parent_parser.add_argument(
        "-j", "--cores", type=int, help="The number of cores to use."
    )

    subparser = parent_parser.add_subparsers(
        help="Commands to run different flavors of ReCo",
        title="Subcommands",
        dest="subcmd",
    )

    gui(subparser)
    cli(subparser)

    args = parent_parser.parse_args()
    return args


def main():
    """
    ReCo main function.

    Returns
    -------

    """
    args = parse_args()

    if args.subcmd == "gui":
        print("Running GUI")
    elif args.subcmd == "cli":
        r = ReCo(
            sample_sheet_file=args.sample_sheet,
            output_dir=args.output_dir,
        )
        r.run(remove_unused_files=args.remove_unused_files, cores=args.cores)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Terminated!\n")
        sys.exit(0)
