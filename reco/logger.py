"""Provide logging functionality for ReCo"""
import logging
import sys


def setup_logger(name, logfile, level=logging.INFO, stdout=False):
    """Set up module/function-specific logger,

    Parameters
    ----------
    name : str
        Name of logger for easy reference.
    logfile : str
        Path to the log file.
    level
        Logging level
    stdout : bool, default=False
        Switch on/off printing to stdout.
    Returns
    -------
    logging.Logger

    """
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S"
    )

    handler = logging.FileHandler(logfile)
    handler.setFormatter(formatter)

    logging.getLogger("matplotlib.font_manager").disabled = True
    logging.getLogger("matplotlib.backends.backend_pdf").disabled = True
    logging.getLogger("matplotlib.ticker").disabled = True

    logger = logging.getLogger(name)
    logger.setLevel(level)

    file_handler = logging.FileHandler(filename=logfile, mode="w")
    file_handler.setFormatter(formatter)

    logger.setLevel(level)
    logger.addHandler(file_handler)

    if stdout:
        # Add a stream handler if logging to stdout is desired.
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setFormatter(formatter)
        stdout_handler.setLevel(logging.DEBUG)
        logger.addHandler(stdout_handler)

    return logger
