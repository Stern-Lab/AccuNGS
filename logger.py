import os
import logging
import subprocess


def _logger_already_exists(logger, log_file):
    if logger.hasHandlers(): # If logger exists, just return the existing logger.
        if logger.handlers[1].baseFilename != log_file:
            if log_file is not None:
                logger.warning(f"Logger already exists! Sticking with log file: {logger.handlers[1].baseFilename}")
        return_value = True
    else:
        return_value = False
    return return_value

def _create_new_logger(logger, log_file):
    for handler in logging.root.handlers[:]:  # the logger sometimes doesn't write the log files and this may help...
        logging.root.removeHandler(handler)   # TODO: does it help?
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to info
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create file handler and set level to debug
    if not os.path.exists(log_file): #create file if it doesn't exist.
        open(log_file, 'w').close()
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # add formatter
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add ch & fh to logger
    logger.addHandler(ch)
    logger.addHandler(fh)
    logger.info(f'Log started! Outputing to: {log_file}')
    # get the git hash of this directory.
    git_hash = subprocess.check_output(["git", "describe", "--always"], cwd=os.path.dirname(__file__)).strip().decode()
    logger.debug(f"git hash: {git_hash}")

    return logger

def pipeline_logger(logger_name, log_folder=None):
    # create logger
    if log_folder is not None:
        if not os.path.exists(log_folder):
            os.mkdir(log_folder)
        log_file = os.path.join(log_folder, '.log')
    else:
        log_file = None
    logger = logging.getLogger(logger_name)
    if not _logger_already_exists(logger, log_file):
        if log_file is None:
            raise ValueError("First instance of logger must be initiated with an output file!")
        logger = _create_new_logger(logger, log_file)
    return logger


def aggregate_logs(log_folder, log):
    # TODO: agg all logs in subfolders. stop if log as 'aggregated' in its name.
    pass


def run_with_logger(function_to_run, log_name, log_folder):
    # maybe run mains through something like this to log exceptions and agg logs.
    log = pipeline_logger(logger_name=log_name, log_folder=log_folder)
    try:
        function_to_run
    except Exception as e:
        log.error(e)
    finally:
        aggregate_logs(log_folder, log)
