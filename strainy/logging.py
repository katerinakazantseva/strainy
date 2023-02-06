import os
import logging

logger = logging.getLogger()


def set_thread_logging(log_dir, log_prefix, thread_id):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    #clear all handlers if we rerun this in a new thread
    logger.handlers.clear()
    #thread_id = str(multiprocessing.current_process().name).split("-")[-1]
    #thread_id = str(multiprocessing.current_process().pid)
    if thread_id:
        log_file = os.path.join(log_dir, "{0}-{1}.log".format(log_prefix, thread_id))
        thread_tag = "[Thread {0}]".format(thread_id)
    else:
        log_file = os.path.join(log_dir, "{0}.log".format(log_prefix))
        thread_tag = "[Root]"

    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] " + thread_tag + " %(levelname)s: "
                                          " %(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)
    #logger.setLevel(logging.INFO)
    console_log.setLevel(logging.INFO)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)







