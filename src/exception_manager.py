import time
import psutil
import os


class Error(Exception):
    """Base class for other exceptions"""
    pass


class SamePositionError(Error):
    """Ends of ambiguous and unambiguous groups coincide"""
    pass


class NotSameLengthError(Error):
    """The two deques are not of the same size"""
    pass


class InvalidStart(Error):
    """The starting address is at an invalid spot"""
    pass


class InsufficientArguments(Error):
    """Insufficient number of addresses passed."""
    pass


def get_time(past:time, loud=True, print=print):
    right_now = time.time()
    if loud:
        print('Time Elapsed: %s', str(right_now-past))
    return right_now


PROCESS = psutil.Process(os.getpid())


def print_memory_usage():
    """Prints current memory usage stats.

    inspired by: https://stackoverflow.com/a/15495136

    :return: None
    """
    total, available, percent, used, free = psutil.virtual_memory()
    proc = PROCESS.memory_info()[1]
    print('process = %s total = %s available = %s used = %s free = %s percent = %s'
          % (proc, total, available, used, free, percent))

if __name__ == "__main__":
    try:
        error_message = "Error: do not run " + __file__ + " directly"
        raise Exception(error_message)

    except FileNotFoundError:
        raise Exception("Error: do not run exception_manager.py directly")

