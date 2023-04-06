import time
import os
import psutil
import logging


def _elapsed_since(start):
    return time.strftime("%Mm:%Ss", time.gmtime(time.time() - start))


def _get_process_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss


def profile(func):
    """Decorator for calculate runtime of a function."""

    def wrapper(*args, **kwargs):
        mem_before = _get_process_memory()
        start = time.time()
        result = func(*args, **kwargs)
        elapsed_time = _elapsed_since(start)
        mem_after = _get_process_memory()
        logging.info("{}: \n\tmemory before: {:,}, after: {:,}, consumed: {:,}; \n\texec time: {}".format(
            func.__name__,
            mem_before, mem_after, mem_after - mem_before,
            elapsed_time))
        return result

    return wrapper
