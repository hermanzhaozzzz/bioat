import time
import os
import sys
import psutil
from bioat import get_logger

__module_name__ = 'bioat.lib._dev_tools'


def _elapsed_since(start):
    return time.strftime("%Mm:%Ss", time.gmtime(time.time() - start))


def _get_process_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss


def profile(func):
    """Decorator for calculate runtime of a function."""

    def wrapper(*args, **kwargs):
        logger = get_logger(level='INFO', module_name=__module_name__,
                            func_name=sys._getframe().f_code.co_name)
        mem_before = _get_process_memory()
        start = time.time()
        result = func(*args, **kwargs)
        elapsed_time = _elapsed_since(start)
        mem_after = _get_process_memory()
        logger.info("{}: \n\tmemory before: {:,}, after: {:,}, consumed: {:,}; \n\texec time: {}".format(
            func.__name__,
            mem_before, mem_after, mem_after - mem_before,
            elapsed_time))
        return result

    return wrapper
