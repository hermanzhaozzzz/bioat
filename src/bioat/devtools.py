import os
import time

import psutil

from bioat.exceptions import BioatInvalidParameterError
from bioat.logger import get_logger

__module_name__ = "bioat.devtools"


def _elapsed_since(start):
    return time.strftime("%Mm:%Ss", time.gmtime(time.time() - start))


def _get_process_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss


def profile(num_iterations=1):
    """
    Decorator for calculating the runtime of a function, over multiple iterations.

    :param num_iterations: times to iter, defaults to 1
    :type num_iterations: int, optional
    """

    def inner_decorator(func):
        def wrapper(*args, **kwargs):
            if not isinstance(num_iterations, int) or num_iterations < 1:
                raise BioatInvalidParameterError(
                    "num_iterations must be a positive integer"
                )
            logger = get_logger(
                level="INFO", module_name=__name__, func_name=func.__name__
            )

            total_time = 0
            start_time = 0
            result = None

            mem_before = _get_process_memory()
            for _ in range(num_iterations):
                start_time = time.time()
                result = func(*args, **kwargs)
                elapsed_time = time.time() - start_time
                total_time += elapsed_time

            mem_after = _get_process_memory()
            mem_diff_mb = (mem_after - mem_before) / (1024**2)  # 转换为 MB

            avg_time = total_time / num_iterations
            elapsed_time_formatted = _elapsed_since(start_time)

            text = (
                f"\n[{func.__name__}] performance profiling:\n"
                f"\t↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓\n"
                f"\tIteration num: {num_iterations}x\n"
                f"\tTotal    time: {total_time:.9f}s ({elapsed_time_formatted})\n"
                f"\tAverage  time: {avg_time:.9f}s\n"
                f"\tMemory   diff: {mem_diff_mb:.2f} MB"
            )
            logger.info(text)

            return result

        return wrapper

    return inner_decorator