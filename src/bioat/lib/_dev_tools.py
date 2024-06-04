import os
import sys
import time

import psutil

from bioat.logger import get_logger

__module_name__ = "bioat.lib._dev_tools"


def _elapsed_since(start):
    return time.strftime("%Mm:%Ss", time.gmtime(time.time() - start))


def _get_process_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss


# def profile(func):
#     """Decorator for calculate runtime of a function."""

#     def wrapper(*args, **kwargs):
#         logger = get_logger(
#             level="INFO", module_name=__module_name__, func_name="profile"
#         )
#         mem_before = _get_process_memory()
#         start = time.time()
#         result = func(*args, **kwargs)
#         elapsed_time = _elapsed_since(start)
#         mem_after = _get_process_memory()
#         logger.info("{}: \n\tmemory before: {:,}, after: {:,}, consumed: {:,}; \n\texec time: {}".format(
#             func.__name__,
#             mem_before, mem_after, mem_after - mem_before,
#             elapsed_time))
#         return result


#     return wrapper
def profile(num_iterations=1):
    """Decorator for calculating the runtime of a function, over multiple iterations.

    :param func: function to be decorated.
    :type func: function
    :param num_iterations: times to iter, defaults to 1
    :type num_iterations: int, optional
    """

    def inner_decorator(func):
        def wrapper(*args, **kwargs):
            logger = get_logger(
                level="INFO", module_name=__name__, func_name=func.__name__
            )

            total_time = 0
            for i in range(num_iterations):
                mem_before = _get_process_memory()
                start_time = time.time()
                result = func(*args, **kwargs)
                elapsed_time = time.time() - start_time
                total_time += elapsed_time
                mem_after = _get_process_memory()

            avg_time = total_time / num_iterations
            text = (
                f"\n[{func.__name__}] performance profiling:\n\t↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓\n"
                f"\tIteration num: {num_iterations}x\n"
                f"\tTotal    time: {total_time:.9f}s\n"
                f"\tAverage  time: {avg_time:.9f}s"
            )
            logger.info(text)

            return result

        return wrapper

    return inner_decorator
