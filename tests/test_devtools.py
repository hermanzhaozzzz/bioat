import time

import pytest

from bioat.devtools import profile
from bioat.exceptions import BioatInvalidParameterError


def test_profile_happy_path():
    @profile(num_iterations=5)
    def sample_function():
        time.sleep(0.1)
        return "done"

    result = sample_function()
    assert result == "done"


def test_profile_edge_case_negative_iterations():
    # 使用 pytest.raises 捕获在调用时抛出的异常
    with pytest.raises(BioatInvalidParameterError):

        @profile(num_iterations=-1)
        def sample_function():
            return "done"

        # 执行被装饰的函数以触发异常
        sample_function()


def test_profile_edge_case_large_iterations():
    @profile(num_iterations=1000)
    def sample_function():
        return "done"

    result = sample_function()
    assert result == "done"
