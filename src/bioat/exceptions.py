"""
This module contains custom exception classes for the bioat package.
"""

__module_name__ = "bioat.exceptions"

__all__ = [
    "BioatError",
    "BioatFileFormatError",
    "BioatFileNotCompleteError",
    "BioatFileNotFoundError",
    "BioatInvalidInputError",
    "BioatInvalidOptionError",
    "BioatInvalidParameterError",
    "BioatMissingDependencyError",
    "BioatRuntimeError",
    "BioatValueError",
]


class BioatError(Exception):
    """
    Base class for all custom exceptions in the bioat package.
    """

    def __init__(self, *args, **kwargs):
        # 使用列表来收集消息，避免多次字符串拼接
        messages = []
        if args:
            messages.append(" ".join(map(str, args)))  # 将args转换为字符串并连接
        if kwargs:
            messages.append(
                " ".join(f"{k}={v}" for k, v in kwargs.items())
            )  # 格式化kwargs
        self.msg = "; ".join(messages)

    def __str__(self):
        """
        Return a string representation of the exception.
        """
        return f"{self.msg}"  # 返回当前类的名称


class BioatFileFormatError(BioatError):
    """
    Exception raised for errors in file format.
    """
    pass


class BioatFileNotCompleteError(BioatError):
    """
    Exception raised when a file is not complete.
    """
    pass


class BioatFileNotFoundError(BioatError):
    """
    Exception raised when a required file is not found.
    """
    pass


class BioatInvalidInputError(BioatError):
    """
    Exception raised for invalid input errors.
    """
    pass


class BioatInvalidOptionError(BioatError):
    """
    Exception raised for invalid option errors.
    """
    pass

class BioatInvalidParameterError(BioatError):
    """
    Exception raised for invalid parameter errors.
    """

    pass


class BioatMissingDependencyError(BioatError):
    """
    Exception raised when a required dependency is missing.
    """
    pass


class BioatRuntimeError(BioatError):
    """
    Exception raised for runtime errors.
    """
    pass


class BioatValueError(BioatError):
    """
    Exception raised for errors related to values.
    """
    pass