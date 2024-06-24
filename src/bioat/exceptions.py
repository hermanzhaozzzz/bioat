__module_name__ = 'bioat.exceptions'


class BioatException(ValueError):
    def __init__(self, *args, **kwargs):
        self.msg = ''
        if args:
            self.msg += f'{args};'
        if kwargs:
            self.msg += f'{kwargs};'

    def __str__(self):
        return f"BioatException: {self.msg}"


class BioatFileFormatError(BioatException):
    def __str__(self):
        return f"BioatFileFormatError: {self.msg}"


class BioatFileNotCompleteError(Exception):
    def __str__(self):
        return f"BioatFileNotCompleteError: {self.msg}"


class BioatParameterFormatError(Exception):
    def __str__(self):
        return f"BioatParameterFormatError: {self.msg}"


class BioatFileNameError(Exception):
    def __str__(self):
        return f"BioatFileNameError: {self.msg}"


class BioatCloudflareChallengeFailure(Exception):
    """自定义异常类，表示Cloudflare挑战失败"""

    def __str__(self):
        return f"BioatCloudflareChallengeFailure: {self.msg}"

class BioatRuntimeWarning(UserWarning):
    """自定义警告类，表示运行时出现的警告"""
    def __str__(self):
        return f"BioatRuntimeWarning: {self.msg}"
