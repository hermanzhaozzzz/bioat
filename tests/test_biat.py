import biat


def test_version():
    assert biat.__pysamstats_version__ >= '1.1.2'
    assert biat.__pysam_version__ >= '0.19.1'
    assert biat.__version__ == '0.1.0'
    assert biat.__author__ is not None
    assert biat.__email__ is not None
