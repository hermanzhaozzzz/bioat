from tests import bioat


def test_version():
    assert bioat.__pysamstats_version__ >= '1.1.2'
    assert bioat.__pysam_version__ >= '0.19.1'
    assert bioat.__version__ == '0.1.1.4'


def test_author():
    assert bioat.__author__ is not None
    assert bioat.__email__ is not None
