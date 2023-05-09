import bioat

def test_version_info():
    assert bioat.__version__ == '0.1.5.3'
    assert bioat.__author__ == 'Hua-nan ZHAO @ Tsinghua University'
    assert bioat.__doc_format__ == "restructuredtext"

def test_import():
    from bioat import about, set_logging_level
    from bioat import Bam, Fastx, HiC, Mgi, System, Table, TargetedDeepSequencing
