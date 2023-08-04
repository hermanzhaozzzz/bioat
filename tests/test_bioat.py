import bioat

def test_version_info():
    assert bioat.__version__ == '0.2.9'
    assert bioat.__author__ == 'Hua-nan ZHAO @ Tsinghua University'
    assert bioat.__doc_format__ == "restructuredtext"
