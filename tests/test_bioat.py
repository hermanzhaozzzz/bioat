import bioat


def test_version_info():
    assert bioat.__version__ == "0.10.0"
    assert bioat.__author__ == "Huanan Herman ZHAO @ Tsinghua University"
    assert bioat.__doc_format__ == "restructuredtext"
