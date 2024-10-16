import bioat


def test_version_info():
    assert bioat.__version__ == "0.12.14"
    assert bioat.__author__ == "Huanan Herman ZHAO @ Tsinghua University"
    assert bioat.__doc_format__ == "sphinx"
    # !same with vscode setting
    # !"autoDocstring.docstringFormat": "sphinx",
