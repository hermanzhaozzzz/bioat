import bioat


def test_about():
    assert bioat.about is not None
    assert isinstance(bioat.about, str)
    print(bioat.about)