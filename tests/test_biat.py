import biat


def test_version():
    assert biat.__version__ == '0.1.0'


def main():
    test_version()


if __name__ == '__main__':
    main()