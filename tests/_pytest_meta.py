import tempfile
from pathlib import Path

HOME = Path.home()
MODULE_PATH = Path("src")
DATA_PATH = Path("data")


if __name__ == "__main__":
    print(HOME)
    print(type(HOME))
    print(DATA_PATH)
