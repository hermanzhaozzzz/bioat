import fire
from bioat import __version__


class Calculator(object):

    def add(self, x, y):
        return x + y

    def multiply(self, x, y):
        return x * y


def main() -> int:
    calculator = Calculator()
    fire.Fire(calculator)


if __name__ == "__main__":
    main()
