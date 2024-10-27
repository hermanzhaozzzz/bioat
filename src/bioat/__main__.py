"""Main module of the bioat package.

This module serves as the command-line interface (CLI) for the bioat package.
It uses the Fire library to automatically generate help and command-line parsing capabilities.

Usage:
    To see the available commands and options, run:
        python -m bioat -h
    or simply:
        bioat -h

Attributes:
    calculator (Cli): An instance of the Cli class that handles command-line operations.

Functions:
    main() -> int: The entry point of the module that sets up and runs the CLI.
"""

import fire

from bioat.cli import Cli


def main() -> None:
    calculator = Cli()
    # print in shell stdout instead of view in `less`
    fire.core.Display = lambda lines, out: print(*lines, file=out)
    #
    fire.Fire(calculator, name='bioat')


if __name__ == '__main__':
    print(type(main()))
