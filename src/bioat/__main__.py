"""main module of bioat.

The cli interface of bioat.
    - `python -m bioat -h` or `bioat -h` to see the exposed interface.
"""
import fire
from bioat.cli import Cli


def main() -> int:
    calculator = Cli()
    # print in shell stdout instead of view in `less`
    fire.core.Display = lambda lines, out: print(*lines, file=out)
    #
    fire.Fire(calculator, name='bioat')


if __name__ == '__main__':
    main()
