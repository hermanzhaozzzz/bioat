#!/usr/bin/env python

from biat.console.commands.bam_tools import BamToolsCommand

from cleo.application import Application


def main() -> int:
    application = Application()
    application.add(BamToolsCommand())
    exit_code: int = application.run()
    return exit_code


if __name__ == "__main__":
    main()