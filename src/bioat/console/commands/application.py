#!/usr/bin/env python
from bioat import __version__
from cleo.application import Application as BaseApplication
from bioat.console.commands.bam_tools import BamToolsCommand
from bioat.console.commands.bed_tools import BedToolsCommand


class Application(BaseApplication):
    def __init__(self) -> None:
        super().__init__(name="bioat", version=__version__)

        for command in self.get_default_commands():
            self.add(command)

    def get_default_commands(self) -> list:
        commands = [
            BamToolsCommand(),
            BedToolsCommand(),
        ]
        return commands


def main() -> int:
    application = Application()
    exit_code: int = application.run()
    return exit_code


if __name__ == "__main__":
    main()
