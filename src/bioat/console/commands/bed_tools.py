from cleo.commands.command import Command
from cleo.helpers import argument, option


class BedToolsCommand(Command):
    name = "bedtools"
    description = "Toolkit for bed file"
    arguments = [
        argument(
            "sort",
            description="Sort a bed file",
            optional=True
        )
    ]
    options = [
        option(
            long_name="--input",
            short_name="-i",
            description="Input bed file",
            flag=True,
            value_required=True,
            multiple=False,
            default=None,
        )
    ]

    def handle(self):
        name = self.argument("sort")
        self.line(name)

