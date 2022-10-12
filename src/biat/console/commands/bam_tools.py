from cleo.commands.command import Command
from cleo.helpers import argument, option


class BamToolsCommand(Command):
    name = "bamtools"
    description = "Toolkit for bam file"
    arguments = [
        argument(
            "to-pmat",
            description="Convert a coordinated bam to a pmat.tsv file",
            optional=True
        )
    ]
    options = [
        option(
            long_name="--input",
            short_name="-i",
            description="Input bam file, sorted and indexed",
            flag=True,
            value_required=True,
            multiple=False,
            default=None,
        )
    ]

    def handle(self):
        name = self.argument("bamtools")

        if name:
            text = f"Hello {name}"
        else:
            text = "Hello"

        if self.option("yell"):
            text = text.upper()

        self.line(text)
