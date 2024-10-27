"""
Unit tests for the `bioat about` subcommand.

This module contains tests for the `about` functionality of the
`bioat` command-line interface (CLI). It includes two main tests:

1. `test_about`: Verifies that the `about` method exists on the
   `Cli` object and that it returns a string.

2. `test_cli_about`: Tests the execution of the `bioat about`
   command through the command-line interface, ensuring that it
   runs without any errors.

It utilizes the `subprocess` module to run CLI commands and
asserts expected behaviors using Python's built-in `assert` statements.
"""

import subprocess

from bioat.cli import Cli

cli = Cli()


def test_about():
    """Test the about function of the CLI.

    This function checks:
    1. If the `cli` object has an attribute `about`.
    2. If the `about` method returns a string.
    """
    assert hasattr(cli, "about")  # Checking if the 'about' method exists
    assert isinstance(cli.about(), str)  # Ensuring 'about' method returns a string
    # Test about text attribute
    assert isinstance(cli.about(), str)  # Rechecking that 'about' returns a string


def test_cli_about():
    """Test the 'about' command via command-line interface.

    This function executes the command `bioat about` and checks if it runs without errors.
    """
    args = ["bioat", "about"]  # Command-line arguments for the 'about' command
    assert subprocess.run(
        args, check=True
    )  # Running the command and asserting it completes successfully
