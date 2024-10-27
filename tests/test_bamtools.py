"""Tests for `bamtools` package."""

import os
import random
import string
import subprocess

import pytest

from bioat import __BAM_PARSER_BACKEND__
from bioat.cli import Cli
from bioat.lib.libpath import check_cmd

from ._pytest_meta import DATA_PATH

bioat_cli = Cli()

MPILEUP_FILE = os.path.join(DATA_PATH, "bam/test_sorted.mpileup.gz")
BAM_SORTP_FILE = os.path.join(DATA_PATH, "bam/test_sorted.bam")
BAM_SORTN_FILE = os.path.join(DATA_PATH, "bam/test_sorted_n.bam")
TEMP_DIR = (
    f"/tmp/bioat_{''.join(random.sample(string.ascii_letters + string.digits, 16))}"
)


def test_mpileup_to_table():
    """Tests the mpileup to table conversion functionality.

    This function calls the `mpileup2table` method from the `bioat_cli.bam` module,
    converting a specified mpileup file into a table format.

    Args:
        mpileup (str): The path to the input mpileup file (MPILEUP_FILE).
        output (str): The output file path (set to "/dev/null" for no output in this test).
        threads (int): The number of threads to use for processing (based on available CPU cores).
        mutation_number_threshold (int): The threshold for mutations to consider (set to 3).
        temp_dir (str): The directory to use for temporary files (TEMP_DIR).
        remove_temp (bool): Whether to remove temporary files after processing (set to True).
        log_level (str): The logging level for the operation (set to "WARNING").
    """
    bioat_cli.bam.mpileup2table(
        mpileup=MPILEUP_FILE,
        output="/dev/null",
        threads=os.cpu_count() - 1,
        mutation_number_threshold=3,
        temp_dir=TEMP_DIR,
        remove_temp=True,
        log_level="WARNING",
    )


def test_cli_mpileup_to_table():
    """Test the CLI command for converting mpileup to table format.

    This function constructs a command with specified arguments to
    run the mpileup2table command in the bioat tool, checking its
    execution without creating any output file.
    """

    args = [
        "bioat",  # The main command for the bioat tool
        "bam",  # The subcommand indicating the operation type (bam processing)
        "mpileup2table",  # The specific command to convert mpileup to table
        "--mpileup",  # Option to specify the input mpileup file
        MPILEUP_FILE,  # Placeholder for the path to the mpileup file
        "--output",  # Option to specify the output file location
        "/dev/null",  # Indicates that no output file should be created
        "--threads",  # Option to specify the number of threads to use
        str(os.cpu_count() - 1),  # Uses all available CPU cores minus one
        "--mutation_number_threshold",  # Option to set mutation number threshold
        "0",  # Sets the threshold to zero
        "--temp_dir",  # Option to specify the temporary directory for processing
        TEMP_DIR,  # Placeholder for the path to the temporary directory
        "--remove_temp",  # Option indicating whether to remove temporary files
        "True",  # Specifies that temporary files should be removed
        "--log_level",  # Option to set the logging level
        "WARNING",  # Sets the log level to WARNING
    ]

    subprocess.run(
        args, check=True
    )  # Executes the constructed command and checks for errors


def test_remove_clip():
    """Test the remove_clip function from the bam module.

    This test verifies that the remove_clip function works correctly
    when using a backend that supports it. If the backend is set to
    'bamnostic', the test will be skipped since that backend does
    not support the remove_clip functionality.

    Raises:
        pytest.skip: If the BAM parser backend is 'bamnostic'.
    """

    if __BAM_PARSER_BACKEND__ == "bamnostic":
        pytest.skip("bamnostic backend does not support remove_clip")

    bioat_cli.bam.remove_clip(
        input=BAM_SORTN_FILE,
        output="/dev/null",
        threads=os.cpu_count() - 1,
        output_fmt="SAM",
        remove_as_paired=True,
        max_clip=0,
        log_level="WARNING",
    )


def test_cli_remove_clip_sortn():
    """Test the CLI command for removing clips from BAM files sorted by name.

    This function tests the 'remove_clip' command using the specified
    parameters and ensures that the function behaves as expected. It
    skips the test if the bamnostic backend is being used since it
    does not support the 'remove_clip' functionality.

    Skipped conditions:
        - bamnostic backend does not support remove_clip.

    Args:
        None

    Returns:
        None
    """

    if __BAM_PARSER_BACKEND__ == "bamnostic":
        pytest.skip("bamnostic backend does not support remove_clip")

    # Command line arguments for the subprocess run.
    args = [
        "bioat",  # Command name
        "bam",  # Subcommand for BAM operations
        "remove_clip",  # Specific command to remove clips
        "--input",  # Input file option
        BAM_SORTN_FILE,  # Input BAM file (assumed to be defined elsewhere)
        "--output",  # Output file option
        "/dev/null",  # Discard output (no need to save)
        "--threads",  # Threads option for parallel processing
        str(os.cpu_count() - 1),  # Use all but one CPU core
        "--output_fmt",  # Output format option
        "SAM",  # Specify output format as SAM
        "--remove_as_paired",  # Flag to remove as paired
        "True",  # Set the paired removal to True
        "--max_clip",  # Maximum clip size option
        "0",  # Set to 0 to remove all clips
        "--log_level",  # Log level option
        "WARNING",  # Set log level to warning
    ]

    # Execute the command with the specified arguments
    subprocess.run(args, check=True)


def test_cli_remove_clip_sortn_pipe():
    """Test the command line interface for removing clips and sorting BAM files.

    This test verifies that the `remove_clip` functionality works correctly
    when used in a pipeline with `samtools view` and `head`. It skips the test
    if the bamnostic backend is used or if `samtools` is not installed.

    Raises:
        pytest.SkipException: If the bamnostic backend is in use or
        if samtools is not installed.
    """
    if __BAM_PARSER_BACKEND__ == "bamnostic":
        pytest.skip("bamnostic backend does not support remove_clip")
    if not check_cmd("samtools"):
        pytest.skip("samtools is not installed")

    # Prepare the command to view BAM file with samtools
    args = ["samtools", "view", "-h", BAM_SORTN_FILE]
    p1 = subprocess.Popen(args, stdout=subprocess.PIPE)

    # Prepare the command to remove clips from the BAM data
    args = ["bioat", "bam", "remove_clip"]
    p2 = subprocess.Popen(args, stdin=p1.stdout, stdout=subprocess.PIPE)

    # Prepare the command to output the results using head
    args = ["head"]
    p3 = subprocess.Popen(args, stdin=p2.stdout, stdout=subprocess.PIPE)

    # Communicate with the head process and get the output
    _ = p3.communicate()[0]  # output str

    # Verify that the head command executed successfully
    assert p3.returncode == 0


def test_cli_remove_clip_sortn_pipe2():
    """Test the CLI for removing clips with sorted BAM input through a pipeline.

    This function tests the functionality of the 'remove_clip' command
    from the 'bioat' tool when using the 'samtools' to stream input data.
    It verifies that the command operates correctly and does not raise any errors.

    Raises:
        pytest.skip: If the bamnostic backend is used or if samtools is not installed.
    """
    if __BAM_PARSER_BACKEND__ == "bamnostic":
        pytest.skip("bamnostic backend does not support remove_clip")
    if not check_cmd("samtools"):
        pytest.skip("samtools is not installed")

    # Prepare the first command to get sorted BAM data
    args = ["samtools", "view", "-h", BAM_SORTN_FILE]
    p1 = subprocess.Popen(args, stdout=subprocess.PIPE)

    # Setup the remove_clip command from bioat
    args = ["bioat", "bam", "remove_clip"]
    p2 = subprocess.Popen(args, stdin=p1.stdout, stdout=subprocess.PIPE)

    # Prepare the tail command to read the output
    args = ["tail"]
    p3 = subprocess.Popen(args, stdin=p2.stdout, stdout=subprocess.PIPE)

    # Communicate with the last subprocess and get the output
    _ = p3.communicate()[0]  # output str

    # Assert that the return code of the tail command is 0, indicating success
    assert p3.returncode == 0


def test_cli_remove_clip_sortp_pipe():
    """Test the CLI command for removing clipped reads and sorting using pipes.

    This test specifically checks the integration of `samtools` and `bioat`
    through a pipeline approach, ensuring that clipped reads are removed
    correctly from the BAM file specified by `BAM_SORTP_FILE`.

    Args:
        None

    Raises:
        pytest.SkipException: if the bamnostic backend is used or if
        `samtools` is not installed.
    """
    if __BAM_PARSER_BACKEND__ == "bamnostic":
        pytest.skip("bamnostic backend does not support remove_clip")
    if not check_cmd("samtools"):
        pytest.skip("samtools is not installed")

    # Test the pipe functionality using `samtools` and `bioat`
    args = ["samtools", "view", "-h", BAM_SORTP_FILE]
    p1 = subprocess.Popen(args, stdout=subprocess.PIPE)

    args = ["bioat", "bam", "remove_clip"]
    p2 = subprocess.Popen(args, stdin=p1.stdout, stdout=subprocess.PIPE)

    args = ["tail"]
    p3 = subprocess.Popen(args, stdin=p2.stdout, stdout=subprocess.PIPE)

    # Capture the output, which should be the processed BAM data
    _ = p3.communicate()[0]  # Output as a string

    # Check that the final pipe returned successfully
    assert p3.returncode == 0
