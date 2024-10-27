from bioat.exceptions import (
    BioatError,
    BioatFileFormatError,
    BioatFileNotCompleteError,
    BioatFileNotFoundError,
    BioatInvalidInputError,
    BioatInvalidOptionError,
    BioatMissingDependencyError,
    BioatRuntimeError,
    BioatValueError,
)


def test_bioat_error():
    error = BioatError("An error occurred")
    assert str(error) == "An error occurred"


def test_bioat_error_with_kwargs():
    error = BioatError("An error occurred", code=404)
    assert str(error) == "An error occurred; code=404"


def test_bioat_file_format_error():
    error = BioatFileFormatError("Invalid file format")
    assert str(error) == "Invalid file format"


def test_bioat_file_not_complete_error():
    error = BioatFileNotCompleteError("File is not complete")
    assert str(error) == "File is not complete"


def test_bioat_file_not_found_error():
    error = BioatFileNotFoundError("File not found")
    assert str(error) == "File not found"


def test_bioat_invalid_input_error():
    error = BioatInvalidInputError("Invalid input")
    assert str(error) == "Invalid input"


def test_bioat_invalid_option_error():
    error = BioatInvalidOptionError("Invalid option")
    assert str(error) == "Invalid option"


def test_bioat_missing_dependency_error():
    error = BioatMissingDependencyError("Missing dependency")
    assert str(error) == "Missing dependency"


def test_bioat_runtime_error():
    error = BioatRuntimeError("Runtime error occurred")
    assert str(error) == "Runtime error occurred"


def test_bioat_value_error():
    error = BioatValueError("Invalid value")
    assert str(error) == "Invalid value"
