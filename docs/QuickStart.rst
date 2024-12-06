QuickStart
==========

This is the QuickStart guide for the project.

Installation
------------

To install `bioat`, run the following commands:

Supported platforms: Linux / MacOS (intel & arm64) / Windows.

.. code-block:: shell

    # Install using pip
    pip install --upgrade bioat

    # Install using conda
    conda install -c conda-forge bioat

Usage
-----

As a Python library:

.. code-block:: python

    from bioat import about

    print(about)

Below are some example commands to use `bioat` as a CLI tool:

.. code-block:: shell

    # List commands
    bioat list

    # Check version
    bioat version

    # Check information about bioat
    bioat about

    # Example usage
    bioat bam remove_clip --help
    samtools view -h test_sorted_n.bam | bioat bam remove_clip | tail