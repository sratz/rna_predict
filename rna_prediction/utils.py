# coding: utf-8

import errno
import os


def read_file_line_by_line(filename, skip_empty=True):
    """
    Yields lines in a file while stripping whitespace.

    :param filename: filename to read
    :param skip_empty: True if empty lines should be skipped
    :return: next line in file
    """
    with open(filename) as lines:
        for line in lines:
            line = line.rstrip()
            if line == "" and skip_empty:
                continue
            yield line


def mkdir_p(path):
    """
    Creates directories recursively and does not error when they already exist.

    :param path: directory to create
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
