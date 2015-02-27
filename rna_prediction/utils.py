# coding: utf-8

import errno
import os


def read_file_line_by_line(filename, skip_empty=True):
    with open(filename) as lines:
        for line in lines:
            line = line.rstrip()
            if line == "" and skip_empty:
                continue
            yield line


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
