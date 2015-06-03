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


def comma_separated_ranges_to_list(s):
    """
    Parses a string containing comma separated ranges to a list.

    Example:
    Turns ``1-3,10,20-22`` into ``[1, 2, 3, 10, 20, 21, 22]``

    :param s: comma separated ranges
    :return: list of ints
    """
    l = []
    ranges = s.split(",")
    for r in ranges:
        rs = r.split("-")
        if len(rs) > 2:
            raise ValueError("Invalid range: %s" % r)
        if len(rs) == 1:
            # single number
            l += [int(rs[0])]
        else:
            # start-end
            start = int(rs[0])
            end = int(rs[1])
            if start >= end:
                raise ValueError("Invalid range: %s" % r)
            l += [x for x in xrange(start, end + 1)]
    return l


def comma_separated_entries_to_dict(s, type_key, type_value):
    """
    Parses a string containing comma separated key-value paris (colon-separated) into a dict with fixed key and value types.

    Example:
    Turns ``foo:3,bar:4`` into ``{"foo": 3, "bar": 4}``

    :param s: input string
    :param type_key: type of the keys
    :param type_value: type of the values
    :return: dict
    """
    return {type_key(k): type_value(v) for k, v in map(lambda x: x.split(":"), s.split(","))}
