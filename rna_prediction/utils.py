'''
Created on Dec 19, 2014

@author: sebastian
'''

import errno
import os


def readFileLineByLine(filename, skipEmpty=True):
    with open(filename) as lines:
        for line in lines:
            line = line.rstrip()
            if line == "" and skipEmpty:
                continue
            yield line


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
