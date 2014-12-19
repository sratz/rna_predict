'''
Created on Dec 19, 2014

@author: sebastian
'''


def readFileLineByLine(filename, skipEmpty=True):
    with open(filename) as lines:
        for line in lines:
            line = line.rstrip()
            if line == "" and skipEmpty:
                continue
            yield line
