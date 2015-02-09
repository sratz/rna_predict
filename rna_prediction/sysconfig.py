"""
Created on Sep 22, 2014

@author: sebastian
"""

import re
import ConfigParser
import os
from os.path import expanduser

class SysConfig(object):
    SYSCONFIG_LOCATION = expanduser("~/.rna_predict")
    SYSCONFIG_FILE = SYSCONFIG_LOCATION + os.sep + "sysconfig"

    def __init__(self):
        """
        Load system configuration
        """
        #defaults
        self.rosetta_exe_path = ""
        self.rosetta_exe_suffix = ".linuxgccrelease"
        self.subprocess_buffsize = None

        self.loadSysConfig()

    def loadSysConfig(self):
        config = ConfigParser.RawConfigParser()
        config.read(SysConfig.SYSCONFIG_FILE)
        if config.has_section("rosetta"):
            if config.has_option("rosetta", "exe_path"):
                self.rosetta_exe_path = re.sub("/+$", "", config.get("rosetta", "exe_path")) + "/"
            if config.has_option("rosetta", "exe_suffix"):
                self.rosetta_exe_suffix = config.get("rosetta", "exe_suffix")
        if config.has_section("rna_predict"):
            if config.has_option("rna_predict", "subprocess_buffsize"):
                self.subprocess_buffsize = config.get("rna_predict", "subprocess_buffsize")

    def checkSysConfig(self):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        progs = [self.rosetta_exe_path + "rna_helix" + self.rosetta_exe_suffix]

        if self.subprocess_buffsize is not None:
            progs += ["stdbuf"]

        def isOk(prog):
            fpath, fname = os.path.split(prog)
            if fpath:
                if is_exe(prog):
                    return True
            else:
                for path in os.environ["PATH"].split(os.pathsep):
                    path = path.strip('"')
                    exe_file = os.path.join(path, prog)
                    if is_exe(exe_file):
                        return True

        fail = []
        success = []
        for prog in progs:
            if isOk(prog):
                success += [prog]
            else:
                fail += [prog]
        return {    "fail": fail, "success": success}

    def printSysConfig(self):
        print "System configuration:"
        for key, value in sorted(self.__dict__.items()):
            print "    %s: %s" %(key, "-" if value is None else value)
