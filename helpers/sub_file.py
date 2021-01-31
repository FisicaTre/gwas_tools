#  sub_file.py - this file is part of the asr package,
#  also known as "adaptive scattering recognizer".
#  Copyright (C) 2020- Stefano Bianchi
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see <https://www.gnu.org/licenses/>.


import os


class SubFile(object):
    """Class to create a .sub file to be
    submitted with `condor_submit`.
    
    Parameters
    ----------
    name : str
        name of the .sub file
    universe : str
        universe (default : vanilla)
    description : str
        description string of the sub file (default : empty)
    """

    
    def __init__(self, name, universe="vanilla", description=""):
        if name[-4:] != ".sub":
            raise ValueError("Wrong file extension.")
        self.name = name
        self.sub_text = ["#!/usr/bin/env condor_submit"]
        if description != "":
            self.sub_text.append("#")
            self.sub_text.append("# {}".format(description))
        self.addBlankLine()
        self.sub_text.append("universe = {}".format(universe))


    def addBlankLine(self):
        """Add blank line to .sub file.
        """
        self.sub_text.append("")


    def addExecutable(self, exec_name):
        """Add `executable` to .sub file.
        
        Parameters
        ----------
        exec_name : str
            executable name
        """
        self.addBlankLine()
        self.sub_text.append("executable = {}".format(exec_name))


    def addArguments(self, args):
        """Add `arguments` to .sub file.
        
        Parameters
        ----------
        args : str
            arguments as a single string
        """
        self.addBlankLine()
        self.sub_text.append("arguments = \" {} \"".format(args))


    def addAccountingGroupInfo(self, group, user):
        """Add `accounting_group` and `accounting_group_user` to .sub file.
        
        Parameters
        ----------
        group : str
            accounting group
        user : str
            accounting group user
        """
        self.addBlankLine()
        self.sub_text.append("accounting_group = {}".format(group))
        self.sub_text.append("accounting_group_user = {}".format(user))


    def addSpecs(self, ncpu, memory):
        """Add `request_cpus` and `request_memory` to .sub file.
        
        Parameters
        ----------
        ncpu : int
            number of cpus
        memory : int
            requested memory in Mb
        """
        self.addBlankLine()
        self.sub_text.append("request_cpus = {:d}".format(ncpu))
        self.sub_text.append("request_memory = {:d}".format(memory))


    def addLogs(self, output, error, to_append=[]):
        """Add `output` and `error` to .sub file.
        
        Parameters
        ----------
        output : str
            path to the output folder
        error : str
            path to the error folder
        to_append : list of str
            str to be appended to output and error file names
        """
        out_name = ["out", "$(Process)"] + to_append
        err_name = ["err", "$(Process)"] + to_append
        out_path = os.path.join(output, ".".join(out_name))
        err_path = os.path.join(error, ".".join(err_name))
        self.addBlankLine()
        self.sub_text.append("output = {}".format(out_path))
        self.sub_text.append("error = {}".format(err_path))


    def save(self):
        """Save .sub file to `name`.
        """
        self.addBlankLine()
        self.sub_text.append("notification = never")
        self.sub_text.append("getenv = True")
        self.sub_text.append("queue 1")
        sub_content = "\n".join(self.sub_text)
        f = open(self.name, "w")
        f.write(sub_content)
        f.close()

