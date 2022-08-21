#  dag_file.py - this file is part of the gwadaptive_scattering package.
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


class DagFile(object):
    """Class to create a .dag file to be
    submitted with `condor_submit_dag`.

    Parameters
    ----------
    name : str
        name of the .dag file
    """

    def __init__(self, name):
        if name[-4:] != ".dag":
            raise ValueError("Wrong file extension.")
        self.name = name
        self.dag_text = []

    def add_blank_line(self):
        """Add blank line to .dag file.
        """
        self.dag_text.append("")

    def add_job(self, job_number, sub_file, retry=1, args=None):
        """Add job block to .dag file.

        Parameters
        ----------
        job_number : int
            job number
        sub_file : str
            path to the .sub file
        retry : int, optional
            how many times retry the job if it fails (default : 1)
        args : dict, optional
            variables to be assigned in the .sub file (default : None)
        """
        self.dag_text.append("JOB {:d} {}".format(job_number, sub_file))
        self.dag_text.append("RETRY {} {}".format(job_number, retry))
        if args is not None:
            args_list = ["{}=\"{}\"".format(key, args[key]) for key in args.keys()]
            args_str = " ".join(args_list)
            self.dag_text.append("VARS {} {}".format(job_number, args_str))
        self.add_blank_line()

    def add_post_script(self, job_number, script, args=None):
        """Add post block to .dag file.

        Parameters
        ----------
        job_number : int
            job number
        script : str
            path to the script file
        args : list[str], optional
            `script` arguments (default : [])
        """
        if args is None:
            args = []

        args_list = " ".join(args)
        self.dag_text.append("SCRIPT POST {} {} {}".format(job_number, script, args_list))
        self.add_blank_line()

    def set_parent(self, parents, children):
        """List parent jobs to be followed by children jobs.

        Parameters
        ----------
        parents : list
            parent jobs
        children : list
            child jobs
        """
        parent_list = " ".join([str(p) for p in parents])
        child_list = " ".join([str(c) for c in children])
        self.dag_text.append("PARENT {} CHILD {}".format(parent_list, child_list))
        self.add_blank_line()

    def save(self):
        """Save .dag file to `name`.
        """
        dag_content = "\n".join(self.dag_text)
        f = open(self.name, "w")
        f.write(dag_content)
        f.close()
