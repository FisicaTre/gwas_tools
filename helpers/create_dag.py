#  create_dag.py - this file is part of the asr pagkage,
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

import numpy as np
import argparse
from ..common import defines


ap = argparse.ArgumentParser()
ap.add_argument("--param_file", required=True)
ap.add_argument("--sub_file", required=True)
ap.add_argument("--out_file", required=True)
ap.add_argument("--stop_line", required=False, default=None, type=int)
ap.add_argument("--seconds", required=False, default=2, type=int)
args = vars(ap.parse_args())
param_file = args["param_file"]
sub_file = args["sub_file"]
out_file = args["out_file"]
stop_line = args["stop_line"]
seconds = args["seconds"]

if stop_line is None:
    stop_line = 10 ** 9

f = open(param_file, "r")
gps_flow = [(int(np.floor(float(line.rstrip().split(" ")[0]))),
             int(np.ceil(float(line.rstrip().split(" ")[1]))))
            for i, line in enumerate(f.readlines())
            if line.strip() and i < stop_line]
f.close()

f = open(out_file, "w")
for i, (gps, flow) in enumerate(gps_flow):
    f.write("JOB {} {}\n".format(i + 1, sub_file))
    f.write("RETRY {} 1\n".format(i + 1))
    f.write("VARS {} GPS_start_end=\"{},{}\" freq_lowpass=\"{}\"\n\n".format(i + 1, gps - defines.EXTRA_SECONDS,
                                                                             gps + seconds + defines.EXTRA_SECONDS, flow))
f.write("PARENT 1 CHILD {}\n".format(" ".join([str(n) for n in np.arange(2, len(gps_flow) + 1, dtype=int)])))
f.close()
