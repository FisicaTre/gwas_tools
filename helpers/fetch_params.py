#  fetch_params.py - this file is part of the asr package,
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
from gwpy.table import EventTable


def get_gps_and_freq(glitch_type, gps1, gps2, ifo, save_path=None):
    """Get GPS and peak frequency between `gps1`
    and `gps2` for glitches of type `glitch_type`
    Parameters:
    -----------
    glitch_type : str
        glitch string identifier
    gps1 : str
        starting GPS
    gps2 : str
        ending GPS
    ifo : str
        interferometer identifier
    save_path : str
        where to save output in csv format (default : None)
    Returns:
    --------
    numpy array
        GPS times of the events
    numpy array
        Glitches peak frequencies
    """
    glitches_list = EventTable.fetch("gravityspy", "glitches", #_v2d0",
                                     selection=["ml_label={}".format(glitch_type),
                                                "0.9<=ml_confidence<=1.0",
                                                "12<=snr<=20", "ifo={}".format(ifo),
                                                "{}<event_time<{}".format(gps1, gps2)],
                                     host="gravityspyplus.ciera.northwestern.edu",
                                     user="mla", passwd="gl1tch35Rb4d!")
    
    glitches_list = glitches_list.to_pandas()
    glitches_list.drop_duplicates("event_time", inplace=True)
    peak_freqs = np.array(glitches_list.peak_frequency.values, dtype=float)
    start_times = np.array(glitches_list.event_time.values, dtype=float)
    
    if save_path is not None:
        glitches_list.to_csv(save_path, index=False)

    return start_times, peak_freqs
