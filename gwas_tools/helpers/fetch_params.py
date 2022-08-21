#  fetch_params.py - this file is part of the gwadaptive_scattering package.
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


def get_gps_and_freq(glitch_type, gps1, gps2, ifo, ml_confidence=(0.9, 1.0),
                     snr=(12, 20), save_path=None):
    """Get GPS and peak frequency between `gps1`
    and `gps2` for glitches of type `glitch_type`.
    
    Parameters
    ----------
    glitch_type : str
        glitch string identifier
    gps1 : str
        starting GPS
    gps2 : str
        ending GPS
    ifo : str
        interferometer identifier
    ml_confidence : tuple[float], optional
        ml confidence limits (default : (0.9, 1.0))
    snr : tuple[int], optional
        SNR limits (default : (12, 20))
    save_path : str, optional
        where to save output in csv format (default : None)
        
    Returns
    -------
    numpy array
        GPS times of the events
    numpy array
        Glitches peak frequencies
    """
    ml_low = str(ml_confidence[0])
    ml_high = str(ml_confidence[1])
    snr_low = str(snr[0])
    snr_high = str(snr[1])
    try:
        glitches_list = EventTable.fetch("gravityspy", "glitches_v2d0",
                                         selection=["ml_label={}".format(glitch_type),
                                                    "{}<=ml_confidence<={}".format(ml_low, ml_high),
                                                    "{}<=snr<={}".format(snr_low, snr_high),
                                                    "ifo={}".format(ifo),
                                                    "{}<event_time<{}".format(gps1, gps2)],
                                         host="gravityspyplus.ciera.northwestern.edu",
                                         user="mla", passwd="gl1tch35Rb4d!")

        glitches_list = glitches_list.to_pandas()
        glitches_list.drop_duplicates("peak_time", keep=False, inplace=True)
        peak_freqs = np.array(glitches_list.peak_frequency.values, dtype=float)
        start_times = np.array(glitches_list.peak_time.values, dtype=float)

        if save_path is not None:
            glitches_list.to_csv(save_path, index=False)
    except:
        start_times = []
        peak_freqs = []

    return start_times, peak_freqs


def get_gps_sequence(start, end, step):
    """Get a sequence of GPS times.

    Parameters
    ----------
    start : int
        starting GPS
    end : int
        ending GPS
    step : int
        time between two consecutive GPS times

    Returns
    -------
    numpy array
        consecutive GPS times (every `step` seconds)
    """
    return np.arange(start, end, step, dtype=int)
