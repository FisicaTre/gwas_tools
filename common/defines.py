#  defines.py - this file is part of the gwadaptive_scattering package.
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


# output.yml keys
PARAMS_SECT_KEY = "parameters"
GPS_KEY = "gps"
SECONDS_KEY = "seconds"
EVENT_KEY = "event"
TARGET_CH_KEY = "target_channel"
CH_LIST_KEY = "channels_list"
OUT_PATH_KEY = "opath"
SAMP_FREQ_KEY = "samp_freq"
LOWPASS_FREQ_KEY = "lowpass_freq"
SCATTERING_KEY = "scattering"
SMOOTH_WIN_KEY = "smooth_win"
CORR_SECT_KEY = "correlations"
CORR_2_SECT_KEY = "correlations_2nd_best"
IMF_KEY = "imf"
CHANNEL_KEY = "channel"
CORR_KEY = "corr"
MEAN_FREQ_KEY = "mean_freq"

# common parameters
COMPARISON_FOLDER = "comparison"
EXTRA_SECONDS = 1
EVENT_LOCATION = ["start", "center", "end"]
LOWP_FREQ_OPTS = ["average", "max"]

# lock channels
LCK_CH_VIRGO = "V1:DQ_META_ITF_Mode"

# html
INFO_SECTION = "info"
ENV_NAME = "Environment"
TARGET_CH_NAME = "Target channel"
AUXILIARY_CH_NAME = "Auxiliary channel"
GPS_LIST_NAME = "GPS list"
CH_LIST_NAME = "Channels list"
SINGLE_GPS_SECTION = "Single GPS results"
RESULTS_SECTION = "Results"
WARNING_STR = "Yellow"
ALERT_STR = "Red"
GPS_PARAM = "GPS"
SECONDS_PARAM = "Seconds"
EVENT_PARAM = "Event position"
TARGET_CH_PARAM = "Target channel"
CH_LIST_PARAM = "Channels list"
OUT_PATH_PARAM = "Output path"
SAMP_FREQ_PARAM = "Sampling frequency"
LOWPASS_FREQ_PARAM = "Lowpass frequency"
SCATTERING_PARAM = "Scattering factor"
SMOOTH_WIN_PARAM = "Smoothing window"
CULPRIT_STR = "Culprit"
MEAN_FREQ_STR = "Mean frequency"
OMEGAGRAM_STR = "omegagram"
COMBO_STR = "combo"
PAGE_NAME = "index.html"
