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
LOCK_INFO_KEY = "locked"
CORR_SECT_KEY = "correlations"
CORR_2_SECT_KEY = "correlations_2nd_best"
IMF_KEY = "imf"
CHANNEL_KEY = "channel"
CORR_KEY = "corr"
MEAN_FREQ_KEY = "mean_freq"
COMBO_SECT_KEY = "combos"
SEISMIC_SECT_KEY = "seismic"
LIGO_SEISMIC_CHANNELS = ["L1:ISI-GND_STS_ETMX_X_BLRMS_100M_300M",
                         "L1:ISI-GND_STS_ETMY_Y_BLRMS_100M_300M",
                         "L1:ISI-GND_STS_ETMX_Z_BLRMS_100M_300M",
                         "L1:ISI-GND_STS_ETMX_X_BLRMS_30M_100M",
                         "L1:ISI-GND_STS_ETMX_Y_BLRMS_30M_100M",
                         "L1:ISI-GND_STS_ETMX_Z_BLRMS_30M_100M"]  # ,
                         # "L1:ISI-GND_STS_ETMX_X_DQ",
                         # "L1:ISI-GND_STS_ETMX_Y_DQ",
                         # "L1:ISI-GND_STS_ETMX_Z_DQ"]
HANFORD_SEISMIC_CHANNELS = ["H1:ISI-GND_STS_ETMX_X_BLRMS_100M_300M",
                            "H1:ISI-GND_STS_ETMY_Y_BLRMS_100M_300M",
                            "H1:ISI-GND_STS_ETMX_Z_BLRMS_100M_300M",
                            "H1:ISI-GND_STS_ETMX_X_BLRMS_30M_100M",
                            "H1:ISI-GND_STS_ETMX_Y_BLRMS_30M_100M",
                            "H1:ISI-GND_STS_ETMX_Z_BLRMS_30M_100M"]  # ,
                            # "H1:ISI-GND_STS_ETMX_X_DQ",
                            # "H1:ISI-GND_STS_ETMX_Y_DQ",
                            # "H1:ISI-GND_STS_ETMX_Z_DQ"]
VIRGO_SEISMIC_CHANNELS = ["V1:ENV_WEB_SEIS_N_50Hz_rms_0.1_1Hz",
                          "V1:ENV_WEB_SEIS_W_50Hz_rms_0.1_1Hz",
                          "V1:ENV_WEB_SEIS_V_50Hz_rms_0.1_1Hz",
                          "V1:ENV_WEB_SEIS_N_50Hz_rms_0.03_0.1Hz",
                          "V1:ENV_WEB_SEIS_W_50Hz_rms_0.03_0.1Hz",
                          "V1:ENV_WEB_SEIS_V_50Hz_rms_0.03_0.1Hz"]

# common parameters
COMPARISON_FOLDER = "comparison"
EXTRA_SECONDS = 1
EVENT_LOCATION = ["start", "center", "end"]
LOWP_FREQ_OPTS = ["average", "max"]
SUMMARY_NAME = "summary_table.csv"

# lock channels
LCK_CH_VIRGO = "V1:META_ITF_LOCK_index"  # "V1:DQ_META_ITF_Mode"
VIRGO_SCIENCE_MODE_THR = 170
LCK_CH_LIGO = "L1:DMT-ANALYSIS_READY:1"
LCK_CH_HANFORD = "H1:DMT-ANALYSIS_READY:1"


# summary table
def summary_table_culprit_column_of_imf(i):
    return "culprit_{:d}".format(i)


def summary_table_correlation_column_of_imf(i):
    return "corr_{:d}".format(i)


def summary_table_mean_frequency_column_of_imf(i):
    return "mean_freq_{:d}".format(i)


# summary plots
FREQ_BANDS = [0.03, 0.1, 0.3, 1, 3, 10]
CHAMBERS = {"BSC5": ["ETMY", "TMSY"], "BSC4": ["ETMX", "TMSX"], "BSC1": ["ITMY"], "BSC3": ["ITMX"], "HAM6": ["OM"],
            "HAM2": ["-MC", "-PR"], "HAM3": ["MC2", "-PR2"], "HAM5": ["-SR"], "HAM4": ["-SR2"], "HAM1": ["-RM"]}

# html
INFO_SECTION = "Info"
PAGE_DAY = "Page day"
ENV_NAME = "Environment"
TARGET_CH_NAME = "Target channel"
AUXILIARY_CH_NAME = "Auxiliary channel"
GPS_LIST_NAME = "GPS list"
CH_LIST_NAME = "Channels list"
RES_TABLE = "Results table"
SINGLE_GPS_SECTION = "Single GPS results"
RESULTS_SECTION = "Results"
SUCCESS_STR = "Green"
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
