#  daily_correlations_page.py - this file is part of the gwadaptive_scattering package.
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
import sys
import re
import glob
from astropy.time import Time
from ..html import html_builder as hb
from ..utils import file_utils, signal_utils
from ..common import defines


PIPELINE_SCRIPT_NAME = "./gwas"
SUMMARY_IMFS = 1
COLOR_THRESHOLD_MIN = 0.5
COLOR_THRESHOLD = 0.7
PLOT_EXT = "png"


def generate_web_page(res_path, date, tc_name, aux_ch, prepend_path=""):
    """Output page for daily correlations analysis.

    Parameters
    ----------
    res_path : str
        path containing analyses' results folders
    date : str
        date to be published on the page (year-month-day)
    tc_name : str
        target channel name
    aux_ch : str
        auxiliary channel name
    prepend_path : str, optional
        path to be prepended to the relative path of images and
        files in the current page (default : "")
    """
    res_folders = file_utils.get_results_folders(res_path, must_include=["*.{}".format(PLOT_EXT)])
    page = hb.HtmlBuilder(style="body { background-color: white; }")

    # info
    page.add_section(defines.INFO_SECTION)
    page.open_div(id_="info-list")

    info_dict = {
        defines.PAGE_DAY: date,
        defines.TARGET_CH_NAME: tc_name,
        defines.AUXILIARY_CH_NAME: aux_ch,
        defines.RES_TABLE: page.get_formatted_link(defines.SUMMARY_NAME,
                                                   href=os.path.join(prepend_path,
                                                                     defines.COMPARISON_FOLDER,
                                                                     defines.SUMMARY_NAME),
                                                   download=os.path.join(prepend_path,
                                                                         defines.COMPARISON_FOLDER,
                                                                         defines.SUMMARY_NAME)),
        defines.ENV_NAME: page.get_formatted_code(sys.prefix)
    }
    page.add_bullet_list(info_dict)
    page.close_div()

    # correlation plot
    if os.path.exists(os.path.join(res_path, defines.COMPARISON_FOLDER)):
        page.add_section("Correlation between {} and {}".format(tc_name, aux_ch))

        for i in range(1, SUMMARY_IMFS + 1):
            seismic_plot_name = os.path.join(defines.COMPARISON_FOLDER, file_utils.seismic_plot_name(PLOT_EXT))
            if os.path.exists(os.path.join(res_path, seismic_plot_name)):
                page.open_div(id_="imf-{:d}-summary".format(i))
                page.add_subsection("Imf {:d}".format(i))
                page.add_plot(os.path.join(prepend_path, seismic_plot_name), "imf-{:d}-seismic".format(i))
                page.close_div()

    hour_dict = {}
    for ih in range(24):
        hour_dict[ih] = []

    for gf in res_folders:
        yf = file_utils.YmlFile(gf)
        gps_event = yf.get_gps()
        t1 = Time(gps_event, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        h = int("{}".format(t2).split(" ")[1].split(":")[0])
        hour_dict[h].append(gf)

    # results
    page.add_section(defines.SINGLE_GPS_SECTION)

    page.open_div(id_="info-color")

    info_color_dict = {
        defines.WARNING_STR: "Correlation between {:.1f} and {:.1f}".format(COLOR_THRESHOLD_MIN, COLOR_THRESHOLD),
        defines.ALERT_STR: "Correlation greater than {:.1f}".format(COLOR_THRESHOLD)
    }
    page.add_bullet_list(info_color_dict)
    page.close_div()

    page.open_div(id_="results")

    for hh in hour_dict.keys():
        if len(hour_dict[hh]) > 0:
            date_h = "{} {:02d}:00:00".format(date, hh)

            # main card
            page.open_card(date_h, "primary", "main-card-{:d}".format(hh))

            for gps_path in hour_dict[hh]:
                res_file = file_utils.YmlFile(gps_path)
                parameters = [
                    ("Event {} (UTC)".format(defines.GPS_PARAM),
                     "{:d} (at {}, {:d} seconds total)".format(res_file.get_gps(), res_file.get_event_position(),
                                                               res_file.get_seconds())),
                    (defines.SAMP_FREQ_PARAM, "{:.3f} Hz".format(res_file.get_sampling_frequency())),
                    (defines.LOWPASS_FREQ_PARAM, "{:.3f} Hz".format(res_file.get_lowpass_frequency())),
                    ("Plots and info", res_file.get_output_path())
                ]
                imfs_data = {}
                above_thr = False
                for i in range(1, SUMMARY_IMFS + 1):
                    if res_file.get_imfs_count() >= i:
                        imfs_data[i] = {}
                        imfs_data[i][defines.TARGET_CH_NAME] = tc_name
                        imfs_data[i][defines.COMBO_STR] = ""
                        regex = "[_+]{:d}[_+]".format(i)
                        for cf in glob.glob(os.path.join(gps_path, file_utils.combo_plot_name(["*"], PLOT_EXT))):
                            if re.search(regex, cf):
                                imfs_data[i][defines.COMBO_STR] = cf
                                imfs_data[i][defines.CULPRIT_STR] = res_file.get_channel_of_combo_with_imf(i)
                                break

                        try:
                            if res_file.get_corr_of_combo_with_imf(i) >= COLOR_THRESHOLD:
                                above_thr = True
                        except:
                            pass

                gps_event = res_file.get_gps()
                res_id = str(gps_event)
                t1 = Time(gps_event, format="gps")
                t2 = Time(t1, format="iso", scale="utc")
                gps_date = "{} UTC (GPS: {:d})\n".format(t2, gps_event)

                if above_thr:
                    color_key = "danger"
                else:
                    color_key = "warning"

                # card
                page.open_card(gps_date, color_key, res_id)

                # parameters table
                gps_start, gps_end = signal_utils.get_gps_interval_extremes(gps_event, res_file.get_seconds(),
                                                                            res_file.get_event_position())
                page.parameters_table(parameters, int(gps_start), int(gps_end))

                # plots
                for i in range(1, SUMMARY_IMFS + 1):
                    if i in imfs_data.keys():
                        combo_file = imfs_data[i].pop(defines.COMBO_STR)

                        page.open_div(id_="imf-{:d}-{}-sect".format(i, res_id))
                        page.add_subsection("Combo Imf {:d}".format(i))
                        page.open_div(id_="imf-{:d}-{}-info".format(i, res_id))
                        page.add_bullet_list(imfs_data[i])
                        page.close_div()

                        if combo_file != "":
                            combo_file = os.path.join(res_id, os.path.split(combo_file)[1])
                            page.open_div(id_="combo-{:d}-{}-plot".format(i, res_id))
                            page.add_plot(os.path.join(prepend_path, combo_file), "combo-{:d}-{}".format(i, res_id))
                            page.close_div()

                        page.close_div()

                page.close_card()

            page.close_card()

    page.close_div()

    html_file = os.path.join(res_path, defines.PAGE_NAME)
    page.save_page(html_file)
