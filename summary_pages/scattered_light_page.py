#  scattered_light_page.py - this file is part of the gwadaptive_scattering package.
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
import glob
import re
import sys
from astropy.time import Time
from ..html import html_builder as hb
from ..utils import file_utils, signal_utils
from ..common import defines


PIPELINE_SCRIPT_NAME = "./gwas_glitch"
COLOR_THRESHOLD_MIN = 0.5
COLOR_THRESHOLD_MAX = 0.7
PLOT_EXT = "png"


def generate_web_page(res_path, date, tc_name, ch_list_file, gps_file, summary_imfs=2):
    """Output page for scattered light analysis.

    Parameters
    ----------
    res_path : str
        path containing analyses' results folders
    date : str
        date to be published on the page (year-month-day)
    tc_name : str
        target channel name
    ch_list_file : str
        channels list
    gps_file : str
        gps glitches list
    summary_imfs : int, optional
        imfs to show in the page
    """
    res_folders = file_utils.get_results_folders(res_path, must_include=["imf_*_culprit.{}".format(PLOT_EXT)])

    # title
    title = "Scattered light analysis ({})".format(date)
    page = hb.HtmlBuilder(title=title, style="body { background-color: white; }")

    # pipeline code
    code = [
        PIPELINE_SCRIPT_NAME,
        "--target_channel {}".format(tc_name),
        "--channels_list {}".format(ch_list_file),
        "--date {}".format("".join(date.split("-")))
    ]
    description = "This page can be reproduced with the following command line:"
    page.add_command_line_block(" ".join(code), description, "pipeline-command-line")

    # info
    page.add_section(defines.INFO_SECTION)
    page.open_div(id_="info-list")

    # copy channels file to local folder
    os.system("cp {} {}".format(ch_list_file, res_path))
    info_dict = {
        defines.ENV_NAME: page.get_formatted_code(sys.prefix),
        defines.TARGET_CH_NAME: tc_name,
        defines.GPS_LIST_NAME: page.get_formatted_link(os.path.basename(gps_file),
                                                       href=os.path.basename(gps_file),
                                                       download=os.path.basename(gps_file)),
        defines.CH_LIST_NAME: page.get_formatted_link(os.path.basename(ch_list_file),
                                                      href=os.path.basename(ch_list_file),
                                                      download=os.path.basename(ch_list_file))
    }
    page.add_bullet_list(info_dict)
    page.close_div()

    # summary plots
    if os.path.exists(os.path.join(res_path, defines.COMPARISON_FOLDER)):
        page.add_section("Summary")

        for i in range(1, summary_imfs + 1):
            freq_range_plot_name = os.path.join(defines.COMPARISON_FOLDER, file_utils.summary_freq_plot_name(PLOT_EXT))
            chamber_plot_name = os.path.join(defines.COMPARISON_FOLDER, file_utils.summary_chamber_plot_name(PLOT_EXT))
            if os.path.exists(os.path.join(res_path, freq_range_plot_name)) or os.path.exists(os.path.join(res_path, chamber_plot_name)):
                page.open_div(id_="imf-{:d}-summary".format(i))
                page.add_subsection("Imf {:d}".format(i))
                if os.path.exists(os.path.join(res_path, freq_range_plot_name)):
                    page.add_plot(freq_range_plot_name, "imf-{:d}-freq-range".format(i))
                if os.path.exists(os.path.join(res_path, chamber_plot_name)):
                    page.add_plot(chamber_plot_name, "imf-{:d}-chamber".format(i))
                page.close_div()

    # results
    page.add_section(defines.RESULTS_SECTION)
    page.open_div(id_="results")
    for gps_folder in res_folders:
        gps_path = os.path.join(res_path, gps_folder)
        res_file = file_utils.YmlFile(gps_path)
        parameters = [
            ("Event {} (UTC)".format(defines.GPS_PARAM),
             "{:d} (at {}, {:d} seconds total".format(res_file.get_gps(), res_file.get_event_position(),
                                                      res_file.get_seconds())),
            (defines.SAMP_FREQ_PARAM, "{:.3f} Hz".format(res_file.get_sampling_frequency())),
            (defines.LOWPASS_FREQ_PARAM, "{:.3f} Hz".format(res_file.get_lowpass_frequency())),
            ("Plots and info", res_file.get_output_path())
        ]
        imfs_data = {}
        above_thr_max = False
        above_thr_min = False
        for i in range(1, summary_imfs + 1):
            if res_file.get_imfs_count() >= i:
                imfs_data[i] = {}
                imfs_data[i][defines.TARGET_CH_NAME] = tc_name
                imfs_data[i][defines.CULPRIT_STR] = res_file.get_channel_of_imf(i)
                imfs_data[i][defines.MEAN_FREQ_STR] = "{:.4f} Hz".format(res_file.get_mean_freq_of_imf(i))
                imfs_data[i][defines.OMEGAGRAM_STR] = os.path.exists(os.path.join(gps_path,
                                                                                  file_utils.omegagram_plot_name(i, PLOT_EXT)))
                # imfs_data[i][defines.COMBO_STR] = ""
                # regex = "[_+]{:d}[_+]".format(i)
                # for cf in glob.glob(os.path.join(gps_path, file_utils.combo_plot_name(["*"], PLOT_EXT))):
                #    if re.search(regex, cf):
                #        imfs_data[i][defines.COMBO_STR] = cf
                #        break
                imf_i_corr = res_file.get_corr_of_imf(i)
                if imf_i_corr >= COLOR_THRESHOLD_MAX:
                    above_thr_max = True
                elif COLOR_THRESHOLD_MIN <= imf_i_corr < COLOR_THRESHOLD_MAX:
                    above_thr_min = True

        gps_event = res_file.get_gps()
        res_id = "{:d}".format(gps_event)
        t1 = Time(gps_event, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        gps_date = "{} UTC (GPS: {:d})\n".format(t2, gps_event)

        if above_thr_max:
            color_key = "danger"
        elif above_thr_min:
            color_key = "warning"
        else:
            color_key = "success"

        # card
        page.open_card(gps_date, color_key, res_id)

        # parameters table
        gps_start, gps_end = signal_utils.get_gps_interval_extremes(gps_event, res_file.get_seconds(),
                                                                    res_file.get_event_position())
        page.parameters_table(parameters, int(gps_start), int(gps_end))

        # plots
        for i in range(1, summary_imfs + 1):
            if i in imfs_data.keys():
                imf_plot_name = os.path.join(res_id, file_utils.culprit_plot_name(i, PLOT_EXT))
                # combo_file = imfs_data[i].pop(defines.COMBO_STR)
                omegagram = imfs_data[i].pop(defines.OMEGAGRAM_STR)

                page.open_div(id_="imf-{:d}-{}-sect".format(i, res_id))
                page.add_subsection("Imf {:d}".format(i))
                page.open_div(id_="imf-{:d}-{}-info".format(i, res_id))
                page.add_bullet_list(imfs_data[i])
                page.close_div()

                if os.path.exists(os.path.join(res_path, imf_plot_name)):
                    page.open_div(id_="imf-{:d}-{}-plot".format(i, res_id))
                    page.add_plot(imf_plot_name, "imf-{:d}-{}".format(i, res_id))
                    page.close_div()

                # if combo_file != "":
                #    combo_file = os.path.join(res_id, os.path.split(combo_file)[1])
                #    page.open_div(id_="combo-{:d}-{}-plot".format(i, res_id))
                #    page.add_plot(combo_file, "combo-{:d}-{}".format(i, res_id))
                #    page.close_div()

                if omegagram:
                    omegagram_plot_name = os.path.join(res_id, file_utils.omegagram_plot_name(i, PLOT_EXT))
                    page.open_div(id_="omegagram-{:d}-{}-plot".format(i, res_id))
                    page.add_plot(omegagram_plot_name, "omegagram-{:d}-{}".format(i, res_id))
                    page.close_div()

                page.close_div()

        page.close_card()

    page.close_div()

    html_file = os.path.join(res_path, defines.PAGE_NAME)
    page.save_page(html_file)
