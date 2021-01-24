#  asr_html.py - this file is part of the asr pagkage,
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
import sys
import argparse
import glob
import re
from .utils import file_utils
from datetime import datetime
from astropy.time import Time
from .html import scattered_light_page as slp
from .common import defines


SCRIPT_NAME = "./scattered_light_tool.py"
PLOTS_SCRIPT_NAME = "./scattered_light_plots.py"
HTML_SCRIPT_NAME = "./scattered_light_html.py"
SUMMARY_IMFS = 2
OMEGAGRAM_THR = 0.5

def generate_html(res_path):    
    res_folders = []
    for folder in os.listdir(res_path):
        curr_dir = os.path.join(res_path, folder)
        if os.path.isdir(curr_dir):
            if file_utils.yml_exists(curr_dir):
                res_folders.append(folder)
    res_folders.sort()
    
    # title
    curr_date = datetime.today().strftime('%d-%m-%Y')
    title = "Scattered light analysis ({})".format(curr_date)
    page = slp.ScatteredLightPage(title=title)
    
    # code
    page.openDiv(id_="page-command-line-section")
    page.addParagraph("This page can be reproduced with the following command line:", class_="mb-2")
    page.openDiv(id_="page-command-line")
    code = [HTML_SCRIPT_NAME]
    code.append("--ipath {}".format(res_path))
    page.addCommandLine(" ".join(code))
    page.closeDiv()
    page.closeDiv()
    
    # info
    page.addSection("Info", id_="info-section")
    page.openDiv(id_="info-list")
    
    if len(res_folders) > 0:
        tc_name = file_utils.load_yml(os.path.join(res_path, res_folders[0]))[defines.PARAMS_SECT_KEY][defines.TARGET_CH_KEY]
        chl_name = file_utils.load_yml(os.path.join(res_path, res_folders[0]))[defines.PARAMS_SECT_KEY][defines.CH_LIST_KEY]
    else:
        tc_name = ""
        chl_name = ""
    
    info_dict = {"Target channel": tc_name,
                 "Flags": "",
                 "GPS list": "",
                 "Channels list": page.getFormattedLink(os.path.basename(chl_name), **{"href": chl_name,
                                                                                       "download": chl_name})}
    page.addBulletList(info_dict)
    page.closeDiv()

    # results   
    page.addSection("Results", id_="results-section") 
    page.openDiv(id_="results")
    for gps_folder in res_folders:
        gps_path = os.path.join(res_path, gps_folder)
        res_file = file_utils.load_yml(gps_path)
        parameters = res_file[defines.PARAMS_SECT_KEY]
        imfs_data = {}
        above_thr = False
        for i in range(1, SUMMARY_IMFS + 1):
            try:
                imfs_data[i] = {}
                imfs_data[i]["Culprit"] = res_file[defines.CORR_SECT_KEY][i - 1][defines.CHANNEL_KEY]
                #imfs_data[i]["Correlation"] = "{:.4f}".format(res_file[defines.CORR_SECT_KEY][i - 1][defines.CORR_KEY])
                imfs_data[i]["Mean frequency"] = "{:.4f} Hz".format(res_file[defines.CORR_SECT_KEY][i - 1][defines.MEAN_FREQ_KEY])
                imfs_data[i]["omegagram"] = os.path.exists(os.path.join(gps_path, "imf_{:d}_omegagram.png".format(i)))
                imfs_data[i]["combo"] = ""
                
                regex = "[_+]{:d}[_+]".format(i)
                for cf in glob.glob(os.path.join(gps_path, "combo_imf_*_culprit.png")):
                    if re.search(regex, cf):
                        imfs_data[i]["combo"] = cf
                if res_file[defines.CORR_SECT_KEY][i - 1][defines.CORR_KEY] >= 0.5:
                    above_thr = True
            except:
                continue
        
        gps_start = gps_folder.split("_")[0]
        gps_end = gps_folder.split("_")[1]
        res_id = "{}-{}".format(gps_start, gps_end)
        t1 = Time(gps_start, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        gps_date = "{} UTC (GPS: {})\n".format(t2, gps_start)
        
        if above_thr:
            color_key = "danger"
        else:
            color_key = "warning"
        
        page.openDiv(class_="card border-{} mb-1 shadow-sm".format(color_key))
        
        # header
        page.openDiv(class_="card-header text-white bg-{}".format(color_key))
        page.addLink("{}".format(gps_date),
                     **{"class": "btn card-link cis-link",
                        "data-toggle": "collapse",
                        "data-target": "#{}".format(res_id)})
        page.closeDiv()
        
        # body
        page.openDiv(id_=res_id, class_="collapse")#,
                     #**{"data-parent": "#{}".format(parent_id)})
        page.openDiv(class_="card-body")
        page.parametersTable(parameters, int(gps_start), int(gps_end))
        
        page.openDiv(id_="command-line-{}".format(res_id))
        page.addParagraph("This analysis can be reproduced with the following command line:", class_="mb-2")
        code = [SCRIPT_NAME]
        for key in parameters.keys():
            code.append("--" + key)
            code.append(str(parameters[key]))
        page.addCommandLine(" ".join(code))
        page.closeDiv()
        
        page.openDiv(id_="plots-command-line-{}".format(res_id))
        page.addParagraph("Plots can be reproduced with the following command line:", class_="mb-2")
        code = [PLOTS_SCRIPT_NAME]
        code.append("--ipath {}".format(gps_path))
        code.append("--imfs_to_plot {}".format(",".join([str(i) for i in range(1, SUMMARY_IMFS + 1)])))
        code.append("--omegagram_thr {:.2f}".format(OMEGAGRAM_THR))
        page.addCommandLine(" ".join(code))
        page.closeDiv()
        
        for i in range(1, SUMMARY_IMFS + 1):
            if i in imfs_data.keys():
                page.openDiv(id_="imf-{}-{}".format(i, res_id))
                page.addSubsubsection("Imf {}".format(i), id_="imf-{}-{}-section".format(i, res_id))
                
                omegagram = imfs_data[i].pop("omegagram")
                combo_file = imfs_data[i].pop("combo")
                
                page.openDiv(id_="imf-{}-{}-info".format(i, res_id))
                page.addBulletList(imfs_data[i])
                page.closeDiv()
            
                imf_plot_name = os.path.join(gps_path, "imf_{}_culprit.png".format(i))
                page.openDiv(id_="imf-{}-{}-plot".format(i, res_id))
                page.addPlot(imf_plot_name, "imf-{}-{}".format(i, res_id))
                page.closeDiv()
                
                if combo_file != "":
                    page.openDiv(id_="combo-{}-{}-plot".format(i, res_id))
                    page.addPlot(combo_file, "combo-{}-{}".format(i, res_id))
                    page.closeDiv()
                
                if omegagram:
                    omegagram_plot_name = os.path.join(gps_path, "imf_{}_omegagram.png".format(i))
                    page.openDiv(id_="omegagram-{}-{}-plot".format(i, res_id))
                    page.addPlot(omegagram_plot_name, "omegagram-{}-{}".format(i, res_id))
                    page.closeDiv()
            
                page.closeDiv()
            
        page.closeDiv()
        page.closeDiv()
        page.closeDiv()
        
    page.closeDiv()
    
    # summary plots
    if os.path.exists(os.path.join(res_path, "comparison")):
        page.addSection("Summary", id_="summary-section")
        
        # command line
        page.openDiv(id_="command-line-summary")
        page.addParagraph("This summary can be reproduced with the following command line:", class_="mb-2")
        
        code = [PLOTS_SCRIPT_NAME]
        code.append("--ipath {}".format(res_path))
        code.append("--imfs_to_plot {}".format(",".join([str(i) for i in range(1, SUMMARY_IMFS + 1)])))
        code.append("--comparison True")
        page.addCommandLine(" ".join(code))
        
        page.closeDiv()
        
        # imfs
        for i in range(1, SUMMARY_IMFS + 1):
            page.openDiv(id_="imf-{}-summary".format(i))
            page.addSubsubsection("Imf {}".format(i), id_="imf-{}-section-summary".format(i))
            summary_name = os.path.join(res_path, "comparison", "imf_{}_summary_{}.png".format(i, tc_name))
            corr_summary_name = os.path.join(res_path, "comparison", "imf_{}_corr_summary_{}.png".format(i, tc_name))
            page.addPlot(summary_name, "imf-{}-summary".format(i))
            page.addPlot(corr_summary_name, "imf-{}-corr-summary".format(i))
            page.closeDiv()
    
    curr_date = datetime.today().strftime('%Y%m%d')
    html_file = os.path.join(res_path, "scattered_light_{}.html".format(curr_date))
    page.savePage(html_file)
    
