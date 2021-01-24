#  scattered_light_page.py - this file is part of the asr pagkage,
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

from gwdetchar.io import html as htmlio
import os
from ..common import defines


CSS_FILES = ["https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.1/css/fontawesome.min.css",
             "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.1/css/solid.min.css",
             "https://cdn.jsdelivr.net/npm/gwbootstrap@1.3.1/lib/gwbootstrap.min.css"]

#JS_FILES = ["https://code.jquery.com/jquery-3.5.1.min.js",
#            "https://code.jquery.com/ui/1.12.1/jquery-ui.min.js",
#            "https://cdnjs.cloudflare.com/ajax/libs/moment.js/2.24.0/moment.min.js",
#            "https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/js/bootstrap.bundle.min.js",
#            "https://cdnjs.cloudflare.com/ajax/libs/fancybox/3.5.7/jquery.fancybox.min.js",
#            "https://cdnjs.cloudflare.com/ajax/libs/bootstrap-datepicker/1.9.0/js/bootstrap-datepicker.min.js",
#            "https://cdn.jsdelivr.net/npm/gwbootstrap@1.3.1/lib/gwbootstrap-extra.min.js"]

JS_FILES = ["https://code.jquery.com/jquery-3.5.1.min.js",
            "https://cdnjs.cloudflare.com/ajax/libs/jquery.lazy/1.7.11/jquery.lazy.min.js",
            "https://cdn.jsdelivr.net/npm/bootstrap@4.5.3/dist/js/bootstrap.bundle.min.js",
            "https://cdnjs.cloudflare.com/ajax/libs/fancybox/3.5.7/jquery.fancybox.min.js",
            "https://cdn.jsdelivr.net/npm/gwbootstrap@1.3.1/lib/gwbootstrap.min.js",
            "https://maxcdn.bootstrapcdn.com/bootstrap/4.1.0/js/bootstrap.min.js"]


class ScatteredLightPage(object):
    
    def __init__(self, title=""):
        self.page = htmlio.new_bootstrap_page(title=title, script=JS_FILES, css=CSS_FILES)#base=path, path=path,
        if title != "":
            self.page.h1(title, class_="mt-2 mb-2")
            
    
    def addParagraph(self, text, **kwargs):
        self.page.p(text, **kwargs)
        
        
    def addLink(self, text, **kwargs):
        self.page.a(text, **kwargs)
    

    def getFormattedLink(self, text, **kwargs):
        link = "<a"
        for k in kwargs.keys():
            link += " {}=\"{}\"".format(k, kwargs[k])
        link += ">{}</a>".format(text)
        
        return link
    

    def openDiv(self, **kwargs):
        self.page.div(**kwargs)
        
        
    def closeDiv(self):
        self.page.div.close()
        
        
    def addSection(self, title, **kwargs):
        self.addParagraph("<strong>" + title + "</strong>", **{"style": "font-size:180%;",
                                                               "class": "mt-2 mb-2"})
        
        
    def addSubsection(self, title, **kwargs):
        self.addParagraph("<strong>" + title + "</strong>", **{"style": "font-size:150%;",
                                                               "class": "mt-2 mb-2"})
        
        
    def addSubsubsection(self, title, **kwargs):
        self.addParagraph("<strong>" + title + "</strong>", **{"style": "font-size:120%;",
                                                               "class": "mt-2 mb-2"})
        
        
    def addBulletList(self, items):
        items_list = ("<strong>" + str(k) + "</strong>: " + str(items[k]) for k in items.keys())
        self.page.ul()
        self.page.li(items_list)
        self.page.ul.close()
        
        
    def parametersTable(self, parameters, start, end):
        content = []
        content.append(("GPS (<code>--{}</code>)".format(defines.GPS_KEY),
                        parameters[defines.GPS_KEY]))
        content.append(("Target channel (<code>--{}</code>)".format(defines.TARGET_CH_KEY),
                        parameters[defines.TARGET_CH_KEY]))
        content.append(("Channels list (<code>--{}</code>)".format(defines.CH_LIST_KEY),
                        parameters[defines.CH_LIST_KEY]))
        content.append(("Output path (<code>--{}</code>)".format(defines.OUT_PATH_KEY),
                        parameters[defines.OUT_PATH_KEY]))
        content.append(("Sampling frequency (<code>--{}</code>)".format(defines.SAMP_FREQ_KEY),
                        str(parameters[defines.SAMP_FREQ_KEY])))
        content.append(("Lowpass frequency (<code>--{}</code>)".format(defines.LOWPASS_FREQ_KEY),
                        str(parameters[defines.LOWPASS_FREQ_KEY])))
        content.append(("Predictor's harmonics (<code>--{}</code>)".format(defines.SCATTERING_KEY),
                        str(parameters[defines.SCATTERING_KEY])))
        content.append(("Smoothing window (<code>--{}</code>)".format(defines.SMOOTH_WIN_KEY),
                        str(parameters[defines.SMOOTH_WIN_KEY])))
        
        self.addSubsection("Parameters", id_="parameters-section-{}-{}".format(start, end))
        self.openDiv(class_="row")
        self.openDiv(class_="col-md-9 col-sm-12")
        self.page.add(htmlio.parameter_table(content, start, end))
        self.closeDiv()
        self.closeDiv()
        
        
    def addCommandLine(self, code):
        self.openDiv(class_="highlight", style_="background: #f8f8f8")
        self.page.pre(**{"style": "line-height: 125%;"})
        self.page.span(code)
        self.page.pre.close()
        self.closeDiv()
        
        
    def addPlot(self, plot_name, plot_id):
        aparams = {
            "class_": "fancybox",
            "target": "_blank",
            "data-fancybox": "gallery",
            "data-fancybox-group": "images"
        }
        self.page.a(href=plot_name, id_="a-{}".format(plot_id), **aparams)
        self.page.img(id_="img-{}".format(plot_id),
                      **{"src": plot_name, "alt": os.path.basename(plot_name),
                         "height": "400", "width": "500"})
        self.page.a.close()
    
    
    def savePage(self, path):
        htmlio.close_page(self.page, path)

