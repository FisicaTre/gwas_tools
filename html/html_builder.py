#  html_builder.py - this file is part of the asr package,
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

# JS_FILES = ["https://code.jquery.com/jquery-3.5.1.min.js",
#             "https://code.jquery.com/ui/1.12.1/jquery-ui.min.js",
#             "https://cdnjs.cloudflare.com/ajax/libs/moment.js/2.24.0/moment.min.js",
#             "https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/js/bootstrap.bundle.min.js",
#             "https://cdnjs.cloudflare.com/ajax/libs/fancybox/3.5.7/jquery.fancybox.min.js",
#             "https://cdnjs.cloudflare.com/ajax/libs/bootstrap-datepicker/1.9.0/js/bootstrap-datepicker.min.js",
#             "https://cdn.jsdelivr.net/npm/gwbootstrap@1.3.1/lib/gwbootstrap-extra.min.js"]

JS_FILES = ["https://code.jquery.com/jquery-3.5.1.min.js",
            "https://cdnjs.cloudflare.com/ajax/libs/jquery.lazy/1.7.11/jquery.lazy.min.js",
            "https://cdn.jsdelivr.net/npm/bootstrap@4.5.3/dist/js/bootstrap.bundle.min.js",
            "https://cdnjs.cloudflare.com/ajax/libs/fancybox/3.5.7/jquery.fancybox.min.js",
            "https://cdn.jsdelivr.net/npm/gwbootstrap@1.3.1/lib/gwbootstrap.min.js",
            "https://maxcdn.bootstrapcdn.com/bootstrap/4.1.0/js/bootstrap.min.js"]


class HtmlBuilder(object):
    """Class to build HTML page components.
    
        Parameters
        ----------
        title : str
            page title
        kwargs : dict
            general page attributes
        """

    def __init__(self, title="", **kwargs):
        self.page = htmlio.new_bootstrap_page(title=title, script=JS_FILES, css=CSS_FILES, **kwargs)
        if title != "":
            self.page.h1(title, class_="mt-2 mb-2")

    def addParagraph(self, text, **kwargs):
        """Add <p>.
        
        Parameters
        ----------
        text : str
            <p> text
        kwargs : dict
            <p> attributes
        """
        self.page.p(text, **kwargs)

    def addLink(self, text, **kwargs):
        """Add <a>.
        
        Parameters
        ----------
        text : str
            <a> text
        kwargs : dict
            <a> attributes
        """
        self.page.a(text, **kwargs)

    def getFormattedLink(self, text, **kwargs):
        """Return <a> with formatted text.

        Parameters
        ----------
        text : str
            <a> text
        kwargs : dict
            <a> attributes

        Returns
        -------
        str
            <a kwargs>text</a>
        """
        link = "<a"
        for k in kwargs.keys():
            link += " {}=\"{}\"".format(k, kwargs[k])
        link += ">{}</a>".format(text)

        return link

    def getFormattedCode(self, text):
        """Return text in <code> tag.

        Parameters
        ----------
        text : str
            <code> text

        Returns
        -------
        str
            <code>text</code>
        """
        return "<code>{}</code>".format(text)

    def openDiv(self, **kwargs):
        """Add <div>.
        
        Parameters
        ----------
        kwargs : dict
            <div> attributes
        """
        self.page.div(**kwargs)

    def closeDiv(self):
        """Add </div>.
        """
        self.page.div.close()

    def addSection(self, title, **kwargs):
        """Add <p> with bold text magnified to 180%.
        
        Parameters
        ----------
        title : str
            <p> text
        kwargs : dict
            <p> attributes
        """
        self.addParagraph("<strong>" + title + "</strong>", **{"style": "font-size:180%;",
                                                               "class": "mt-2 mb-2"})

    def addSubsection(self, title, **kwargs):
        """Add <p> with bold text magnified to 150%.
        
        Parameters
        ----------                                            
        title : str
            <p> text
        kwargs : dict
            <p> attributes
        """
        self.addParagraph("<strong>" + title + "</strong>", **{"style": "font-size:150%;",
                                                               "class": "mt-2 mb-2"})

    def addSubsubsection(self, title, **kwargs):
        """Add <p> with bold text magnified to 120%.
        
        Parameters
        ----------                                            
        title : str
            <p> text
        kwargs : dict
            <p> attributes
        """
        self.addParagraph("<strong>" + title + "</strong>", **{"style": "font-size:120%;",
                                                               "class": "mt-2 mb-2"})

    def addBulletList(self, items):
        """Add bullet list.
        
        Parameters
        ----------                                            
        items : dict
            list items
        """
        items_list = ("<strong>" + str(k) + "</strong>: " + str(items[k]) for k in items.keys())
        self.page.ul()
        self.page.li(items_list)
        self.page.ul.close()

    def parametersTable(self, parameters, start, end):
        """Add parameters analysis table.
        
        Parameters
        ----------                        
        parameters : list[tuple]
            analysis parameters
        start : int
            start GPS
        end : int
            end GPS
        """
        self.addSubsection("Parameters")
        self.openDiv(**{"class_": "row"})
        self.openDiv(**{"class_": "col-md-9 col-sm-12"})
        self.page.add(htmlio.parameter_table(parameters, start, end))
        self.closeDiv()
        self.closeDiv()

    def addCommandLine(self, code):
        """Add box with command line.
        
        Parameters
        ----------                                            
        code : str
            command line
        """
        self.openDiv(**{"class_": "highlight", "style_": "background: #f8f8f8"})
        self.page.pre(**{"style": "line-height: 125%;"})
        self.page.span(code)
        self.page.pre.close()
        self.closeDiv()

    def addCommandLineBlock(self, code, description, block_id):
        """Add box with command line with a paragraph above.

        Parameters
        ----------
        code : str
            command line
        description : str
            description of the command line
        block_id : str
            div block id
        """
        self.openDiv(**{"id_": "{}-section".format(block_id)})
        self.addParagraph(description, **{"class_": "mb-2"})
        self.openDiv(**{"id_": "{}".format(block_id)})
        self.addCommandLine(" ".join(code))
        self.closeDiv()
        self.closeDiv()

    def openCard(self, title, color, target_id):
        """Open a card element.

        Parameters
        ----------
        title : str
            title shown on the card
        color : str
            card color
        target_id : str
            id to link body to header
        """
        self.openDiv(**{"class_": "card border-{} mb-1 shadow-sm".format(color)})

        # header
        self.openDiv(**{"class_": "card-header text-white bg-{}".format(color)})
        self.addLink(title,
                     **{"class": "btn card-link cis-link",
                        "data-toggle": "collapse",
                        "data-target": "#{}".format(target_id)})
        self.closeDiv()

        # body
        self.openDiv(**{"id_": target_id, "class_": "collapse"})
        self.openDiv(**{"class_": "card-body"})

    def closeCard(self):
        """Close a card element.
        """
        self.closeDiv()
        self.closeDiv()
        self.closeDiv()

    def addPlot(self, plot_name, plot_id):
        """Add plot.
        
        Parameters
        ----------                                            
        plot_name : str
            path to plot
        plot_id : str
            plot identifier
        """
        aparams = {
            "class_": "fancybox",
            "target": "_blank",
            "data-fancybox": "gallery",
            "data-fancybox-group": "images"
        }
        self.page.a(href=plot_name, id_="a-{}".format(plot_id), **aparams)
        self.page.img(id_="img-{}".format(plot_id),
                      **{"src": plot_name, "alt": os.path.basename(plot_name),
                         "style": "height: 50%; width: 50%; object-fit: contain"})
        self.page.a.close()

    def savePage(self, path):
        """Save page.
        
        Parameters
        ----------                                            
        path : str
            save path
        """
        htmlio.close_page(self.page, path)
