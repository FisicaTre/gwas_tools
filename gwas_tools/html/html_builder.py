#  html_builder.py - this file is part of the gwadaptive_scattering package.
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
from MarkupPy import markup
import os
from getpass import getuser
import datetime
from pytz import reference


REPO_BASE = "https://github.com/stfbnc"
PACKAGE_NAME = "gwadaptive_scattering"

CSS_FILES = ["https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.1/css/fontawesome.min.css",
             "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.1/css/solid.min.css",
             "https://cdn.jsdelivr.net/npm/gwbootstrap@1.3.1/lib/gwbootstrap.min.css"]

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
    kwargs : dict{str}
        general page attributes
    """

    def __init__(self, title="", **kwargs):
        self.page = htmlio.new_bootstrap_page(title=title, script=JS_FILES, css=CSS_FILES, **kwargs)
        if title != "":
            self.page.h1(title, class_="mt-2 mb-2")

    def add_paragraph(self, text, **kwargs):
        """Add <p>.
        
        Parameters
        ----------
        text : str
            <p> text
        kwargs : dict{str}
            <p> attributes
        """
        self.page.p(text, **kwargs)

    def add_link(self, text, **kwargs):
        """Add <a>.
        
        Parameters
        ----------
        text : str
            <a> text
        kwargs : dict{str}
            <a> attributes
        """
        self.page.a(text, **kwargs)

    def get_formatted_link(self, text, **kwargs):
        """Return <a> with formatted text.

        Parameters
        ----------
        text : str
            <a> text
        kwargs : dict{str}
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

    def get_formatted_code(self, text):
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

    def open_div(self, **kwargs):
        """Add <div>.
        
        Parameters
        ----------
        kwargs : dict{str}
            <div> attributes
        """
        self.page.div(**kwargs)

    def close_div(self):
        """Add </div>.
        """
        self.page.div.close()

    def add_section(self, title):
        """Add <p> with bold text magnified to 180%.
        
        Parameters
        ----------
        title : str
            <p> text
        """
        self.add_paragraph("<strong>" + title + "</strong>", **{"style": "font-size:180%;",
                                                                "class": "mt-2 mb-2"})

    def add_subsection(self, title):
        """Add <p> with bold text magnified to 150%.
        
        Parameters
        ----------                                            
        title : str
            <p> text
        """
        self.add_paragraph("<strong>" + title + "</strong>", **{"style": "font-size:150%;",
                                                                "class": "mt-2 mb-2"})

    def add_subsubsection(self, title):
        """Add <p> with bold text magnified to 120%.
        
        Parameters
        ----------                                            
        title : str
            <p> text
        """
        self.add_paragraph("<strong>" + title + "</strong>", **{"style": "font-size:120%;",
                                                                "class": "mt-2 mb-2"})

    def add_bullet_list(self, items):
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

    def parameters_table(self, parameters, start, end):
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
        self.add_subsubsection("Parameters")
        self.open_div(**{"class_": "row"})
        self.open_div(**{"class_": "col-md-9 col-sm-12"})
        self.page.add(htmlio.parameter_table(parameters, start, end))
        self.close_div()
        self.close_div()

    def add_command_line(self, code):
        """Add box with command line.
        
        Parameters
        ----------                                            
        code : str
            command line
        """
        self.open_div(**{"class_": "highlight", "style_": "background: #f8f8f8"})
        self.page.pre(**{"style": "line-height: 125%;"})
        self.page.span(code)
        self.page.pre.close()
        self.close_div()

    def add_command_line_block(self, code, description, block_id):
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
        self.open_div(**{"id_": "{}-section".format(block_id)})
        self.add_paragraph(description, **{"class_": "mb-2"})
        self.open_div(**{"id_": "{}".format(block_id)})
        self.add_command_line(code)
        self.close_div()
        self.close_div()

    def open_card(self, title, color, target_id):
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
        self.open_div(**{"class_": "card border-{} mb-1 shadow-sm".format(color)})

        # header
        self.open_div(**{"class_": "card-header text-white bg-{}".format(color)})
        self.add_link(title,
                      **{"class": "btn card-link cis-link",
                         "data-toggle": "collapse",
                         "data-target": "#{}".format(target_id)})
        self.close_div()

        # body
        self.open_div(**{"id_": target_id, "class_": "collapse"})
        self.open_div(**{"class_": "card-body"})

    def close_card(self):
        """Close a card element.
        """
        self.close_div()
        self.close_div()
        self.close_div()

    def add_plot(self, plot_name, plot_id):
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

    def add_footer(self):
        """Add footer with package and page creator information.
        """
        source = os.path.join(REPO_BASE, PACKAGE_NAME)
        issues = os.path.join(source, "issues")

        self.page.twotags.append("footer")
        markup.element("footer", case=self.page.case, parent=self.page)(class_="footer")
        self.open_div(**{"class_": "container"})
        self.open_div(**{"class_": "row"})
        self.open_div(**{"class_": "col-sm-3 icon-bar"})
        self.page.a(markup.oneliner.i('', class_='fas fa-code'), href=source,
                    title='View {} source code'.format(PACKAGE_NAME), target='_blank')
        self.page.a(markup.oneliner.i('', class_='fas fa-ticket-alt'), href=issues,
                    title='Report a problem', target='_blank')
        self.close_div()

        self.open_div(**{"class_": "col-sm-6"})
        now = datetime.datetime.now()
        tz = reference.LocalTimezone().tzname(now)
        date = now.strftime('%H:%M {} on %d %B %Y'.format(tz))
        self.add_paragraph("Created by {0} at {1}".format(getuser(), date))
        self.close_div()
        self.close_div()
        self.close_div()
        markup.element("footer", case=self.page.case, parent=self.page).close()

    def save_page(self, path, add_footer=False):
        """Save page.
        
        Parameters
        ----------                                            
        path : str
            save path
        add_footer : bool, optional
            add footer (default : False)
        """
        self.close_div()
        if add_footer:
            self.add_footer()
        if not self.page._full:
            self.page.body.close()
            self.page.html.close()
        with open(path, "w") as f:
            f.write(self.page())
