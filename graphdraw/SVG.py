import os
import pdb
from graphdraw.utils import hex_to_rgb, style_to_str

COLORS = {"Red": "FF0000",
          "White": "FFFFFF",
          "Cyan": "00FFFF",
          "Silver": "C0C0C0",
          "Blue": "0000FF",
          "Gray": "808080",
          "Grey": "808080",
          "DarkBlue": "00008B",
          "Black": "000000",
          "LightBlue": "ADD8E6",
          "Orange": "FFA500",
          "Purple": "800080",
          "Brown": "A52A2A",
          "Yellow": "FFFF00",
          "Maroon": "800000",
          "Lime": "00FF00",
          "Green": "008000",
          "Magenta": "FF00FF",
          "Olive": "808000",
          "Pink": "FFC0CB",
          "Aquamarine": "7FFFD4"}


def check_style(style=None):
    """
    check the style given for drawing a line
    """
    if not style:  # no style given, create one
        style = dict()
        style["stroke"] = "Black"
        style["stroke-width"] = 3
        style["stroke-opacity"] = 1
    elif not isinstance(style, dict):
        raise ValueError(
            "If you want to add style, it has to be a dict of name of style and value, e.g. 'stroke':'#000000'")
    else:
        if not style['stroke'] in COLORS:
            style['stroke'] = COLORS["Black"]
    return style


class SVG:
    def __init__(self, height=300, width=600, legend=10):
        self.height = height
        self.width = width
        self.legend = legend
        self.stack = []
        self.init_stack()

    def init_stack(self):
        """
        Initialise the SVG stack
        """
        self.stack.append(f'<svg height="{self.height}" width="{self.width}" style="background-color:white;">')
        self.stack.append('<rect width="100%" height="100%" fill="white"/>')

    def add_ref(self, reference):
        """
        Add the reference line to the SVG stack
        reference is a dictionary with keys start, end, length, start_px, end_px, color, thick, chromosome
        """
        color = hex_to_rgb(COLORS[reference['color']])
        y1 = y2 = int(self.height * 0.07)  # top margin for reference line (7%)
        self.stack.append(
            f'<line x1="{reference["start_px"]}" y1="{y1}" x2="{reference["end_px"]}" y2="{y2}" style="stroke:rgb{color};stroke-width:{reference["thick"]}"/>')

        y = int(self.height * 0.04)  # top margin for text
        self.stack.append(f'<text x="{reference["start_px"]}" y="{y}">Start: {reference["start"]}</text>')

        x = int(self.width * 0.15)  # end margins
        self.stack.append(f'<text x="{reference["end_px"]}" y="{y}" text-anchor="end">End: {reference["end"]}</text>')

        mid_point = int((reference["end_px"] - reference["start_px"]) / 2)
        self.stack.append(f'<text x="{mid_point}" y="{y}" text-anchor="middle">{reference["chromosome"]}</text>')

    def add_several_refs(self, references):
        """
        Add the reference line to the SVG stack
        reference is a dictionary with keys start, end, length, start_px, end_px, color, thick, chromosome
        """
        color = hex_to_rgb(COLORS[reference['color']])
        y1 = y2 = int(self.height * 0.07)  # top margin for reference line (7%)
        self.stack.append(
            f'<line x1="{reference["start_px"]}" y1="{y1}" x2="{reference["end_px"]}" y2="{y2}" style="stroke:rgb{color};stroke-width:{reference["thick"]}"/>')

        y = int(self.height * 0.04)  # top margin for text
        self.stack.append(f'<text x="{reference["start_px"]}" y="{y}">Start: {reference["start"]}</text>')

        x = int(self.width * 0.15)  # end margins
        self.stack.append(f'<text x="{reference["end_px"]}" y="{y}" text-anchor="end">End: {reference["end"]}</text>')

        mid_point = int((reference["end_px"] - reference["start_px"]) / 2)
        self.stack.append(f'<text x="{mid_point}" y="{y}" text-anchor="middle">{reference["chromosome"]}</text>')

    def add_txt(self, x, y, text, color="Black", size=5, rotation=-30, anchor="middle"):
        self.stack.append(
            f'<text fill="#{COLORS[color]}" transform="translate({x}, {y}) rotate({rotation})" font-size="{size}" stoke="none" text-anchor="{anchor}">{text}</text>')

    def add_line(self, x1, y1, x2, y2, style=None):
        style = check_style(style)

        style_str = style_to_str(style)
        # style_str = ";".join([f"{key}:{value}" for key, value in style.items()])
        self.stack.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" style="{style_str}"/>')

    def add_curve(self, x1, y1, x2, y2, len1, len2, dir1, dir2, style=None):
        """
        I am drawing edges as Bezier curves, the coordinates are for the start and end
        I am taking the length of the node and using that to calculate the pivot points for the curve
        I need to know from which direction the edge is to calculate the pivot points correctly
        len1 and len2 should be the length in pixels

        example of Bezier curve
        <path d="M x1 y1 C p_x1 p_y1, p_x2 p_y2, x2 y2" stroke="black" fill="transparent"/>

        <path d="M 10 10 C 50 10, 50 100, 100 100" fill='transparent' style="stroke:rgb(0,0,0);stroke-width:0.7;stroke-opacity:1"/>
        """
        style = check_style(style)

        p_y1 = y1
        p_y2 = y2

        # calculating pivot point for the curve based on the length of the node
        if dir1 == 0:  # from start
            p_x1 = x1 - int(len1 * 0.15)
        else:  # from end
            p_x1 = x1 + int(len1 * 0.15)

        if dir2 == 0:  # to start
            p_x2 = x2 - int(len2 * 0.15)
        else:  # to end
            p_x2 = x2 + int(len2 * 0.15)

        style_str = style_to_str(style)
        self.stack.append(
            f'<path d="M {x1} {y1} C {p_x1} {p_y1}, {p_x2} {p_y2}, {x2} {y2}" fill="transparent" style="{style_str}"/>')

    def add_arc(self):
        """
        adds an arc to the SVG stack
        this is an example of some arcs
        The first M y x
        is the start position of the arc, the A starts the arc line
        with the first two numbers being the radius, if same then circle, different then elips
        the following 3 numbers can be kept as 0 0 1, for 1 being clockwise or counterclockwise
        the last two numbers are the end point of the arc.
        Good website to understand them https://www.nan.fyi/svg-paths/arcs

        <svg width="100" height="100" xmlns="http://www.w3.org/2000/svg">
        <path d="M 20 50
               A 30 30 0 0 1 50 20
        " stroke="black" fill="transparent" stroke-width="2" fill-opacity="0.5"/>

        <path d="M 50 20
               A 30 30 0 0 1 80 50
        " stroke="red" fill="transparent" stroke-width="2" fill-opacity="0.5"/>

        <path d="M 80 50
               A 30 30 0 0 1 50 80
        " stroke="#000000" fill="transparent" stroke-width="2" fill-opacity="0.5"/>

        <path d="M 50 80
               A 30 30 0 0 1 20 50
        " stroke="blue" fill="transparent" stroke-width="2" fill-opacity="0.5"/>

        <path d="M 25 40
               A 25 25 0 0 1 40 25
        " stroke="blue" fill="transparent" stroke-width="2" fill-opacity="0.5"/>
        </svg>
        """
        style = check_style()
        pass

    def out_svg(self, out_file="graph.svg"):
        with open(out_file, 'w') as outfile:
            for c in self.stack:
                outfile.write(c + "\n")
            outfile.write("</svg>")
