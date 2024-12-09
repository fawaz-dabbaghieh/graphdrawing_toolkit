import logging
import random
import gzip
import os
from colorsys import rgb_to_hsv
from math import sqrt


def read_fasta_gen(fasta_file_path):
    """
    A generator function that reads one read at a time
    Can be used for big FASTA files to not keep them in memory

    :param fasta_file_path: path to fasta file
    :yield: a tuple of sequence id and sequence
    """

    if not os.path.exists(fasta_file_path):
        logging.error("file {} does not exist".format(fasta_file_path))
        sys.exit()

    if fasta_file_path.endswith("gz"):
        fasta_file = gzip.open(fasta_file_path, "rt")
    else:
        fasta_file = open(fasta_file_path, "r")

    seqs = []
    seq_name = ""
    for line in fasta_file:
        line = line.strip()
        if not line:  # empty line
            continue

        if line.startswith(">"):
            if len(seqs) != 0:  # there was a sequence before
                yield seq_name, "".join(seqs)
                seq_name = line[1:]
                seqs = []
            else:
                seq_name = line[1:]
        else:
            seqs.append(line)

    # last sequence
    if seqs:
        yield seq_name, "".join(seqs)


def is_valid_hex(hex_string):
    """
    checks if the hex color that is represented in a string is valid or not
    """
    if not hex_string.startswith("#"):
        return False

    if len(hex_string) != 7:
        return False

    for c in hex_string[1:]:
        if not (('0' <= c <= '9') or ('a' <= c <= 'f') or ('A' <= c <= 'F')):
            return False
    return True


def hex_to_rgb(hex_color):
    if hex_color.startswith("#"):
        hex_color = hex_color[1:]
    rgb = tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return rgb


def rgb_to_hex(r, g, b):
    def check_bounds(x): return 0 <= x <= 255

    if check_bounds(r) and check_bounds(g) and check_bounds(b):
        return "#{:02x}{:02x}{:02x}".format(r, g, b)
    else:
        logging.warning("Invalid RGB value {}".format(hex_color))
        return "#000000"


def hex_distance(hex1, hex2):
    """
    Calculate the distance between two hex colors
    """
    (r1, g1, b1) = hex_to_rgb(hex1)
    (r2, g2, b2) = hex_to_rgb(hex2)
    return sqrt((r1 - r1)**2 + (g1 - g2)**2 + (b1 - b2)**2)


def rgb_picker(previous_hex, main_color):
    """
    Pick a random RGB color and return the hex color code.
    previous_hex: hex color that was previously chosen
    main_color: is an int with values 0, 1, 2 which represent the main color category red, green or blue
    """
    while True:
        # r = random.randint(0, 255)
        # g = random.randint(0, 255)
        # b = random.randint(0, 255)
        if main_color == 0:
            r = 255
            g = random.randint(0, 255)
            b = random.randint(0, 255)
        elif main_color == 1:
            r = random.randint(0, 255)
            g = 255
            b = random.randint(0, 255)
        else:
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = 255

        # darkness = (r*299 + g*587 + b*114)/1000
        # luma = 0.2126 * r + 0.7152 * g + 0.0722 * b
        h, s, v = rgb_to_hsv(r, g, b)
        # if darkness > 100 and luma > 100:
        if h >= 0.3 and s >= 0.15:
            if hex_distance(previous_hex, rgb_to_hex(r, g, b)) > 100:
                break
    return rgb_to_hex(r, g, b)


def style_to_str(style):
    style_str = []
    for key, value in style.items():
        style_str.append(f'{key}:{value}')
        # if key == "stroke":
        #     style_str.append(f'{key}:"{value}"')
        # else:
        #     style_str.append(f'{key}:{value}')
    return ";".join(style_str)
