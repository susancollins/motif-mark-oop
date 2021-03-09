#!/usr/bin/env python

import argparse
import re
import cairo
import math
import seaborn as sns


def get_args():
    '''to run script from the command line'''
    parser = argparse.ArgumentParser(description="find motifs and exons in fasta file, outputs a diagram of motif/exon locations")
    parser.add_argument("-f", "--fasta", help="input fasta file")
    parser.add_argument("-m", "--motifs", help="txt file of motifs")
    return parser.parse_args()

args = get_args()

input_file = args.fasta
motif_file = args.motifs

#create output filename from input file name
output_file = input_file.split(".")[0] + ".svg"

#create starting surface for pycairo
surface = cairo.SVGSurface(output_file, 1000, 1000)
context = cairo.Context(surface)

######CLASSES#######

class Gene:
    '''Draws gene as a black line'''
    
    def __init__(self, sequence):

        self.length = len(sequence)
        self.width = 1

    def draw(self, context, gene_number):
        y = HEIGHT_GENE_GROUP * gene_number + Y_OFFSET_GENE
        context.set_line_width(self.width)
        context.set_source_rgba(0, 0, 0)
        context.move_to(LEFT_MARGIN, y)
        context.line_to(LEFT_MARGIN + self.length, y)
        context.stroke()


class Exon:
    '''Draws exon as black box'''

    def __init__(self, sequence):

        self.width = 30
        self.sequence = sequence
        self.start = self.ExonFinder()[0]
        self.length = self.ExonFinder()[1]
    
    def ExonFinder(self):
        '''Takes a sequences as input and returns index of exon and exon length'''
        exon_location = []
        for n, i in enumerate(self.sequence):
            if i.isupper():
                exon_location.append(n)
        #creates tuple containing location of start of exon and exon length
        exon_info = (exon_location[0], len(exon_location))
        return exon_info

    def draw(self, context, gene_number):
        y = HEIGHT_GENE_GROUP * gene_number + Y_OFFSET_GENE
        context.set_line_width(self.width)
        context.set_source_rgba(0, 0, 0)
        context.move_to(LEFT_MARGIN + self.start, y)
        context.line_to(LEFT_MARGIN + self.start + self.length, y)
        context.stroke()


class Motif:
    '''Draws motifs as color-coded vertical lines and key'''

    def __init__(self, sequence):
        self.sequence = sequence
        self.locations = self.MotifFinder()

    def MotifFinder(self):
        '''Find start locations of motifs in a read and returns dictionary of start indexes'''
        motif_locations = {}
        motif_num = 0
        for motif in MOTIF_DICT:
            motif_num += 1
            motif_length = len(MOTIF_DICT[motif])
            locations = [i.start() for i in re.finditer(motif, self.sequence)]
            motif_locations.setdefault(MOTIF_DICT[motif], (motif_num, motif_length, locations))
        return motif_locations

    def draw(self, context, gene_number):
        y = HEIGHT_GENE_GROUP * gene_number + Y_OFFSET_GENE
        for motif in self.locations:
            y1, y2 = y+15, y-15
            context.set_line_width(self.locations[motif][1])
            motif_color = MOTIF_COLOR_PAL[self.locations[motif][0]-1]
            context.set_source_rgba(motif_color[0], motif_color[1], motif_color[2], 0.75)
            for j in self.locations[motif][2]:
                context.move_to(LEFT_MARGIN + j, y1)
                context.line_to(LEFT_MARGIN + j, y2)
                context.stroke()

    def draw_legend(self, context):
        #draw legend
        for motif in self.locations:
            y0 = 25*self.locations[motif][0]
            context.set_line_width(20)
            motif_color = MOTIF_COLOR_PAL[self.locations[motif][0]-1]
            context.set_source_rgba(motif_color[0], motif_color[1], motif_color[2])
            context.move_to(LEFT_MARGIN, y0)
            context.set_font_size(12)
            context.line_to(LEFT_MARGIN + 20, y0)
            context.stroke()
            context.move_to(LEFT_MARGIN + 30, y0)
            context.set_source_rgba(0,0,0)
            context.show_text(motif.upper())


class FastaHeader:
    '''Draws gene name, as given in Fasta header line'''

    def __init__(self, gene_name):
        self.gene_name = gene_name

    def draw(self, context, gene_number):
        y = HEIGHT_GENE_GROUP * gene_number + Y_OFFSET_GENE - 20
        context.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_source_rgba(0, 0, 0)
        context.set_font_size(14)
        context.move_to(LEFT_MARGIN, y)
        context.show_text(self.gene_name)


class GeneGroup:
    '''Tells children to draw themselves'''

    def __init__(self, exon, gene, motif, gene_name):
        self.exon = exon
        self.gene = gene
        self.motif = motif
        self.gene_name = gene_name
    
    def draw(self, context, gene_number):
        exon.draw(context, gene_number)
        gene.draw(context, gene_number)
        motif.draw(context, gene_number)
        gene_name.draw(context, gene_number)
        motif.draw_legend(context)
        



######MAIN#######

##parse motif file to produce dictionary of motifs
def DegNucleoConverter(motif):
    '''Take a motif as input and returns new motif containing regex characters for any degenerate nucleotides'''
    IUPAC_dict = {'U': '[TU]', 'T': '[TU]', 'W': '[ATU]', 'S': '[CG]', 'M': '[AC]', 'K': '[GTU]', 'R': '[AG]', 'Y': '[CTU]', 'B': '[CGTU]', 'D': '[AGTU]', 'H': '[ACTU]', 'V': '[ACG]', 'N': '[ACGTU]'}
    motif = motif.upper()
    motif_split = list(motif)
    for n, i in enumerate(motif_split):
        if i in IUPAC_dict:
            motif_split[n] = IUPAC_dict[i]
    motif = ''.join(motif_split)
    return motif

def MotifFileParser(motif_file):
    '''Takes a motif file as input and returns a dictionary of converted motifs with their original degenerate sequences as values'''
    with open (motif_file, "r") as mf:
        MOTIF_DICT = {}
        NUM_MOTIFS = 0
        for line in mf:
            NUM_MOTIFS += 1
            motif = line.strip()
            conv_motif = DegNucleoConverter(motif).lower()
            MOTIF_DICT.setdefault(conv_motif, motif)
    return MOTIF_DICT, NUM_MOTIFS

MOTIF_DICT, NUM_MOTIFS = MotifFileParser(motif_file)

#Define global variables
LEFT_MARGIN = 50
HEIGHT_GENE_GROUP = 80
Y_OFFSET_GENE = 150
MOTIF_COLOR_PAL = sns.color_palette("hls", NUM_MOTIFS)

#Parse FASTA file
def ParseFasta(input_file):
    '''Create svg output file from input fasta file'''
    records = {}
    read_counter = 0
    with open(input_file) as fh:
        for line in fh:
            line = line.strip()
            #check if line is a header line
            if line[0] == '>':
                gene_name = (line.split(" ")[0])[1:]
                records.setdefault(gene_name, "")
            else:
                records[gene_name] += str(line)
        return records

RECORDS = ParseFasta(input_file)

#MAIN
gene_number = 0
motif_number= 0
for i in RECORDS:
    gene_number += 1
    exon = Exon(RECORDS[i])
    gene = Gene(RECORDS[i])
    motif = Motif(RECORDS[i])
    gene_name = FastaHeader(i)
    gene_group = GeneGroup(exon, gene, motif, gene_name)
    gene_group.draw(context, gene_number)


surface.finish()


