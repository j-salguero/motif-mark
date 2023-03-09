#! /usr/bin/env python

import argparse
import re
import math
import cairo

def get_args():
    parser = argparse.ArgumentParser(description="Motif Mark code will parse FASTA file and motif file inputs. It will then locate the positions of these motifs found in each gene sequence, outputting one image per FASTA file.", add_help=False) #PUT HELP MESSAGE HERE
    parser.add_argument("-f", "--fasta", help="Absolute path for fasta file", required=True, type=str)
    parser.add_argument("-m", "--motif", help="Absolute Path for motif file", required=True, type=str)
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Motif mark code assumes that the motif file and FASTA files are correctly formatted.')
    return parser.parse_args()
args = get_args()

#if the sequences in the FASTA file specify that they are reverse complement, then this function will have return the correct sequence
def rev_comp(the_seq):
    '''Input a DNA or RNA sequence, return the reverse complement of the sequence.'''
    comp_seq = ""
    if(("U" in the_seq) == True): #RNA Sequences
        for i in the_seq:
            if i == "A":
                comp_seq = comp_seq + "U"
            elif i == "G":
                comp_seq = comp_seq + "C"
            elif i == "U":
                comp_seq = comp_seq + "A"
            elif i == "C":
                comp_seq = comp_seq + "G"
            elif i == "a":
                comp_seq = comp_seq + "u"
            elif i == "g":
                comp_seq = comp_seq + "c"
            elif i == "u":
                comp_seq = comp_seq + "a"
            elif i == "c":
                comp_seq = comp_seq + "g"
    elif(("T" in the_seq) == True): #DNA Sequences
        for i in the_seq:
            if i == "A":
                comp_seq = comp_seq + "T"
            elif i == "G":
                comp_seq = comp_seq + "C"
            elif i == "T":
                comp_seq = comp_seq + "A"
            elif i == "C":
                comp_seq = comp_seq + "G"
            elif i == "a":
                comp_seq = comp_seq + "t"
            elif i == "g":
                comp_seq = comp_seq + "c"
            elif i == "t":
                comp_seq = comp_seq + "a"
            elif i == "c":
                comp_seq = comp_seq + "g"
    #switch the orientation of the complementary sequence to get the final reverse complement
    rev_seq = comp_seq[::-1]
    return(rev_seq)

class FASTA_Genes:
    def __init__(self, the_file):
        '''This is how a FASTA file is read and stored.'''
        ## Data ##
        self.file_name = the_file
        self.gene_seqs = {} #key = gene name, value = seq
        self.introns = {}
        self.exons = {}
        self.headers = {}

    def oneline_fasta(self, output_file):
        '''Take in string parameters of the input file name where sequences have '\n' in the middle, and the output file name. 
            Function removes '\n' from the middle of sequences in the input file. Output file will have headers on one line, 
            entire sequence on one line. Will return the string 'Successful' when the function is complete'''
        with open(self.file_name) as input:
            with open(output_file, 'w') as output:
                counter = 0 
                for line in input:
                    if ">" in line and counter!=0:
                        output.write('\n')
                    if ">" not in line:
                        line=line.strip('\n')
                    output.write(line)
                    counter+=1
            self.file_name = output_file    
 
    def separate_genes(self):
        '''Separate individual genes and their sequences into a dictionary'''
        with open(self.file_name) as input:
            gene_key = ""
            rev_check = 0
            for line in input:
                if ">" in line: #header lines
                    #get the gene name
                    gene_key = re.search(r'[^(\s)]+', line)
                    gene = gene_key.group()
                    #remove the '>' from the front of the gene name
                    gene_key = gene[1:len(gene)]
                    self.headers[gene_key] = line.strip()
                    #check if the header specifies that the sequence is the reverse complement
                    rev_check = line.find("reverse complement")
                else: #sequence lines
                    if(rev_check != -1):
                        #replace the original sequence with its reverse complement
                        line = rev_comp(line)
                    if(gene_key not in self.gene_seqs.keys()):
                        self.gene_seqs[gene_key] = line
                    #determine the start and stop positions of introns and exons
                    self.intron_exon(line, gene_key)
                    gene_key = "" #reset the gene_key for the next iteration

    def intron_exon(self, the_sequence, the_gene):
        '''Read in sequence and determine the start and end positions of introns and exons. Introns are lowercase letters, exons are uppercase letters.'''
        #search for the start and end positions of the introns
        #introns are in lowercase: [a-z]*
        introns_pos = []
        for matches_I in re.finditer(r'[a-z]*', the_sequence):
            if(matches_I.end() - matches_I.start() != 0):
                intron_add = [matches_I.start(), matches_I.end()]
                introns_pos.append(intron_add)
        self.introns[the_gene] = introns_pos

        #search for the start and end positions of the exons
        #exons are in uppercase: '[A-Z]*'
        exons_pos = []
        for matches_E in re.finditer(r'[A-Z]*', the_sequence):
            if(matches_E.end() - matches_E.start() != 0):
                exon_add = [matches_E.start(), matches_E.end()]
                exons_pos.append(exon_add)
        self.exons[the_gene] = exons_pos

       
class Motifs:
    def __init__(self, the_motif_file, the_fasta):
        '''This is how a Motif file is read and stored.'''
        ## Data ##
        self.seqs = {}
        self.fasta_file = the_fasta.gene_seqs
        self.positions = {} #key = gene name, value = list of start/end positions
        self.motif_file = the_motif_file
        self.degenerate_seqs = {}
        self.read_motifs() #read in the motifs on initialization

    def read_motifs(self):
        '''Read in the motifs file. Each line of the txt file is an individual motif.'''
        with open(self.motif_file) as input: 
            for line in input:
                line=line.strip('\n')
                self.seqs[line] = ""

    def degenerate_ambiguity(self):
        '''Replace degernate bases with ALL posibilities'''
        #Y = [CGUcgu...]
        #replace regular bases with upper and lower case [Aa]
        #If it matches to a T, it shoudl also match to U to make it applicable to both RNA/DNA sequences
        ambig_bases = {'A':'[Aa]', 'C':'[Cc]', 'G':'[Gg]', 'T':'[TtUu]', 'U':'[UuTt]', 'W':'[ATUatu]', 'S':'[CGcg]', 'M':'[ACac]', 'K':'[GTUgtu]', 'R':'[AGag]', 'Y':'[CTUctu]', 'B':'[CGTUcgtu]', 'D':'[AGTUagtu]', 'H':'[ACTUactu]', 'V':'[ACGacg]', 'N':'[ACGTUacgtu]'}
        for key,value in self.seqs.items():
            key_upper = key.upper()
            value_str = ""
            for i in key_upper:
                value_str += ambig_bases[i]
            self.seqs[key] = value_str
            self.search_for(key, value_str)

    def search_for(self, the_ambig, the_possible):
        '''Search FASTA file for all locations where the ambiguous motifs are found. This will create a dictionary with the ambiguous sequence and all start/end positions'''
        #regex equals all of the ambigous possibilities for one specific motif
        regex = the_possible
        for gene, seq in self.fasta_file.items():
           for matches in re.finditer(rf"{regex}", seq):
                if((the_ambig, gene) not in self.positions.keys()):
                    self.positions[the_ambig, gene] = []
                #each individual occurrance of the motif will be stores into a list of lists
                #each occurrance will be stored in a list [start, end position]
                added_list = [matches.start(), matches.end()]
                self.positions[the_ambig, gene].append(added_list)

    def draw_pic(self, the_fasta):
        '''Create a pycairo image depicting introns, exons, and all motifs. This will produce one image per FASTA file'''
        #create the surface and context
        width, height = 1100, (len(the_fasta.gene_seqs.items())*75+75)
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        context = cairo.Context(surface)
        context.save()
        context.set_source_rgb(1, 1, 1)
        context.paint()
        counter = 0.5
        for gene, pos_groups in the_fasta.introns.items():
            #draw a line for introns
            height_line = counter * 75
            context.set_source_rgb(0, 0, 0)
            context.set_line_width(4)
            context.move_to(22, height_line - 20)
            context.show_text(the_fasta.headers[gene])
            context.move_to(22, height_line)        #(x,y)
            context.line_to(int(pos_groups[len(pos_groups)-1][1]) + 22,height_line)
            context.stroke()

            #draw a line for exons
            exon = the_fasta.exons[gene]
            for i in exon:
                context.set_line_width(20)
                context.set_source_rgb(0, 0, 0)
                context.move_to(int(i[0]) + 22, height_line)
                context.line_to(int(i[1]) + 22, height_line)
                context.stroke()

            #draw in the motifs
            #purple, green, red, orange, light blue
            colors = [[0.7, 0.5, 1 ], [0.7,0.9,0.2], [0.8, 0.2, 0.1], [1, 0.6, 0], [0.7, 0.9, 1]]
            color_dict = {}
            m_counter = 0 
            for i in self.seqs.keys():
                color_dict[i] = colors[m_counter]
                m_counter+=1
            for m_seq, m_gene in self.positions.keys():
                if(m_gene == gene):
                    for pos in self.positions[(m_seq, gene)]:
                        motif_color = color_dict[m_seq]
                        context.set_line_width(20)
                        context.set_source_rgb(motif_color[0], motif_color[1], motif_color[2])
                        context.move_to(22 + pos[0], height_line)        #(x,y)
                        context.line_to(pos[1] + 22, height_line)                        
                        context.stroke_preserve()
                        context.fill()
                        context.save()
                        context.stroke()
                        context.restore()
            counter += 1
        #write the legend
        #repeat the show_text commands twice to have a bolder text size
        color_counter = 1
        for gene_name, motif_c in color_dict.items():
            context.set_line_width(15)
            context.move_to(26*color_counter, height - 25)
            context.line_to(26*color_counter + 7, height - 25)
            context.set_source_rgb(motif_c[0], motif_c[1], motif_c[2])
            context.stroke()
            context.move_to(26*color_counter + 20, height - 25)
            context.set_source_rgb(0,0,0)
            context.show_text(gene_name)
            context.move_to(26*color_counter + 20, height - 25)
            context.show_text(gene_name)
            color_counter += 3
        context.move_to(26*(color_counter/2), height - 50)
        context.show_text("Legend:")
        context.move_to(26*(color_counter/2), height - 50)
        context.show_text("Legend:")
        #output the final image to a PNG file
        surface.write_to_png(file_prefix[len(file_prefix) - 2] + ".png")
        surface.finish()


#Where the functioning code begins:
#End goal: single, well-labeled figure per FASTA file
FASTA1 = FASTA_Genes(args.fasta)
file_prefix = args.fasta.split('.')
FASTA1.oneline_fasta(file_prefix[len(file_prefix) - 2] + "_oneline.fasta")
FASTA1.separate_genes()
Motif1 = Motifs(args.motif, FASTA1)
Motif1.degenerate_ambiguity()
Motif1.draw_pic(FASTA1)