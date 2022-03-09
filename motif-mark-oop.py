import argparse
import re
import cairo
import Bio_Module
from operator import itemgetter
from itertools import *

# Creation of Degerante Base Reference 
Degenerate_Bases = {'A':'[A]', 'C':'[C]', 'T':'[T]', 'G':'[G]', 'U':'[U]',
'W':'[A|T]', 'S': '[C|G]',
'M':'[A|C]', 'K':'[G|T]',
'R':'[A|G]', 'Y':'[C|T]',
'B':'[C|G|T]', 'D':'[A|G|T]', 'H':'[A|C|T]', 'V':'[A|C|G]', 'N':'[A|C|G|T]'}

# Get command Line Arguments
def get_args():
    parser = argparse.ArgumentParser(
        description = ""
    )
    parser.add_argument(
        "-f",
        help="Input fasta file", 
        required=True,
        ),
    parser.add_argument(
        "-m",
        help="Input motifs file",
        required=True,
    ),
    parser.add_argument(
        "--Leslie",
        help="Enables squatting commands",
        action ='store_true',
        default=False
    )
    return parser.parse_args()

# Gene Class
class Gene:

    def __init__(self, sequence, name, file_name, motif_targets = {}):
        self.sequence = sequence
        self.name = name
        self.motif_targets = motif_targets
        self.exons = ''
        self.introns = ''
        self.file_name = file_name


    def find_exons_and_introns(self):
        self.exons = [index for index in range(len(self.sequence)) if self.sequence[index].isupper()]
        self.introns = [index for index in range(len(self.sequence)) if self.sequence[index].islower()]
        self.exons = [list(map(itemgetter(1), g)) for k, g in groupby(enumerate(self.exons), lambda x: x[0]-x[1])]
        

    def find_motifs(self,motifs):
        
        for motif in motifs:
            matches = ([(m.start(0), m.end(0)) for m in re.finditer(motif.ambiguous, self.sequence.upper())])
            self.motif_targets[motif.sequence] = matches
        

# Motif Class
class Motif:

    def __init__(self,sequence):
        self.sequence = sequence
        self.ambiguous = ''
        self.lower = any(c.islower() for c in sequence)
    
    def find_ambiguous_sequences(self,sequence,ambiguous_dict):
        # Construct ambiguous regex expression
        self.ambiguous = ''
        for nucleotide in self.sequence:
            self.ambiguous += Degenerate_Bases[nucleotide.upper()]

# Function that creates all Motif objects and initalizes values     
def create_motifs(motif_file):
    motifs = []
    with open (motif_file, 'r') as f:
        for line in f:
            x = line.strip()
            x = Motif(x)
            motifs.append(x)  
    for motif in motifs:
        motif.find_ambiguous_sequences(motif.sequence,motif.ambiguous)
    return motifs

# Function that creates all Gene objects and initalizes values
def create_genes(gene_file, motifs):
    genes = []
    one_line_fasta = Bio_Module.oneline_fasta(gene_file)
    for line in one_line_fasta:
        gene, sequence = line
        name = gene.split(' ')[0][1:]
        motif = {}
        gene = Gene(sequence,name,gene_file,motif)
        genes.append(gene)
    for gene in genes:
        gene.find_exons_and_introns()
        gene.find_motifs(motifs)
    return genes

# Function to Draw PNG
def create_marked_motifs(genes,output_name):
    HEIGHT = 200 + len(genes) * 200
    WIDTH = 2000
    ims = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    cr = cairo.Context(ims)

    # Draw Title
    cr.set_source_rgb(1, 1, 1)
    cr.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, 
        cairo.FONT_WEIGHT_NORMAL)
    cr.set_font_size(80)
    (x, y, width, height, dx, dy) = cr.text_extents("Motif Mappings")
    cr.move_to(WIDTH/2 - width/2,100)
    cr.show_text("Motif Mapping")

    # Draw Motif thing

    cr.rectangle(1400, HEIGHT/2, 8*len(genes[0].motif_targets), 80)
    cr.fill()

    # Draw Genes
    cr.set_font_size(40)
    for count,gene in enumerate(genes):
        cr.set_source_rgb(0.9, 0.9, 0.9)
        target_height = (count + 1) * 200
        cr.move_to(100, target_height)
        cr.show_text(gene.name)
        cr.move_to(100, target_height + 50)
        cr.set_line_width(5)
        cr.line_to(len(gene.sequence),target_height + 50)
        cr.stroke()

        # Loop through exon list
        for i in gene.exons:
            
            cr.rectangle(i[0], target_height + 35, len(i), 30)
            cr.fill()

        # Add Motifs 
        for item in gene.motif_targets.items():
            name,sequence = item
            if sequence:
                cr.set_source_rgb(0.3, 0.3, 1)
                for hit in sequence:
                    print(hit[1] - hit[0])
                    cr.rectangle(hit[0] + 100, target_height + 40, hit[1] - hit[0], 20)
                    cr.fill()

    output = output_name.split('.')[0]
    ims.write_to_png(f'{output}.png')
    
def main():
    args = get_args()
    motifs = create_motifs(args.m)
    genes = create_genes(args.f,motifs)
    create_marked_motifs(genes,args.f)
    print(genes[0].motif_targets)
    print(len(genes[0].sequence))
   

if __name__ == "__main__":    
    main()
    