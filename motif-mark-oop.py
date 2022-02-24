import argparse
import re
import cairo
import Bio_Module

# Creation of Degerante Base Reference 
Degenerate_Bases = {'A':'[A]', 'C':'[C]', 'T':'[T]', 'G':'[G]', 'U':'[U]',
'W':'[A|T]', 'S': '[C|G]',
'M':'[A|C]', 'K':'[G|T]',
'R':('[A|G]'), 'Y':'[C|T]',
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

    def __init__(self, sequence, file_name, motif_targets = {}):
        self.sequence = sequence
        self.motif_targets = motif_targets
        self.exons = ''
        self.introns = ''
        self.file_name = file_name


    def find_exons_and_introns(self):
        self.exons = [index for index in range(len(self.sequence)) if self.sequence[index].isupper()]
        self.introns = [index for index in range(len(self.sequence)) if self.sequence[index].islower()]

    def find_motifs(self,motifs):
        
        for motif in motifs:
            matches = ([(m.start(0), m.end(0)) for m in re.finditer(motif.ambiguous, self.sequence.upper())])
            self.motif_targets[motif.sequence] = matches
        
    def draw(self):
        # one graph per fasta file
        # Should have the same name as the input fasta file 
        pass

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

# Function taht creates all Motif objects and initalizes values     
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
        gene = Gene(sequence,gene_file)
        genes.append(gene)
    for gene in genes:
        gene.find_exons_and_introns()
        gene.find_motifs(motifs)
    return genes

# Function to Draw PNG
def create_marked_motifs(genes,output_name):
    HEIGHT = 2000
    WIDTH = 2000
    ims = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    cr = cairo.Context(ims)

    for gene in genes:
        pass
    # cr.set_source_rgb(0.6, 0.6, 0.6)

    # cr.rectangle(20, 20, 120, 80)
    # cr.fill()

    # cr.rectangle(180, 20, 80, 80)
    # cr.fill()
    
    # cr.set_line_width(1.5)
    # cr.move_to(120,150)
    # cr.line_to(360, 350)
    # cr.stroke()

    output = output_name.split('.')[0]
    ims.write_to_png(f'{output}.png')
    
def main():
    args = get_args()
    motifs =create_motifs(args.m)
    genes = create_genes(args.f,motifs)
    create_marked_motifs(genes,args.f)
   

if __name__ == "__main__":    
    main()
    