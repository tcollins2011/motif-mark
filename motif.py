import cairo

#As of the moment I’m thinking of two major classes for this assignment.  
# The Gene class will contain most of the functionality and will have as many instances as there are genes in the input fasta file. 
# The motif class will have an instance for each motif passed into the program through a text file. 

#GENE:
#•	Gene will be initialized with an input sequence coming from the fasta file that we will be working with. It will have a sequence variable, found_motif_dictionary, and some undecided data structure containing exon and intron location data.
#•	It will contain a function that takes in a motif and records all of the matches that the motif has within its sequence. This will be recorded as a start and end position. These will be saved within a dictionary inside of GENE.
#•	It will also contain a function that returns the locations of the introns and exons of the gene
#•	It will contain a draw function that will use pycairo to draw a png with the exons, introns, and motif hits marked.

#Motif:
#•	Motif will initialize with a motif sequence and save it as a variable within the class
#•	Motifs will also be designated with some RGB color and save that variable
#•	It will contain a function that determines which nucleotide sequences will work with that motif. Ie) if both ACGA and ACGT could possibly be valid for some motif

class Gene:
    pass

class Motif:
    pass


HEIGHT = 400
WIDTH = 400
       
def main():
    

    ims = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    cr = cairo.Context(ims)


    cr.set_source_rgb(0.6, 0.6, 0.6)

    cr.rectangle(20, 20, 120, 80)
    cr.fill()

    cr.rectangle(180, 20, 80, 80)
    cr.fill()
    
    cr.set_line_width(1.5)
    cr.move_to(120,150)
    cr.line_to(360, 350)
    cr.stroke()

    ims.write_to_png("image.png")
        
        
if __name__ == "__main__":    
    main()