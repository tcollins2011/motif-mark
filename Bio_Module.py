#populate this file with any fuctions you have written so far
#We suggest convert_phred() in particular :)
#Import this into your PS4.ipynb by issuing the command: import Bioinfo


# Constants 
DNA_bases = ['G','C','T','A','N']
RNA_bases = ['G','C','U','A','N']


# Functions 
def convert_phred(letter):
    """Converts a single character into a phred score"""
    return((ord(letter))-33)

def qual_score(phred_score):
    """Determins the mean quality score"""
    total = 0
    for i in phred_score:
        total += convert_phred(i)
    total = total / len(phred_score)
    return total

def validate_base_seq(seq, RNAflag):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts, Gs, and Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()         
    Gs = DNA.count("G")       
    Cs = DNA.count("C")       
    return (Gs+Cs)/len(DNA)

def oneline_fasta(fasta):
    """Takes in a multline fasta file and coverts it to a single line fasta file"""
    header = ''
    sequence = ''

    if fasta == '':
        return 
    
    # Open and read thorugh file
    target_file = open(fasta,'r') 
    line = target_file.readline()
    while not line.startswith('>'):
        line = target_file.readline()
    header = line.strip()

    for line in target_file:
        if line.startswith ('>'):
            # Returns header and sequence resets line and adds new header
            yield header,sequence
            header = line.strip()
            sequence = ''

            # Append sequence to sequence variable 
        else:
            sequence += ''.join(line.strip().split())

    # Final yield
    yield header,sequence
    target_file.close()


if __name__=="__main__":
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")


    assert qual_score(["I","C"]) == 37, "wrong average phred score"
    print("You calcluated the correct average phred score")

    assert validate_base_seq("TATUC",False) == False
    assert validate_base_seq("UCUGCU", False) == False
    print("Passed DNA and RNA tests")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")

    test_fasta_a = dict()
    test_fasta_a.update(oneline_fasta('fasta_test.fa'))
    assert len(test_fasta_a) == 1
    print('combined fasta a to one line in dict')