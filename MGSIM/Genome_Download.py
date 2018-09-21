"""Class for simulating fragments from genomes"""

# import
## batteries
import sys, os
import time
from functools import partial
## 3rd party
from multiprocessing import Pool
import pandas as pd
from Bio import SeqIO
from Bio import Entrez


def main(args):
    """Main entry point"""

    # loading genome accession file
    acc_tbl = parse_acc_table(args['<accession_table>'])    

    # creating put directory (if needed)
    if not os.path.isdir(args['-d']):
        os.makedirs(args['-d'])
    
    # query via efetch & download
    nprocs = int(args['-n'])
    pool = Pool(nprocs)
    func = partial(e_query, email=args['-e'], outdir=args['-d'])
    if args['--debug']:
        map(func, acc_tbl)
    else:
        pool.map(func, acc_tbl)
            

def e_query(x, email, outdir):
    genome_id,acc = x[:2]

    # querying 
    Entrez.email = email
    search = acc + '[Accession]'

    while 1:
        try:
            handle = Entrez.read(Entrez.esearch(db="nucleotide",
                                                term=search,
                                                retmode="xml",
                                                retmax=1))
            genomeIds = handle['IdList']
            records = Entrez.efetch(db="nucleotide",
                                    id=genomeIds,
                                    rettype="fasta",
                                    retmode="text")
        except IOError:
            msg = 'Network error! Trying again for accession: {}\n'
            sys.stderr.write(msg.format(acc))
            time.sleep(10)
            continue
        break
        
    # writing out file
    outfile = os.path.join(outdir, genome_id + '.fna')
    outF = open(outfile, 'w')
    outF.write(''.join([x for x in records]) + '\n')
    outF.close()
    records.close()

    # status
    msg = 'File written: {}\n'
    sys.stderr.write(msg.format(outfile))

    
def parse_acc_table(infile):
    """Parsing tab-delim accession table (genome_name<tab>accession)
    """
    if infile == '-':
        inF = sys.stdin
    else:
        inF = open(infile)

    tbl = []
    for line in inF:
        line = line.rstrip().split('\t')
        tbl.append(line)

    return tbl
        
