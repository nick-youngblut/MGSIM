"""Class for simulating fragments from genomes"""
from __future__ import print_function
# import
## batteries
import os
import sys
import re
import time
from functools import partial
## 3rd party
from multiprocessing import Pool
import pandas as pd
from Bio import SeqIO
from Bio import Entrez


def main(args):
    """Main entry point
    """
    # loading genome accession file
    acc_tbl = parse_acc_table(args['<accession_table>'])    
    
    # creating put directory (if needed)
    if not os.path.isdir(args['-d']):
        os.makedirs(args['-d'])
    
    # query via efetch & download
    nprocs = int(float(args['-n']))
    pool = Pool(nprocs)
    func = partial(e_query, email=args['-e'], outdir=args['-d'])
    if args['--debug'] is True:
        acc_tbl = list(map(func, acc_tbl))
    else:
        acc_tbl = pool.map(func, acc_tbl)

    # writing out new acc_tbl to STDOUT
    print('\t'.join(['Taxon', 'Accession', 'Fasta']))
    for x in acc_tbl:
        # tidy taxon names
        x[0] = re.sub(r'[()\/:;, ]+', '_', str(x[0]))
        # writing line
        print('\t'.join(x))
        

def e_query(x, email, outdir):
    """Efetch query for one genome
    x : [taxon, accession]
    """
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
    out_file = os.path.join(outdir, genome_id + '.fna')
    outF = open(out_file, 'w')
    outF.write(''.join([x for x in records]) + '\n')
    outF.close()
    records.close()

    # status
    msg = 'File written: {}\n'
    sys.stderr.write(msg.format(out_file))
    return [genome_id, acc, out_file]
    
def parse_acc_table(infile):
    """Parsing tab-delim accession table (genome_name<tab>accession)
    """
    # load file
    if infile == '-':
        inF = sys.stdin
    else:
        inF = open(infile)
    df = pd.read_csv(inF, sep='\t')
    ## check headers
    diff = set(['Taxon','Accession']) - set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns:{}'.format(diff))

    # convert df to list of lists
    tbl = []
    for i,x in df.iterrows():
        tbl.append([x['Taxon'], x['Accession']])
    return tbl
        
