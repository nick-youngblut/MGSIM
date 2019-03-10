"""Class for simulating fragments from genomes"""
from __future__ import print_function
# import
## batteries
import os
import sys
import re
import time
import uuid
import logging
import tempfile
from functools import partial
## 3rd party
from multiprocessing import Pool
import pandas as pd
from Bio import SeqIO
from Bio import Entrez

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


# functions
def main(args):
    """Main entry point
    """
    # loading genome accession file
    acc_tbl = parse_acc_table(args['<accession_table>'])    
    
    # creating put directory (if needed)
    if not os.path.isdir(args['-d']):
        os.makedirs(args['-d'])
    
    # query via efetch & download
    logging.info('Downloading genomes via efetch...')
    nprocs = int(float(args['-n']))
    pool = Pool(nprocs)
    func = partial(e_query, email=args['-e'], outdir=args['-d'], rename=args['-r'])
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
        
def e_query(x, email, outdir, rename=False):
    """Efetch query for one genome
    x : [taxon, accession]
    """
    genome_id,acc = x[:2]
    logging.info('Downloading genome: {}'.format(genome_id))

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
    if rename:
        out_file = os.path.join(outdir, '.' + str(uuid.uuid4()) + '.fna')
    else:
        out_file = os.path.join(outdir, genome_id + '.fna')
    with open(out_file, 'w') as outF:
        outF.write(''.join([x for x in records]) + '\n')
    records.close()

    if rename:
        in_file = out_file
        out_file = os.path.join(outdir, genome_id + '.fna')
        genome_rename(in_file, out_file, genome_id)
    
    # status
    msg = 'File written: {}\n'
    logging.info(msg.format(out_file))
    return [genome_id, acc, out_file]

def genome_rename(in_file, out_file, taxon):
    """Renaming sequences in genome file.
    All will be named based on `taxon`
    """
    logging.info('Renaming genome seqs to: {}'.format(taxon))

    # number of sequences in the genome fasta file
    seq_cnt = 0
    with open(in_file) as inF:
        for line in inF:
            if line.startswith('>'):
                seq_cnt += 1
                
    # renaming sequence headers
    with open(in_file) as inF, open(out_file, 'w') as outF:
        idx = 1
        for line in inF:
            if line.startswith('>'):
                if seq_cnt < 2:
                    outF.write('>{}\n'.format(taxon))
                else:
                    outF.write('>{}_{}\n'.format(taxon, idx))
                    idx += 1
            else:
                outF.write(line)

    # removing temp file
    os.remove(in_file)

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
    # removing non-safe characters from taxon labels
    df['Taxon'] = df['Taxon'].apply(lambda x: re.sub('[^A-Za-z0-9-]+', '_', x))
    df['Taxon'] = df['Taxon'].apply(lambda x: re.sub('_+$', '', x))
    # convert df to list of lists
    tbl = []
    for i,x in df.iterrows():
        tbl.append([x['Taxon'], x['Accession']])
    return tbl
        
