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
import urllib
from urllib.error import HTTPError
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
    func = partial(e_query, email=args['-e'], outdir=args['-d'],
                   rename=args['-r'], ambig_cutoff=int(args['-a']),
                   tries=int(args['-t']))
    if args['--debug'] is True or nprocs < 2:
        acc_tbl = list(map(func, acc_tbl))
    else:
        pool = Pool(nprocs)
        acc_tbl = pool.map(func, acc_tbl)

    # writing out new acc_tbl to STDOUT
    print('\t'.join(['Taxon', 'Accession', 'Fasta']))
    for x in acc_tbl:
        # writing line
        print('\t'.join([str(y) for y in x]))

def query_assembly(acc, genome_id, outdir, rename=False, ambig_cutoff=0, tries=10):
    """Querying assembly db
    """
    msg = '  Trying to query the "assembly" db with accession: {}'
    logging.info(msg.format(acc))
    
    # query assembly
    record = Entrez.read(Entrez.esearch(db="assembly", term=acc))
    accs = []
    for ID in record['IdList']:
        rec = Entrez.read(Entrez.esummary(db="assembly", id=ID, report="full"), validate=False)
        accession_id = rec['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
        for attempt in range(tries):
            try:
                refseq_id = Entrez.read(Entrez.esearch(db="nucleotide", term=accession_id), validate=False)['IdList'][0]
                break
            except HTTPError as e:
                logging.info('HTTPError! Error code: {}. Retrying...'.format(e.code))
                time.sleep((attempt + 1) ** 2 + 1)
        accs.append(refseq_id)
    if len(accs) < 1:
        return None
    if len(accs) > 1:
        msg = '  WARNING: >1 refseq accession for assembly: {}'
        logging.warning(msg.format(acc))

    # querying nucleotide
    acc = accs[0]
    records = None
    for attempt in range(tries):
        try:
            genomeIds = Entrez.read(Entrez.esearch(db="nucleotide",
                                                   term=acc,
                                                   retmode="xml",
                                                   retmax=1), validate=False)['IdList']
            records = Entrez.efetch(db="nucleotide",
                                    id=genomeIds,
                                    rettype="fasta",
                                    retmode="text")
        except IOError:
            msg = 'WARNING: Network error! Trying again for accession: {} (attempt {})'
            logging.warning(msg.format(acc, attempt+1))
            records = None
            time.sleep(10 + attempt * 2)
            continue
        if records is not None:
            records = list(records)
            if len(records) == 0 or records[0] == '' or \
               records[0] == 'Supplied id parameter is empty.':
                records = None
            else:
                break
    if records is None:
        sys.exit('Could not download accession: {}'.format(acc))
        
    # writing out file (renaming sequence headers if needed)
    if rename:
        out_file = os.path.join(outdir, '.' + str(uuid.uuid4()) + '.fna')
    else:
        out_file = os.path.join(outdir, genome_id + '.fna')
    with open(out_file, 'w') as outF:
        outF.write(''.join(records) + '\n')

    # checking genome fasta
    out_file = check_genome(out_file, ambig_cutoff)
    
    return out_file
        
def e_query(x, email, outdir, rename=False, ambig_cutoff=0, tries=10):
    """Efetch query for one genome
    x : [taxon, accession]
    """
    genome_id,acc = x[:2]
    logging.info('Downloading genome: {}'.format(genome_id))

    # output fasta file name
    if rename:
        out_file = os.path.join(outdir, '.' + str(uuid.uuid4()) + '.fna')
    else:
        out_file = os.path.join(outdir, genome_id + '.fna')
    
    # querying 
    Entrez.email = email

    ## query nucleotide db
    logging.info('  Querying "nucleotide" db with accession: {}'.format(acc))
    search = acc + '[Accession]'
    records = None
    for attempt in range(tries):
        try:
            genomeIds = Entrez.read(Entrez.esearch(db="nucleotide",
                                                   term=acc,
                                                   retmode="xml",
                                                   retmax=1))['IdList']
            records = Entrez.efetch(db="nucleotide",
                                    id=genomeIds,
                                    rettype="fasta",
                                    retmode="text")
        except IOError:
            msg = 'WARNING: Network error! Trying again for accession: {} (attempt {})'
            logging.warning(msg.format(acc, attempt+1))
            records = None
            time.sleep(10 + attempt * 2)
            continue
        if records is not None:
            records = list(records)
            if len(records) == 0 or records[0] == '' or \
               records[0] == 'Supplied id parameter is empty.':
                records = None
            else:
                break
    if records is None:
        sys.exit('Could not download accession: {}'.format(acc))
        
    # writing out file (renaming sequence headers if needed)
    with open(out_file, 'w') as outF:
        outF.write(''.join(records) + '\n')

    # checking that genome is formatted correctly formatted (and complete)
    out_file = check_genome(out_file, ambig_cutoff)
    if out_file is None:
        out_file = query_assembly(acc, genome_id, outdir, rename=rename,
                                  ambig_cutoff=ambig_cutoff)

    # renaming sequence headers (if needed)    
    if out_file is not None and rename:
        in_file = out_file
        out_file = os.path.join(outdir, genome_id + '.fna')
        genome_rename(in_file, out_file, genome_id)
        
    # status
    msg = '  File written: {}'
    logging.info(msg.format(out_file))
    return [genome_id, acc, out_file]

def check_genome(genome_fasta, ambig_cutoff=0):
    """Checking that genome fasta is formatted correctly
    """
    ambig_chars = re.compile(r'[RYSWKMBVDH]')
    N_cnt = 0
    ambig_cnt = 0
    len_cnt = 0
    for seq_record in SeqIO.parse(genome_fasta, 'fasta'):
        len_cnt += len(seq_record)
        seq = seq_record.seq.upper()
        N_cnt += seq.count('N')
        ambig_cnt += len(ambig_chars.findall(str(seq)))
            
    if ambig_cnt > 0:
        if ambig_cnt > ambig_cutoff:
            msg = '  WARNING: no. of ambig nucleotides is {}; removing genome: {}'
            logging.warning(msg.format(ambig_cnt, genome_fasta))
            return(None)
        else:
            msg = '  WARNING: no. of ambig nucleotides is {}'
            msg += ', but ambig-cutoff is {}, so keeping genome: {}'
            logging.warning(msg.format(ambig_cnt, ambig_cutoff, genome_fasta))            

    if len_cnt < 1000:
        msg = '  WARNING: Genome length < 1000bp; removing genome: {}'
        logging.warning(msg.format(genome_fasta))
        os.remove(genome_fasta)
        return(None)
    if N_cnt / float(len_cnt) > 0.5:
        msg = '  WARNING: Genome sequences are >50% N\'s; removing genome: {}'
        logging.warning(msg.format(genome_fasta))
        os.remove(genome_fasta)
        return(None)

    return genome_fasta

def genome_rename(in_file, out_file, taxon):
    """Renaming sequences in genome file.
    All will be named based on `taxon`
    """
    logging.info('  Renaming genome seqs to: {}'.format(taxon))

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
    # checking for duplicates
    if any(df.duplicated(subset=['Taxon'])):
        raise ValueError('Duplicate taxon names in the accession table')
    
    # convert df to list of lists
    tbl = []
    for i,x in df.iterrows():
        tbl.append([x['Taxon'], x['Accession']])
    return tbl
        
