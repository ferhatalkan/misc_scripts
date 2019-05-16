import sys,gzip,argparse
from Bio import SeqIO, Seq

parser = argparse.ArgumentParser(description='Preprocess fasta_file for Ribo-Seq data.')

parser.add_argument('-u','--up_file', nargs='+', help='<Required> Upstream sequence file.', required=True)
parser.add_argument('-o','--out_fasta', nargs='+', help='<Required> Output fasta sequence file.', required=True)
parser.add_argument('-i','--in_fasta', nargs='+', help='<Required> Input fasta sequence file.', required=True)

label_args = parser.add_mutually_exclusive_group()
label_args.add_argument('-ex','--exclude_labels', nargs='+', help='<Required> Exclude all transcripts with these labels.', required=True)
label_args.add_argument('-in','--include_labels', nargs='+', help='<Required> Include transcripts only with these labels.', required=True)


# File opener function
def f_open(filename, mode='r'):
    if filename.endswith(".gz"):
        if mode=="r":
            return(gzip.open(filename,"rt"))
        else:
            return(gzip.open(filename, mode))
    else:
        return(open(filename,mode))

# Parse upstream biomart fasta file
TX_up = {}
with f_open(args.up_file) as handle:
    for record in SeqIO.parse(handle, "fasta") :
        # Ensembl transcript IDs are 18characters
        TX_up[str(record.id)[:18]] = str(record.seq)

# Parse input fasta that also includes CDS info / write into the output file
with f_open(args.out_fasta,'w') as out_f:
    with f_open(args.in_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta") :
            id_cols = record.id.split('|')
            for col in id_cols:
                if col.startswith('CDS:') and len(id_cols) >= 8 and ((len(args.include_labels)==0 or (id_cols[7] in args.include_labels)) and (id_cols[7] not in args.exclude_labels)):
                    upstream_seq = TX_up[str(record.id)[:18]]
                    out_f.write('>'+str(record.id)+)
                    break
        

