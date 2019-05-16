#!/usr/bin/env python3

import sys,gzip,argparse
from Bio import SeqIO, Seq

parser = argparse.ArgumentParser(description='Preprocess fasta_file for Ribo-Seq data.')

parser.add_argument('-u','--up_file', help='<Required> Upstream sequence file.', required=True)
parser.add_argument('-d','--down_file', help='<Required> Downstream sequence file.', required=True)
parser.add_argument('-o','--out_fasta', help='<Required> Output fasta sequence file.', required=True)
parser.add_argument('-i','--in_fasta', help='<Required> Input fasta sequence file.', required=True)
parser.add_argument('-k', help='<Required> Number of nts to complete the utr5 and utr3 to.', required=True)


label_args = parser.add_mutually_exclusive_group()
label_args.add_argument('-ex','--exclude_labels', nargs='+', help='<Required> Exclude all transcripts with these labels.')
label_args.add_argument('-in','--include_labels', nargs='+', help='<Required> Include transcripts only with these labels.')


# File opener function
def f_open(filename, mode='r'):
    if filename.endswith(".gz"):
        if mode=="r":
            return(gzip.open(filename,"rt"))
        else:
            return(gzip.open(filename, mode))
    else:
        return(open(filename,mode))


args = parser.parse_args()

# Parse upstream biomart fasta file
TX_up = {}
with f_open(args.up_file) as handle:
    for record in SeqIO.parse(handle, "fasta") :
        # Ensembl transcript IDs are 18characters
        TX_up[str(record.id)[:18]] = str(record.seq)

TX_down = {}
with f_open(args.down_file) as handle:
    for record in SeqIO.parse(handle, "fasta") :
        # Ensembl transcript IDs are 18characters
        TX_down[str(record.id)[:18]] = str(record.seq)


# Parse input fasta that also includes CDS info / write into the output file
with f_open(args.out_fasta,'w') as out_f:
    with f_open(args.in_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta") :
            id_cols = record.id.split('|')
            cds = None
            utr5 = None
            utr3 = None

            seq = str(record.seq)
            for i in range(len(id_cols)):
                if id_cols[i].startswith('CDS'):
                    cds = (i,(int(id_cols[i].split(':')[1].split('-')[0]),int(id_cols[i].split(':')[1].split('-')[1])))
                elif id_cols[i].startswith('UTR5'):
                    utr5 = (i,(int(id_cols[i].split(':')[1].split('-')[0]),int(id_cols[i].split(':')[1].split('-')[1])))
                elif id_cols[i].startswith('UTR3'):
                    utr3 = (i,(int(id_cols[i].split(':')[1].split('-')[0]),int(id_cols[i].split(':')[1].split('-')[1])))


            if len(id_cols) > 7 and ((args.include_labels==None) or (id_cols[7] in args.include_labels)) and ((args.exclude_labels==None) or (id_cols[7] not in args.exclude_labels)):
                up = False

                if cds[1][0]<(int(args.k)+1):
                    upstream_seq = TX_up[str(record.id)[:18]]
                    if cds!=None:
                        id_cols[cds[0]] = 'CDS:'+str((int(args.k)+1))+'-'+str(cds[1][1]+((int(args.k)+1)-cds[1][0]))
                    if utr3!=None:
                        id_cols[utr3[0]] = 'UTR3:'+str(utr3[1][0]+((int(args.k)+1)-cds[1][0]))+'-'+str(utr3[1][1]+((int(args.k)+1)-cds[1][0]))
                    if utr5!=None:
                        id_cols[utr5[0]] = 'UTR5:1-'+str(int(args.k)+1)
                    else:
                        id_cols = id_cols[:cds[0]]+['UTR5:1-'+args.k]+id_cols[cds[0]:]
                        cds = (cds[0]+1,cds[1])

                    if id_cols[-1] == '':
                        id_cols[-1] = 'up'+str((int(args.k)+1)-cds[1][0])
                    else:
                        id_cols.append = 'up'+str((int(args.k)+1)-cds[1][0])

                    seq = upstream_seq[((int(args.k)+1)-cds[1][0]):] + seq
                    up = True

                if (len(record.seq)-cds[1][1]) < int(args.k):
                    downstream_seq = TX_down[str(record.id)[:18]]

                    if utr3!=None:
                        if cds[1][1]+1 != utr3[1][0]:
                            print('PROBLEM !!!')
                            print(record.id)
                        id_cols[utr3[0]] = 'UTR3:'+str(utr3[1][0])+'-'+str(utr3[1][0]+(int(args.k)-1))
                    else:
                        id_cols = id_cols[:(cds[0]+1)]+['UTR3:'+str(cds[1][1]+1)+'-'+str(cds[1][1]+int(args.k))]+ id_cols[(cds[0]+1):]

                    seq = seq + downstream_seq[:(int(args.k)-(len(record.seq)-cds[1][1]))]
                    if up:
                        id_cols[-1] += 'down'+str(int(args.k)-(len(record.seq)-cds[1][1]))
                    elif id_cols[-1] == '':
                        id_cols[-1] = 'down'+str(int(args.k)-(len(record.seq)-cds[1][1]))
                    else:
                        id_cols.append = 'down'+str(int(args.k)-(len(record.seq)-cds[1][1]))

                out_f.write('>'+str('|'.join(id_cols))+'|\n'+seq+'\n')

            else:
                print('PROBLEM !!!')
                print(record.id)
                exit()
