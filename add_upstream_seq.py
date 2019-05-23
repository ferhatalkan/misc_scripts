#!/usr/bin/env python3

import sys,gzip,argparse
from Bio import SeqIO, Seq

parser = argparse.ArgumentParser(description='Preprocess fasta_file for Ribo-Seq data.')

parser.add_argument('-u','--up_file', help='<Required> Upstream sequence file.', required=True)
parser.add_argument('-d','--down_file', help='<Required> Downstream sequence file.', required=True)
parser.add_argument('-udk','--nt_updown', help='<Required> No of nts in upstream, downstream sequence files.', required=True)
parser.add_argument('-p','--peptide_file', help='<Required> Peptide sequence file.', required=True)
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

def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))
def align(chaine1, chaine2):
    res = ''
    for c1, c2 in zip(chaine1, chaine2):
        if c1==c2:
            res+='-'
        else:
            res+='*'
    return(res)

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

TX_peptide = {}
with f_open(args.peptide_file) as handle:
    for record in SeqIO.parse(handle, "fasta") :
        # Ensembl transcript IDs are 18characters
        TX_peptide[str(record.id)[:18]] = str(record.seq)

# CODON TABLE
COD_TO_AA = { # T
'TTT': 'Phe', 'TCT': 'Ser', 'TAT': 'Tyr', 'TGT': 'Cys', # TxT
'TTC': 'Phe', 'TCC': 'Ser', 'TAC': 'Tyr', 'TGC': 'Cys', # TxC
'TTA': 'Leu', 'TCA': 'Ser', 'TAA': '---', 'TGA': '---', # TxA
'TTG': 'Leu', 'TCG': 'Ser', 'TAG': '---', 'TGG': 'Trp', # TxG
# C
'CTT': 'Leu', 'CCT': 'Pro', 'CAT': 'His', 'CGT': 'Arg', # CxT
'CTC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', # CxC
'CTA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', # CxA
'CTG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', # CxG
# A
'ATT': 'Ile', 'ACT': 'Thr', 'AAT': 'Asn', 'AGT': 'Ser', # AxT
'ATC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', # AxC
'ATA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', # AxA
'ATG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', # AxG
# G
'GTT': 'Val', 'GCT': 'Ala', 'GAT': 'Asp', 'GGT': 'Gly', # GxT
'GTC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', # GxC
'GTA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', # GxA
'GTG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'  # GxG
}
AA_TO_SYM = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
'Trp': 'W', 'Asn': 'N', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Ala': 'A',
'Gly': 'G', 'Ile': 'I', 'Leu': 'L', 'His': 'H', 'Arg': 'R', 'Met': 'M',
'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', '---': '*'}
SYM_TO_COD = {sym:[cod for cod,aaa in COD_TO_AA.items() if aaa==aa] for aa,sym in AA_TO_SYM.items()}


# Parse input fasta that also includes CDS info / write into the output file
with f_open(args.out_fasta,'w') as out_f:
    with f_open(args.in_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta") :
            # Hold the current info
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

            #Ignore transcripts with CDS shorter than 30
            if cds[1][1]-cds[1][0]>28:
                # FIRST SOLVE THE FUCKED UP ANNOTATION in terms of CDS peptide identity
                # complete the cds length to multiple of 3, by adding or removing from downstream
                while (cds[1][1]-cds[1][0])%3 != 2:
                    if (1+cds[1][1]-cds[1][0]) < 3*len(TX_peptide[record.id[:18]]):
                        if utr3 is not None:
                            cds =  (cds[0], (cds[1][0], cds[1][1]+1))
                            utr3 = (utr3[0], (cds[1][1]+1, utr3[1][1]))
                        else:
                            seq = seq + downstream_seq[0]
                            downstream_seq = downstream_seq[1:]
                            TX_down[str(record.id)[:18]] =  downstream_seq
                            cds =  (cds[0], (cds[1][0], cds[1][1]+1))
                    else:
                        if utr3 is not None:
                            cds =  (cds[0], (cds[1][0], cds[1][1]-1))
                            utr3 = (utr3[0], (cds[1][1]+1, utr3[1][1]))
                        else:
                            downstream_seq = seq[-1] + downstream_seq
                            TX_down[str(record.id)[:18]] = downstream_seq
                            seq = seq[:-1]
                            cds =  (cds[0], (cds[1][0], cds[1][1]-1))

                # Look into cases where CDS translation and peptide sequence are highly dissimilar
                if str(Seq.translate(seq[cds[1][0]-1:cds[1][1]])) != TX_peptide[record.id[:18]]:
                    if hamming_distance(str(Seq.translate(seq[cds[1][0]-1:cds[1][1]])), TX_peptide[record.id[:18]]) > (3*len(TX_peptide[record.id[:18]])*0.1):
                        print(record.id)
                        print(utr5,cds,utr3)
                        upstream_seq = TX_up[str(record.id)[:18]]
                        downstream_seq = TX_down[str(record.id)[:18]]
                        full_seq = upstream_seq+seq+downstream_seq
                        peptide_seq = TX_peptide[str(record.id)[:18]]
                        cds_length = 3*len(peptide_seq)
                        # i is the number of nts we remove at 5'end each step
                        for i in range(4):
                            if i==3:
                                print("PROBLEM of not being able to find where CDS start")
                                exit(0)

                            k = str(Seq.translate(full_seq[i:])).find(peptide_seq[1:-1])
                            if k != -1:
                                print(i,k)
                                new_cds_start = ((k-1)*3) + i + 1 - int(args.nt_updown)
                                new_cds_end = new_cds_start + cds_length - 1
                                print(new_cds_start,new_cds_end)
                                # UPDATE CDSEND
                                if cds[1][1] < new_cds_end:
                                    if utr3 is not None:
                                        cds =  (cds[0], (cds[1][0], new_cds_end))
                                        utr3 = (utr3[0], (new_cds_end+1, utr3[1][1]))
                                    else:
                                        diff = new_cds_end - cds[1][1]
                                        seq = seq + downstream_seq[:diff]
                                        downstream_seq = downstream_seq[diff:]
                                        TX_down[str(record.id)[:18]] =  downstream_seq
                                        cds =  (cds[0], (cds[1][0], new_cds_end))
                                elif cds[1][1] > new_cds_end:
                                    if utr3 is not None:
                                        cds =  (cds[0], (cds[1][0], new_cds_end))
                                        utr3 = (utr3[0], (new_cds_end+1, utr3[1][1]))
                                    else:
                                        diff = cds[1][1] - new_cds_end
                                        downstream_seq = seq[-diff:] + downstream_seq
                                        TX_down[str(record.id)[:18]] =  downstream_seq
                                        seq = seq[:-diff]
                                        cds =  (cds[0], (cds[1][0], new_cds_end))
                                # UPDATE CDSSTART
                                if cds[1][0] > new_cds_start:
                                    if utr5 is not None:
                                        utr5 = (utr5[0], (utr5[1][0], new_cds_start-1))
                                        cds =  (cds[0], (new_cds_start, cds[1][1]))
                                    else:
                                        diff = cds[1][0] - new_cds_start
                                        seq = upstream_seq[-diff:] + seq
                                        upstream_seq = upstream_seq[:-diff]
                                        TX_up[str(record.id)[:18]] = upstream_seq
                                        cds = (cds[0], (1, cds[1][1]+diff))
                                        if utr3 is not None:
                                            utr3 = (utr3[0], (utr3[1][0]+diff, utr3[1][1]+diff))
                                elif cds[1][0] < new_cds_start:
                                    if utr5 is not None:
                                        utr5 = (utr5[0], (utr5[1][0], new_cds_start-1))
                                        cds =  (cds[0], (new_cds_start, cds[1][1]))
                                    else:
                                        diff = new_cds_start - cds[1][0]
                                        upstream_seq = upstream_seq + seq[:diff]
                                        TX_up[str(record.id)[:18]] = upstream_seq
                                        seq = seq[diff:]
                                        cds = (cds[0], (1, cds[1][1]-diff))
                                        if utr3 is not None:
                                            utr3 = (utr3[0], (utr3[1][0]-diff, utr3[1][1]-diff))
                                # update id_cols with updated cds utr5 and utr3
                                id_cols[cds[0]] = 'CDS:'+str(cds[1][0])+'-'+str(cds[1][1])
                                if utr5 is not None:
                                    id_cols[utr5[0]] = 'UTR5:'+str(utr5[1][0])+'-'+str(utr5[1][1])
                                if utr3 is not None:
                                    id_cols[utr3[0]] = 'UTR3:'+str(utr3[1][0])+'-'+str(utr3[1][1])

                                break
                        print(utr5,cds,utr3,'new')

                # ASSERT CDS TRANSLATION == PEPTIDE
                if str(Seq.translate(seq[cds[1][0]-1:cds[1][1]])) != TX_peptide[record.id[:18]]:
                    if hamming_distance(str(Seq.translate(seq[cds[1][0]-1:cds[1][1]])), TX_peptide[record.id[:18]]) > (3*len(TX_peptide[record.id[:18]])*0.1):
                        print(record.id)
                        print(len(seq))
                        print(TX_peptide[record.id[:18]][:10])
                        print(str(seq[cds[1][0]-1:cds[1][1]])[:30])
                        print(str(Seq.translate(seq[cds[1][0]-1:cds[1][1]]))[:10])
                        print("############### PEP #################")
                        print(len(TX_peptide[record.id[:18]]),3*len(TX_peptide[record.id[:18]]))
                        print(TX_peptide[record.id[:18]])
                        print("############### TRANSLATE #################")
                        print(str(Seq.translate(seq[cds[1][0]-1:cds[1][1]])))
                        print("############### ALIGN #################")
                        print(align(TX_peptide[record.id[:18]],str(Seq.translate(seq[cds[1][0]-1:cds[1][1]]))))

                        exit(0)


                ### DO the actual completion ###
                if len(id_cols) > 7 and ((args.include_labels==None) or (id_cols[7] in args.include_labels)) and ((args.exclude_labels==None) or (id_cols[7] not in args.exclude_labels)):
                    up = False
                    add = 0

                    if cds[1][0]<(int(args.k)+1):
                        add = (int(args.k)+1)-cds[1][0]
                        upstream_seq = TX_up[str(record.id)[:18]]
                        if utr3!=None:
                            utr3 = (utr3[0], ((utr3[1][0]+add),(utr3[1][1]+add)))
                        if cds!=None:
                            cds = (cds[0], ((int(args.k)+1), (cds[1][1]+add)))
                        if utr5!=None:
                            utr5 = (utr5[0],(1,int(args.k)))
                        else:
                            # update id_cols & cds & utr5 & utr3
                            id_cols = id_cols[:cds[0]]+['UTR5:1-'+args.k]+id_cols[cds[0]:]
                            utr5 = (cds[0], (1,int(args.k)))
                            cds = (cds[0]+1, cds[1])
                            if utr3!=None: # only the position in id_cols
                                utr3 = (utr3[0]+1, utr3[1])

                        if id_cols[-1] == '':
                            id_cols[-1] = 'up'+str(add)
                        else:
                            id_cols.append = 'up'+str(add)

                        seq = upstream_seq[-add:] + seq
                        up = True

                    if (len(seq)-cds[1][1]) < int(args.k):
                        downstream_seq = TX_down[str(record.id)[:18]]
                        add = (int(args.k)-(len(seq)-cds[1][1]))
                        if utr3!=None:
                            if cds[1][1]+1 != utr3[1][0]:
                                print('PROBLEM !!!')
                                print(record.id)
                                print(utr5,cds,utr3)
                                exit()
                            utr3 = (utr3[0],(utr3[1][0], utr3[1][0]+(int(args.k)-1)))
                        else:
                            id_cols = id_cols[:(cds[0]+1)]+['UTR3:'+str(cds[1][1]+1)+'-'+str(cds[1][1]+int(args.k))]+ id_cols[(cds[0]+1):]
                            utr3 = (cds[0]+1, (cds[1][1]+1, cds[1][1]+int(args.k)))

                        seq = seq + downstream_seq[:add]
                        if up:
                            id_cols[-1] += 'down'+str(add)
                        elif id_cols[-1] == '':
                            id_cols[-1] = 'down'+str(add)
                        else:
                            id_cols.append = 'down'+str(add)

                    # update id_cols with updated cds utr5 and utr3
                    id_cols[cds[0]] = 'CDS:'+str(cds[1][0])+'-'+str(cds[1][1])
                    if utr5 is not None:
                        id_cols[utr5[0]] = 'UTR5:'+str(utr5[1][0])+'-'+str(utr5[1][1])
                    if utr3 is not None:
                        id_cols[utr3[0]] = 'UTR3:'+str(utr3[1][0])+'-'+str(utr3[1][1])

                    out_f.write('>'+str('|'.join(id_cols))+'|\n'+seq+'\n')

                else:
                    print('PROBLEM !!!')
                    print(record.id)
                    exit()
