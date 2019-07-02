#!/usr/bin/env python3
import sys

## pass the fasta file with CDS,UTR information in headers as sys.argv[1] and output prefix as sys.argv[2]
genome = open(sys.argv[2]+'_CustomGenome.fa','w')
annot = open(sys.argv[2]+'_CustomGenome_annot.gff3','w')
alt_annot = open(sys.argv[2]+'_CustomGenome_annot_alt.gff3','w')
alias = open(sys.argv[2]+'_CustomGenome.alias','w')

annot.write('##gff-version 3\n')
annot.write('#description: CUSTOM TX-specific annotation of the mouse genome (GRCm38), version M21 (Ensembl 96)\n')
annot.write('#provider: Faller Lab based on GENCODE\n')
annot.write('#contact: f.alkan@nki.nl and gencode-help@ebi.ac.uk\n')
annot.write('#format: gff3\n')
annot.write('#date: 2019-07-01 2019-03-27\n')

with open(sys.argv[1]) as inf:
    for line in inf:
        if line[0]=='>':
            cols = line[1:].rstrip().split('|')
            tid,gid,gid2,tid2,ufTid,ufGid,wrong_l,utr5,cds,utr3 = cols[:10]
            genome.write('>'+ufTid+'\n')
            alias.write(line[1:].rstrip()+'\t'+ufTid+'\n')
            alias.write(ufTid+'\t'+ufTid+'\n')

            utr5spos,utr5epos = utr5.split(':')[1].split('-')
            utr3spos,utr3epos = utr3.split(':')[1].split('-')
            CDSspos,CDSepos = cds.split(':')[1].split('-')

            annot.write('\t'.join([ufTid,'CUSTOM','gene',utr5spos,utr3epos,'.','+','.',
                        ';'.join(['ID=TX:'+ufTid,'gene=TX_'+ufGid,'gene_id=TX_'+ufTid,'gene_type=ProteinCoding','EnsemblGeneID='+gid,'EnsemblTxID='+tid])])+'\n')
            annot.write('\t'.join([ufTid,'CUSTOM','five_prime_UTR',utr5spos,utr5epos,'.','+','.',
                        ';'.join(['ID=five_prime:'+ufTid,'Parent=TX_'+ufTid,'gene='+ufGid,'gene_id=TX_'+ufTid,'gene_type=ProteinCoding','EnsemblGeneID='+gid,'EnsemblTxID='+tid])])+'\n')
            annot.write('\t'.join([ufTid,'CUSTOM','CDS',CDSspos,CDSepos,'.','+','.',
                        ';'.join(['ID=CDS:'+ufTid,'Parent=TX_'+ufTid,'gene='+ufGid,'gene_id=TX_'+ufTid,'gene_type=ProteinCoding','EnsemblGeneID='+gid,'EnsemblTxID='+tid])])+'\n')
            annot.write('\t'.join([ufTid,'CUSTOM','three_prime_UTR',utr3spos,utr3epos,'.','+','.',
                        ';'.join(['ID=three_prime:'+ufTid,'Parent=TX_'+ufTid,'gene='+ufGid,'gene_id=TX_'+ufTid,'gene_type=ProteinCoding','EnsemblGeneID='+gid,'EnsemblTxID='+tid])])+'\n')

        else:
            genome.write(line)
genome.close()
annot.close()
alt_annot.close()
alias.close()
