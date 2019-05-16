
# Script to seperate transcripts from given gencode fasta file
# Takes 4 arguments (see below)

import sys,os,gzip
from Bio import SeqIO

def file_open(x, mode='rf'):
	return gzip.open(x,mode) if x.endswith(".gz") else open(x,mode)

print sys.argv
if len(sys.argv)!=4:
	sys.stderr.write('Problem in given arguments\n')
else:
	given_fasta=sys.argv[1]
	rRNA_fasta=sys.argv[2]
	others_fasta=sys.argv[3]
	#tRNA_fasta=sys.argv[4]

	rRNA_keys = set(['Mt_rRNA','rRNA'])
	#tRNA_keys = set(['Mt_tRNA','tRNA'])
	rRNAs=[]
	#tRNAs=[]
	others=[]

	with file_open(given_fasta, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			gt = ""
			cols = record.id.split('|')
			if cols[7] in rRNA_keys:
				rRNAs.append(record)
			#elif cols[7] in tRNA_keys:
			#	tRNAs.append(record)
			else:
				others.append(record)

		SeqIO.write(rRNAs,rRNA_fasta,"fasta")
		SeqIO.write(others,others_fasta,"fasta")
		#SeqIO.write(tRNAs,tRNA_fasta,"fasta")


