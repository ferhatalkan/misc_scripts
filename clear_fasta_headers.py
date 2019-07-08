#!/usr/bin/env python3

import sys

from Bio import SeqIO, Seq, SeqRecord
from Bio.Alphabet import generic_dna, generic_rna

with open(sys.argv[1]) as handle:
    SeqIO.write([SeqRecord.SeqRecord(record.seq, id=record.id.split('|')[0],description='') for record in SeqIO.parse(handle, "fasta")], sys.argv[2], "fasta")
