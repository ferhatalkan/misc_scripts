#!/usr/bin/env python3

import sys

class Transcript:
    def  __init__(self, data, inds):
        self.id = data[inds["Transcript stable ID"]]
        self.length =  int(data[inds["Transcript length (including UTRs and CDS)"]])
        self.type = data[inds["Transcript type"]]
        self.name = data[inds["Transcript name"]]
        if inds["Protein stable ID"] < len(data):
            self.pid = data[inds["Protein stable ID"]]
        else:
            self.pid = None
        self.otherIDs = {}
        for otherID in ["RefSeq mRNA ID", "RefSeq mRNA predicted ID", "CCDS ID","APPRIS annotation", "Transcript support level (TSL)"]:
            if inds[otherID] < len(data):
                if data[inds[otherID]] != '':
                    self.otherIDs[otherID] = set([data[inds[otherID]]])
    def add_otherIDS(self,tag,oid):
        if oid!='':
            if tag in self.otherIDs:
                self.otherIDs[tag].add(oid)
            else:
                self.otherIDs[tag] = set([oid])

class Gene:
    def  __init__(self, data, inds):
        self.id = data[inds["Gene stable ID"]]
        self.name = data[inds["Gene name"]]
        self.TXs = {data[inds["Transcript stable ID"]]:Transcript(data, inds)}
    def add_TX(self, data, inds):
        TXid = data[inds["Transcript stable ID"]]
        if TXid in self.TXs.keys():
            for otherID in ["RefSeq mRNA ID", "RefSeq mRNA predicted ID", "CCDS ID","APPRIS annotation", "Transcript support level (TSL)"]:
                if inds[otherID] < len(data):
                    self.TXs[TXid].add_otherIDS(otherID, data[inds[otherID]])
        else:
            self.TXs[TXid] = Transcript(data, inds)

gene_dic = {}
with open(sys.argv[1]) as in_f:
    colnames = {}
    for line in in_f:
        if len(colnames)==0:
            cols = line.rstrip().split('\t')
            for i in range(len(cols)):
                colnames[cols[i]] = i
        else:
            cols = line.rstrip().split('\t')
            if cols[colnames["Gene stable ID"]] in gene_dic:
                gene_dic[cols[colnames["Gene stable ID"]]].add_TX(cols,colnames)
            else:
                gene = Gene(cols,colnames)
                gene_dic[gene.id] = gene
            #
            # if cols[colnames["Gene stable ID"]] == "ENSMUSG00000078695":
            #     print(line)
            #     for TX in gene_dic["ENSMUSG00000078695"].TXs.values():
            #         print(TX.otherIDs)

print(len(gene_dic),"genes")

n201 = 0
ntot = 0
npro = 0
with open(sys.argv[2],'w') as outf:
    for gene_id, gene in gene_dic.items():
        selectedTX = None
        tx201 = None

        for TXid,transcript in gene.TXs.items():
            if transcript.type=="protein_coding":
                if transcript.name.endswith("-201"):
                    tx201 = transcript
                npro += 1
                if selectedTX is None:
                    selectedTX=transcript
                else:
                    if len(selectedTX.otherIDs.keys())==0:
                        if transcript.length > selectedTX.length:
                            selectedTX = transcript
                            continue
                    for tag in ["APPRIS annotation", "CCDS ID", "RefSeq mRNA ID", "RefSeq mRNA predicted ID", "Protein stable ID"]:

                        if tag == "APPRIS annotation":
                            selected_APPRIS = None if tag not in selectedTX.otherIDs else list(selectedTX.otherIDs[tag])[0]
                            tx_APPRIS = None if tag not in transcript.otherIDs else list(transcript.otherIDs[tag])[0]
                            if tx_APPRIS is None:
                                if selected_APPRIS is not None:
                                    if selected_APPRIS.startswith("alternative"):
                                        continue
                                    if selected_APPRIS.startswith("principal"):
                                        break
                                else:
                                    continue
                            elif selected_APPRIS is None:
                                if tx_APPRIS is not None:
                                    if tx_APPRIS.startswith("alternative"):
                                        continue
                                    if tx_APPRIS.startswith("principal"):
                                        selectedTX=transcript
                                        break
                            else:
                                if tx_APPRIS.startswith("principal") and selected_APPRIS.startswith("alternative"):
                                    selectedTX=transcript
                                    break
                                elif tx_APPRIS.startswith("alternative") and selected_APPRIS.startswith("principal"):
                                    break

                        if tag in selectedTX.otherIDs and tag not in transcript.otherIDs:
                            break

                        if tag in transcript.otherIDs:
                            if tag in selectedTX.otherIDs:
                                if transcript.length > selectedTX.length:
                                    selectedTX=transcript
                                    break
                            else:
                                selectedTX=transcript
                                break

        if selectedTX is not None:
            outf.write(selectedTX.id+'\n')
            ntot+=1
            if selectedTX.name.endswith("-201"):
                n201+=1
            else:
                if tx201 is not None:
                    print(gene_id,selectedTX.name,selectedTX.length,tx201.name,tx201.length)

print(npro, "npro", ntot, "ntot", n201,"n201")
