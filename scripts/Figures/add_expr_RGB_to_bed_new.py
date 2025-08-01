# create 4 output bed files with corresponding rgb values for expression levels

# !/usr/bin/env python3

fhand = open("tama2geneid_transcriptome.filtered.bed")
tx_ids = open("expression_colours_txs.txt")

samples=['fetal','adult','hippo','caudate']
for sample in samples:
    fout=open(sample+"_coloured_by_expr.bed","w")
    txDict = dict()
    for line in tx_ids:
        lineList = line.split('\t')
        txDict[lineList[0]] = lineList[samples.index(sample)+1].strip('\n') # use index of sample in samples list to get corresponding column in tx_ids file
    for line in fhand:
        if not line.startswith('#'):
                chr,start,end,name,score,strand,a,b,rgb,blockCount,blocksizes,blockstarts = line.split('\t')
                tx_id=name.split(';')[1]
                for key in txDict:
                    if tx_id == key:
                            line = line.replace(rgb,txDict[key])
                            fout.write(line)
    fhand.seek(0)
    tx_ids.seek(0)