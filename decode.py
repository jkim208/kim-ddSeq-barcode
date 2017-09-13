import re
import editdistance # Find levenshtein distance
import distance  # Find hamming distance
from Bio import SeqIO
from collections import Counter

# Use a small subset of random records to determine the true linker 
linker1s = []
linker2s = []
for record in SeqIO.parse("example.fastq", "fastq"):
    seqString = "%s" % (record.seq)
    matchObj = re.match(r'.*([ACGT]{6})([ACGT]{15})([ACGT]{6})([ACGT]{15})([ACGT]{6})ACG([AGCT]{8})GACT', seqString)
    linker1s.append(matchObj.group(2))  # set of linker 1s
    linker2s.append(matchObj.group(4))  # set of linker 2s

linker1True = Counter(linker1s).most_common(1)[0][0]
linker2True = Counter(linker2s).most_common(1)[0][0]


# read in all 96 possible cell barcode blocks in 'allBarcodeBlocks'
with open('barcodeBlocks.txt') as barcodeBlocksFile:
    allBarcodeBlocks = barcodeBlocksFile.read().splitlines()


# Function "correct_bc_blocks" takes each barcode_block from a sequence, checks the hamming distances
# between the variable block and all 96 possible blocks, and makes corrections depending on the smallest
# ED(edit distance) found. Barcode block is unchanged if smallest ED is 2.
# This function is used in the parsing loop below.

def correct_bc_blocks(barcode_block):
    hamming_dict = {}
    print('In function "correct_bc_blocks" for ' + barcode_block)
    # determine hamming distances between barcode block and all 96 possible blocks
    for i in range(95):
        hamming_dict[i] = (distance.hamming(barcode_block, allBarcodeBlocks[i]))

    # find the blocks with the smallest hamming distances
    shd_index = min(hamming_dict, key=hamming_dict.get)
    print('Smallest hamming dist:', hamming_dict[shd_index])  # value of smallest hamming distance

    if hamming_dict[shd_index] == 1:
        print('Fix barcode to the one which is 1 ED apart')
        barcode_block = allBarcodeBlocks[shd_index]
    print('Adjusted barcode is:' + allBarcodeBlocks[shd_index])
    return barcode_block


linker1s = []  # redefine linker lists for reuse
linker2s = []
# Loop through all records and follow decoding algorithm
# For linkers, allow only 1 edit distance
for record in SeqIO.parse("example.fastq", "fastq"):
    print(record.id)
    seqString = "%s" % (record.seq)
    matchObj = re.match(r'.*([ACGT]{6})([ACGT]{15})([ACGT]{6})([ACGT]{15})([ACGT]{6})ACG([AGCT]{8})GACT', seqString)
    linker1s.append(matchObj.group(2))
    linker2s.append(matchObj.group(4))
    if matchObj:
        # print "matchObj.group() : ", matchObj.group()
        # print "BC1       : ", matchObj.group(1)
        # print("Linker1   : ", matchObj.group(2))

        # evaluate edit distance. If >1 ED, skip this sequence
        print('Edit distance for linker1: ', editdistance.eval(linker1True, matchObj.group(2)))
        if editdistance.eval(linker1True, matchObj.group(2)) > 1:
            print('linker 1 has >1 edit distance')
            continue
        # print "BC2       : ", matchObj.group(3)
        # print("Linker2   : ", matchObj.group(4))
        print('Edit distance for linker2: ', editdistance.eval(linker2True, matchObj.group(4)))
        if editdistance.eval(linker2True, matchObj.group(4)) > 1:
            print('linker 2 has >1 edit distance')
            continue
        # print "BC3       : ", matchObj.group(5)
        # print "UMI       : ", matchObj.group(6)
        bc1 = correct_bc_blocks(matchObj.group(1))
        bc2 = correct_bc_blocks(matchObj.group(3))
        bc3 = correct_bc_blocks(matchObj.group(5))
        bcFull = bc1 + bc2 + bc3
        print('Full cell barcode is ' + bcFull + '\n')

    else:
        print("No match")
