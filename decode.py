#!/usr/bin/env python
import re
import editdistance # Find levenshtein distance
import distance  # Find hamming distance
from Bio import SeqIO
from collections import Counter

# Function "get_all_barcode_blocks" reads from a file containing all possible barcode blocks.
# Returns a list with these barcode blocks.
def get_all_barcode_blocks():
    # read in all 96 possible cell barcode blocks in 'all_barcode_blocks'
    with open('barcodeBlocks.txt', 'r') as barcode_blocks_file:
        all_barcode_blocks = barcode_blocks_file.read().splitlines()
        return all_barcode_blocks

# Function "find_true_linkers" uses a small subset of random records to find the
# true linker sequences.
def find_true_linkers(sample_records):
    linker1s = []
    linker2s = []
    for record in SeqIO.parse(sample_records, "fastq"):
        seq_string = "%s" % record.seq
        match_obj = re.match(r'.*([ACGT]{6})([ACGT]{15})([ACGT]{6})([ACGT]{15})([ACGT]{6})ACG([AGCT]{8})GACT',
                             seq_string)
        linker1s.append(match_obj.group(2))  # set of linker 1s
        linker2s.append(match_obj.group(4))  # set of linker 2s

    linker1_true = Counter(linker1s).most_common(1)[0][0]
    linker2_true = Counter(linker2s).most_common(1)[0][0]
    return linker1_true, linker2_true

# Function "correct_bc_blocks" takes each barcode_block from a sequence, checks the hamming distances
# between the variable block and all 96 possible blocks, and makes corrections depending on the smallest
# ED(edit distance) found. Barcode block is unchanged if smallest ED is 2.
def correct_bc_blocks(all_barcode_blocks, barcode_block):
    hamming_dict = {}
    print('In function "correct_bc_blocks" for ' + barcode_block)
    # determine hamming distances between barcode block and all 96 possible blocks
    for i in range(95):
        hamming_dict[i] = (distance.hamming(barcode_block, all_barcode_blocks[i]))

    # find the blocks with the smallest hamming distances
    shd_index = min(hamming_dict, key=hamming_dict.get)
    print('Smallest hamming dist:', hamming_dict[shd_index])  # value of smallest hamming distance

    if hamming_dict[shd_index] == 1:
        print('Fix barcode to the one which is 1 ED apart')
        barcode_block = all_barcode_blocks[shd_index]
    print('Adjusted barcode is:' + all_barcode_blocks[shd_index])
    return barcode_block


def demultiplex_barcodes(all_records, all_barcode_blocks, linker1_true, linker2_true):
    linker1s = []  # define linker and barcode lists
    linker2s = []
    full_barcodes = []
    full_barcodes_file = open('fullBarcodes.txt', 'w')  # store full barcodes in new file
    # Loop through all records and follow decoding algorithm
    # For linkers, allow only 1 edit distance
    for record in SeqIO.parse(all_records, "fastq"):
        print(record.id)
        seq_string = "%s" % (record.seq)
        match_obj = re.match(r'.*([ACGT]{6})([ACGT]{15})([ACGT]{6})([ACGT]{15})([ACGT]{6})ACG([AGCT]{8})GACT',
                             seq_string)
        linker1s.append(match_obj.group(2))
        linker2s.append(match_obj.group(4))
        if match_obj:
            # print "match_obj.group() : ", match_obj.group()
            # print "BC1       : ", match_obj.group(1)
            # print("Linker1   : ", match_obj.group(2))

            # evaluate edit distance. If >1 ED, skip this sequence
            print('Edit distance for linker1: ', editdistance.eval(linker1_true, match_obj.group(2)))
            if editdistance.eval(linker1_true, match_obj.group(2)) > 1:
                print('linker 1 has >1 edit distance')
                continue
            # print "BC2       : ", match_obj.group(3)
            # print("Linker2   : ", match_obj.group(4))
            print('Edit distance for linker2: ', editdistance.eval(linker2_true, match_obj.group(4)))
            if editdistance.eval(linker2_true, match_obj.group(4)) > 1:
                print('linker 2 has >1 edit distance')
                continue
            # print "BC3       : ", match_obj.group(5)
            # print "UMI       : ", match_obj.group(6)
            bc1 = correct_bc_blocks(all_barcode_blocks, match_obj.group(1))
            bc2 = correct_bc_blocks(all_barcode_blocks, match_obj.group(3))
            bc3 = correct_bc_blocks(all_barcode_blocks, match_obj.group(5))
            umi = match_obj.group(6)
            print('UMI is ' + umi)
            bc_full_cell = bc1 + bc2 + bc3 + umi
            full_barcodes.append(bc_full_cell)
            print('Full cell barcode is ' + bc_full_cell + '\n')
            full_barcodes_file.write("%s\n" % bc_full_cell)

        else:
            print("No match")

    full_barcodes_file.close()
    return full_barcodes


def main():
    # obtain all possible barcode block combinations
    all_barcode_blocks = get_all_barcode_blocks()

    # determine the true linker. Supply the function a subset of records in fastq format.
    sample_records = "example.fastq"
    (linker1_true, linker2_true) = find_true_linkers(sample_records)

    # construct full cell barcodes from every sequence records. Supply the records in fastq format.
    all_records = "example.fastq"
    full_barcodes = demultiplex_barcodes(all_records, all_barcode_blocks, linker1_true, linker2_true)
    # i = 0
    # for record in SeqIO.parse("example2.fastq", "fastq"):
    #     print(record.seq + full_barcodes[i])
    #     i = i + 1
    #     if i == len(full_barcodes):
    #         break


if __name__ == "__main__":
    main()