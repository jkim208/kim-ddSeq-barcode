#!/usr/bin/env python
import argparse  # command line options
import regex  # regular expressions with edit distance functions
import distance  # Find hamming distance

# Function "get_all_barcode_blocks" reads from a file containing all possible barcode blocks.
# Returns a list with these barcode blocks.
def get_all_barcode_blocks(barcode_blocks_file):
    # read in all 96 possible cell barcode blocks in 'all_barcode_blocks'
    with open(barcode_blocks_file, 'r') as barcode_blocks_fh:
        all_barcode_blocks = barcode_blocks_fh.read().splitlines()
        return all_barcode_blocks

# Function "correct_bc_blocks" takes each barcode_block from a sequence, checks the hamming distances
# between the variable block and all 96 possible blocks, and makes corrections depending on the smallest
# ED(edit distance) found. Barcode block is unchanged if smallest ED is 2.
def correct_bc_blocks(all_barcode_blocks, barcode_block):
    hamming_dict = {}

    # determine hamming distances between barcode block and all 96 possible blocks
    for i in range(95):
        hamming_dict[i] = (distance.hamming(barcode_block, all_barcode_blocks[i]))

    # find the blocks with the smallest hamming distances
    shd_index = min(hamming_dict, key=hamming_dict.get)

    if hamming_dict[shd_index] == 1:
    #    print('Fix barcode to the one which is 1 ED apart')
        barcode_block = all_barcode_blocks[shd_index]
    #elif hamming_dict[shd_index] == 2:
    #    print('Smallest ED is 2. Barcode will be not be changed.')
    #elif hamming_dict[shd_index] > 2:
    #    print('Smallest ED is > 2.')

    # print('Adjusted barcode is:' + all_barcode_blocks[shd_index])
    return barcode_block

# Function "demultiplex_barcodes" accounts for edit distance while extracting barcodes
# Includes the correct_bc_blocks function in order to return full barcode
def demultiplex_barcodes(all_records, all_barcode_blocks):
    originalSAM = open(all_records, 'r')
    barcodedRead2File = open('writeSamTest.sam', 'w')

    # write first two lines of sam file (header lines) to new file.
    # start the loop once the pointer is on the actual records
    barcodedRead2File.write(originalSAM.readline())
    barcodedRead2File.write(originalSAM.readline())
    # declare variable to hold one full barcode in one scope higher than the loop
    # so that the barcode from read1 is 'saved' for appending to read2
    cell_bc = ''
    umi = ''
    match_obj1 = ''
    # loop through read1 records and apply decoding algorithm
    for count, line in enumerate(originalSAM, start=0):
        # every even line refers to read1. Decode read1 for barcodes. Do not write read1 to new SAM file
        if count % 2 == 0:
            # match to where the sequence should be
            match_obj1 = regex.match(r".*\t([NACGT]*)\t.*\t.*", line)
            if match_obj1:
                seq_string = match_obj1.group(1)

                # allow up to 1 edit distance(ED) when matching for Linkers
                # ?e flag forces best matches (as few EDs as possible)
                linker1 = regex.findall("(?e)(TAGCCATCGCATTGC){e<=1}", seq_string)
                linker2 = regex.findall("(?e)(TACCTCTGAGCTGAA){e<=1}", seq_string)

                if linker1 and linker2:
                    # match blocks accordingly with linkers
                    # consider potential insertions/deletions in barcodes/anchors
                    # only keep reads where ACG and GACT anchors are not mutated
                    # main unaccounted case is if insertions occur before/after barcode and before ACG
                    regex_search = r".*([ACGT]{6})" + str(linker1[0]) + r"([ACGT]{6})" \
                                   + str(linker2[0]) + r"([ACGT]{6})ACG([AGCT]{8}).?GACT"
                    # match to find the barcodes
                    match_obj2 = regex.match(regex_search, seq_string)

                    if match_obj2:
                        bc1 = correct_bc_blocks(all_barcode_blocks, match_obj2.group(1))
                        bc2 = correct_bc_blocks(all_barcode_blocks, match_obj2.group(2))
                        bc3 = correct_bc_blocks(all_barcode_blocks, match_obj2.group(3))
                        umi = match_obj2.group(4)
                        cell_bc = bc1 + bc2 + bc3

                    else:
                        #print("Regex match failed. Either sequence was mutated beyond acceptable measures"
                              #"or the wrong read was processed")
                        match_obj1 = ''

                else:
                    #print("Linkers could not be found in good condition (ED > 1)")
                    match_obj1 = ''

            else:
                #print('Did not match to a SAM record')
                match_obj1 = ''

        # every odd line is read2. Read2 will have its sequence appended by the previous read1 barcode
        if count % 2 == 1 and match_obj1:
                read2 = line.rstrip().split('\t')

                # find where the sequence is in the list. Then, append the barcode
                read2[9] = read2[9] + cell_bc + umi
                # also, add tags to the SAM/BAM file line to represent what the barcodes are
                read2.append('XC:Z:' + cell_bc)
                read2.append('XM:Z:' + umi)

                # rejoin the split line for writing to the new Sam file
                barcodedRead2 = '\t'.join(read2)
                barcodedRead2File.write(barcodedRead2 + '\n')

    originalSAM.close()
    barcodedRead2File.close()
    return

def main():
    # grab SAM filename/path from command line arguments
    parser = argparse.ArgumentParser(description='Process paired reads from SAM file by extracting '
                                                 'barcodes from read 1 and tagging them onto read 2')
    parser.add_argument('filename1', help='(.sam file)')
    parser.add_argument('filename2', help='(file containing barcode blocks)')
    args = parser.parse_args()

    # obtain all possible barcode block combinations
    barcode_blocks_file = args.filename2
    all_barcode_blocks = get_all_barcode_blocks(barcode_blocks_file)

    # construct full cell barcodes from every sequence record. Supply the records in SAM format
    sam_records_file = args.filename1
    demultiplex_barcodes(sam_records_file, all_barcode_blocks)

if __name__ == "__main__":
    main()
