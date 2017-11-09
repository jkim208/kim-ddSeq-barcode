#!/usr/bin/env python3
import argparse  # command line options
import regex  # regular expressions with edit distance functions
import pyximport  # install Cython to access pyximport
pyximport.install(build_in_temp=False)
from editDistance import edit_distance
import sys
import distance

# Updates:
# ACGGAC must be correctly positioned. There MUST be >=1 base after the anchor
# Any indels in linker 1 + 2 are not accepted.
# No quality score of < 10 is allowed
# #REMOVED# Read 1 sequence MUST end in GAC[T]+. Used distance.levenshtein to track missing T's.
# Phase blocks accept indels/substitutions (doing so prevents room for future mutations in that sequence)
# The sequence tail (ACGGACT) can accept mutations (up to 1 ED assuming no other source of EDs)

# Summary:
# parseBarcodes.py was written for the BioRad ddSeq procedure (scRNA).
# This program applies a decoding algorithm recommended by Illumina to extract, parse and filter
# barcodes in read 1 of a paired set. The barcode information is then tagged onto the SAM file record
# of read 2.

# Decoding logic:
# The bar code component has 6 read structures, 5 of which have distinct phase blocks. Room for mutations
# should be minimized to 1 edit distance (insertion, deletion, substitution) and indels will be only allowed to
# occur in the linkers (when it does occur in the linker, a counter will keep track of the change for the
# downstream base positions to assure that only one indel exists per read). An insertion or deletion at
# the bar code will lead to a dropped read. Any substitution mutations on the bar codes will be corrected
# if possible. The ACG and GAC anchors are used to catch non-linker insertions/deletions and confirm
# correct positioning of the read.


# Keep track of success and failures
bad_phase = 0
bad_block = 0
low_quality = 0
bad_linker = 0
matches = 0

def append_barcode(line, cell_bc, umi):
    # Function 6 "append_barcode" takes the barcode and adds it to the record as a separate tag
    line = "%s\t%s\t%s" % (line.rstrip(), 'XC:Z:' + cell_bc, 'XM:Z:' + umi)
    return line


def check_bc_quality(q_seq, bc_index):
    # Function 5 'check_bc_quality' accepts a quality sequence string, a barcode index, and current low
    # quality count. Via the barcode index, the function checks if the corresponding quality string
    # indicates any low quality bases. If any are detected, loq_q_count is incremented by 1.
    # Function 4 is used in function 3 (demultiplexing) and returns the new low_q_count.
    low_q_count = 0

    for element in q_seq[bc_index:bc_index + 6]:
        # Break out of loop immediately if low_q_count threshold has been passed
        if low_q_count > 0:
            break

        # Check the ASCII code if the base quality is low (q-score 10 = ASCII score 43)
        if ord(element) < 43:
            low_q_count += 1

    return low_q_count


def correct_bc_blocks(ref_barcode_blocks, barcode_block):
    # Function 4 "correct_bc_blocks" takes each barcode_block from a sequence, checks the hamming distances
    # between the variable block and all 96 possible blocks, and makes corrections depending on the smallest
    # ED(edit distance) found. Barcode block is unchanged if smallest ED is 2.

    # Barcode blocks must be 6 bases long to be valid
    if len(barcode_block) != 6:
        barcode_block = None
        return barcode_block

    # start above acceptable threshold
    lowest_hamming = 6
    lh_reference_block = ''

    # end function immediately if barcode block matches list of known barcode blocks
    if barcode_block in ref_barcode_blocks:
        return barcode_block

    # determine hamming distances between barcode block and all 96 possible blocks
    for reference_block in ref_barcode_blocks:

        # use Cython to rapidly conduct calculation. .pyx script will perform type conversion
        hamming_dist = (edit_distance(barcode_block, reference_block))

        if hamming_dist < lowest_hamming:
            lowest_hamming = hamming_dist
            lh_reference_block = reference_block
        if lowest_hamming == 1:
            # print('Fix barcode to the one which is 1 ED apart')
            return lh_reference_block

    # find the blocks with the smallest hamming distances. If that ED is >1, skip the read
    barcode_block = None

    return barcode_block


def demultiplex(match_obj1, mod, linker1, linker2, ref_barcode_blocks, read1):
    global bad_phase
    global bad_block
    global low_quality
    global bad_linker
    global matches
    linker_ed = 0

    bc1 = match_obj1[0 + mod:6 + mod]

    bc2 = match_obj1[linker1.end(1): linker2.start(1)]

    bc3 = match_obj1[linker2.end(1): linker2.end(1) + 6]
    ACGGAC = match_obj1[linker2.end(1) + 6:linker2.end(1) + 9] + \
             match_obj1[linker2.end(1) + 17:linker2.end(1) + 20]
    postBase = match_obj1[linker2.end(1) + 20:]  # need a base after the GAC anchor

    if edit_distance(ACGGAC, 'ACGGAC') > 0:
        # mutations in these two anchors are not tolerated
        bad_block += 1
        return None, None, None

    if not postBase:
        bad_block += 1
        return None, None, None

    umi = match_obj1[linker2.end(1) + 9:linker2.end(1) + 17]

    bc1_n = correct_bc_blocks(ref_barcode_blocks, bc1)
    bc2_n = correct_bc_blocks(ref_barcode_blocks, bc2)
    if bc2_n and bc1_n:
        bc3_n = correct_bc_blocks(ref_barcode_blocks, bc3)
    else:
        bc3_n = None

    if bc3_n:
        # for bc1, beginning is at index 0.
        # for bc2, beginning is after at least bc1 (6b) and linker1 (15b)
        # for bc3, beginning is after at least bc2 (12b) and linker2 (30b)
        # 68 is the length of the sequence. Hard code it for speed. Otherwise, remove.
        bc1_index = match_obj1.find(bc1, 0, 20)
        bc2_index = match_obj1.find(bc2, 20, 41)
        bc3_index = match_obj1.find(bc3, 41, 55)
        umi_index = match_obj1.find(umi, 48, 68)

        # count number of low quality bases.
        low_quality_count = 0
        # break the loop and remove the read combo if count is > 1
        low_quality_count += check_bc_quality(read1[10], bc1_index)
        low_quality_count += check_bc_quality(read1[10], bc2_index)
        low_quality_count += check_bc_quality(read1[10], bc3_index)
        low_quality_count += check_bc_quality(read1[10], umi_index)

        # No low quality barcode bases allowed
        if low_quality_count == 0:
            cell_bc = bc1_n + bc2_n + bc3_n
            matches += 1
            return match_obj1, cell_bc, umi

        else:
            low_quality += 1
            pass
    else:
        bad_block += 1
        pass

    return None, None, None


def extract_barcode(line, ref_barcode_blocks):
    global bad_phase
    global bad_linker

    # Function 3: "extract_barcode" uses regex to extract barcode blocks and return complete barcodes
    # split read 1 to extract relevant parameters
    read1 = line.rstrip().split('\t')

    # match to where the sequence should be
    match_obj1 = read1[9]
    if 'N' in match_obj1:
        return None, None, None

    phase_blocks = ['', 'A','CT','GCA','TGCG','ATCGA']

    if match_obj1:
        # match blocks accordingly with linkers, leaving room for 1 edit distance
        # only keep reads where ACG and GACT anchors are not mutated
        # main unaccounted case is if insertions occur before/after barcode and before ACG
        linker1 = regex.search(r"(TAGCCATCGCATTGC){s<=1}", match_obj1)
        linker2 = regex.search(r"(?er)(TACCTCTGAGCTGAA){s<=1}", match_obj1)

        if linker1 and linker2:

            pb = match_obj1[0:linker1.start(1) - 6]

            if pb in phase_blocks:  # same as if ED = 0 between pb and corresponding reference block
                mod = len(pb)
                return demultiplex(match_obj1, mod, linker1, linker2, ref_barcode_blocks, read1)

            lowest_dist = 2
            for phase_block in phase_blocks:

                dist = distance.levenshtein(pb, phase_block)

                if dist < lowest_dist:
                    lowest_dist = dist
                    pb_reference = phase_block

            if lowest_dist <= 2:  # ??? Might be too much
                mod = len(pb)
                return demultiplex(match_obj1, mod, linker1, linker2, ref_barcode_blocks, read1)
            else:  # lowest levenshtein distance is > 1
                bad_phase += 1
                pass
        else:
            bad_linker += 1
            pass
    else:
        # print('Did not match to a SAM record')
        pass

    umi = None
    match_obj1 = None
    cell_bc = None

    return match_obj1, cell_bc, umi


def read_and_write_sam(all_records, ref_barcode_blocks, output):
    # Function 2 "read_and_write_sam" accounts for edit distance while extracting barcodes
    # Includes the correct_bc_blocks function in order to return full barcode

    try:
        originalSAM = open(all_records, 'r')
    except IOError:
        print("Could not open SAM file for reading. Ending program...")
        sys.exit()

    barcodedRead2File = open(output, 'w')

    # write first two lines of sam file (header lines) to new file.
    # start the loop once the pointer is on the actual records
    barcodedRead2File.write(originalSAM.readline())
    barcodedRead2File.write(originalSAM.readline())
    # declare variable to hold one full barcode in one scope higher than the loop
    # so that the barcode from read1 is 'saved' for appending to read2

    # loop through read1 records and apply decoding algorithm
    for count, line in enumerate(originalSAM, start=0):

        # every even line refers to read1. Decode read1 for barcodes. Do not write read1 to new SAM file
        if count % 2 == 0:
            # Make sure that extract_barcode always returns these 3 variables, even if they're empty
            match_obj1, cell_bc, umi = extract_barcode(line, ref_barcode_blocks)

        # every odd line is read2. Read2 will have its sequence appended by the previous read1 barcode
        if count % 2 == 1 and match_obj1:

            barcodedRead2 = append_barcode(line, cell_bc, umi)
            barcodedRead2File.write(barcodedRead2 + '\n')

    originalSAM.close()
    barcodedRead2File.close()
    print("Bad phases: " + str(bad_phase))
    print("Bad blocks: " + str(bad_block))
    print("Low quality blocks: " + str(low_quality))
    print("Bad linkers: " + str(bad_linker))

    return


def get_ref_barcode_blocks(barcode_blocks_file):
    # Function 1 "get_ref_barcode_blocks" reads from a file containing all possible barcode blocks.
    # Returns a list with these reference barcode blocks.

    # read in all 96 possible cell barcode blocks in 'ref_barcode_blocks'
    with open(barcode_blocks_file, 'r') as barcode_blocks_fh:
        ref_barcode_blocks = barcode_blocks_fh.read().splitlines()
        return ref_barcode_blocks


def main():
    # Main function
    # grab SAM filename/path from command line arguments
    parser = argparse.ArgumentParser(description='Process paired reads from SAM file by extracting '
                                                 'barcodes from read 1 and tagging them onto read 2. '
                                                 'Read 1 is removed. Result is an unmapped, unpaired '
                                                 'SAM file.')
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument("-blocks", help='file containing barcode blocks', required=True, metavar='')
    required_group.add_argument("-input", help='.sam input file', required=True, metavar='')
    required_group.add_argument("-output", help='.sam output file', required=True, metavar='')
    args = parser.parse_args()
    # obtain all possible barcode block combinations
    ref_barcode_blocks = get_ref_barcode_blocks(barcode_blocks_file=args.blocks)

    # construct full cell barcodes from every sequence record. Supply the records in SAM format
    read_and_write_sam(all_records=args.input, ref_barcode_blocks=ref_barcode_blocks, output=args.output)

    return


if __name__ == "__main__":
    main()
