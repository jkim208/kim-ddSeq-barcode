#!/usr/bin/env python3
import argparse  # command line options
import regex  # regular expressions with edit distance functions
import pyximport  # install Cython to access pyximport
pyximport.install(build_in_temp=False)
from editDistance import edit_distance
import sys
import distance
#import profile

### LOW QUALITY CHECK ON UMI
### Lets assume that linkers have the highest chance of mutation error
### Use anchors/umi to track
### Compared to parseBarcodes5, this script most uses matching. Minimal indexing for linker checks

# parseBarcodes.py was written for the BioRad ddSeq procedure (scRNA).
# This program applies a decoding algorithm recommended by Illumina to extract, parse and filter
# barcodes in read 1 of a paired set. The barcode information is then tagged onto the SAM file record
# of read 2.


def append_barcode(line, cell_bc, umi):
    # Function 6 "append_barcode" takes the barcode and adds it to the record as a separate tag
    line = "%s\t%s\t%s" % (line.rstrip(), 'XC:Z:' + cell_bc, 'XM:Z:' + umi)
    return line


def check_bc_quality(q_seq, bc_index, low_q_count):
    # Function 5 'check_bc_quality' accepts a quality sequence string, a barcode index, and current low
    # quality count. Via the barcode index, the function checks if the corresponding quality string
    # indicates any low quality bases. If any are detected, loq_q_count is incremented by 1.
    # Function 4 is used in function 3 (demultiplexing) and returns the new low_q_count.

    for element in q_seq[bc_index:bc_index + 6]:
        # Break out of loop immediately if low_q_count threshold has been passed
        if low_q_count > 2:
            break

        # Check the ASCII code if the base quality is low (q-score 10 = ASCII score 43)
        if ord(element) < 43:
            low_q_count += 1

    return low_q_count


def correct_bc_blocks(ref_barcode_blocks, barcode_block):
    # Function 4 "correct_bc_blocks" takes each barcode_block from a sequence, checks the hamming distances
    # between the variable block and all 96 possible blocks, and makes corrections depending on the smallest
    # ED(edit distance) found. Barcode block is unchanged if smallest ED is 2.

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
            barcode_block = lh_reference_block
            return barcode_block


    # find the blocks with the smallest hamming distances. If that ED is >1, skip the read
    barcode_block = None

    return barcode_block


def extract_barcode(line, ref_barcode_blocks, bad_phase, bad_block, low_quality, bad_linker):
    # Function 3: "extract_barcode" uses regex to extract barcode blocks and return complete barcodes

    # split read 1 to extract relevant parameters
    read1 = line.rstrip().split('\t')
    # match to where the sequence should be
    match_obj1 = read1[9]
    phase_blocks = ['', 'A','CT','GCA','TGCG','ATCGA']
    shift = 0
    linker_ed = 0

    if match_obj1:
        # match blocks accordingly with linkers, leaving room for 1 edit distance
        # only keep reads where ACG and GACT anchors are not mutated
        # main unaccounted case is if insertions occur before/after barcode and before ACG
        #ACGAGTACGGTGAC
        # Consider 1 deletion or insertion
        # Assume no insertion/deletion at UMI and barcodes (if it happens, better to drop read)
        capture = regex.search(r"^([AGCT]{0,5})([AGCT]{6})(TAG[AGCT]{11,13})([AGCT]{6})(TAC[AGCT]{11,13})"
                               r"([AGCT]{6})ACG([AGCT]{8})GAC[T]{0,}$", match_obj1)
        if not capture:
            print('problems with capture. Drop read.')
            return None, None, None, bad_phase, bad_block, low_quality, bad_linker

        linker1 = capture.group(3)
        linker2 = capture.group(5)

        if linker1 and linker2:
            pb = capture.group(1)

            if pb in phase_blocks:

                bc1 = capture.group(2)

                if distance.levenshtein(linker1, 'TAGCCATCGCATTGC') > 0:
                    # that's one edit distance
                    linker_ed += 1

                bc2 = capture.group(4)

                if distance.levenshtein(linker2, 'TACCTCTGAGCTGAA') == 1 and linker_ed > 0:
                    # problem. Too many EDs. Drop the read.
                    # bad linkers
                    bad_linker += 1
                    return None, None, None, bad_phase, bad_block, low_quality, bad_linker

                bc3 = capture.group(6)
                umi = capture.group(7)

                bc1 = correct_bc_blocks(ref_barcode_blocks, bc1)
                bc2 = correct_bc_blocks(ref_barcode_blocks, bc2)
                bc3 = correct_bc_blocks(ref_barcode_blocks, bc3)

                if bc1 and bc2 and bc3:
                    # for bc1, beginning is at index 0.
                    # for bc2, beginning is after at least bc1 (6b) and linker1 (15b)
                    # for bc3, beginning is after at least bc2 (12b) and linker2 (30b)
                    # 68 is the length of the sequence. Hard code it for speed. Otherwise, remove.
                    bc1_index = match_obj1.find(bc1, 0, 11)
                    bc2_index = match_obj1.find(bc2, 20, 33)
                    bc3_index = match_obj1.find(bc3, 41, 54)
                    umi_index = match_obj1.find(umi, 50, 65)

                    # count number of low quality bases.
                    low_quality_count = 0
                    # break the loop and remove the read combo if count is > 2
                    low_quality_count = check_bc_quality(read1[10], bc1_index, low_quality_count)
                    low_quality_count = check_bc_quality(read1[10], bc2_index, low_quality_count)
                    low_quality_count = check_bc_quality(read1[10], bc3_index, low_quality_count)
                    low_quality_count = check_bc_quality(read1[10], umi_index, low_quality_count)

                    # Up to 2 low quality barcode bases allowed
                    if low_quality_count < 2:
                        cell_bc = bc1 + bc2 + bc3
                        return match_obj1, cell_bc, umi, bad_phase, bad_block, low_quality, bad_linker

                    else:
                        low_quality += 1
                        pass
                else:
                    bad_block += 1
                    pass

            else:
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

    return match_obj1, cell_bc, umi, bad_phase, bad_block, low_quality, bad_linker


def read_and_write_sam(all_records, ref_barcode_blocks, output):
    # Function 2 "read_and_write_sam" accounts for edit distance while extracting barcodes
    # Includes the correct_bc_blocks function in order to return full barcode
    bad_phase = 0
    bad_block = 0
    low_quality = 0
    bad_linker = 0

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

            match_obj1, cell_bc, umi, bad_phase, bad_block, low_quality, bad_linker = extract_barcode(
                line, ref_barcode_blocks, bad_phase, bad_block, low_quality, bad_linker)

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
    # start = timeit.default_timer()
    # obtain all possible barcode block combinations
    ref_barcode_blocks = get_ref_barcode_blocks(barcode_blocks_file=args.blocks)

    # construct full cell barcodes from every sequence record. Supply the records in SAM format
    read_and_write_sam(all_records=args.input, ref_barcode_blocks=ref_barcode_blocks, output=args.output)

    # stop = timeit.default_timer()
    # print stop - start
    return


if __name__ == "__main__":
    #profile.run("main()")
    main()
