#!/usr/bin/env python3
import argparse  # command line options
import sys
import re

###
# compareSam.py reads in two aligned+merged+tagged SAM files, one from Illumina and one from a custom
# pipeline, and locates differences in bar codes. The experiment goal is to determine where the discrepancies
# in gene and transcript counts are coming from when comparing these two pipelines. Consistent differences
# in the bar codes and umis' may suggest key differences in the pipelines' demultiplexing algorithms.

# Note: these indices assume list starts at 0 (10th element has index of 9)
s1_id_i = 0  # indices for illumina pipeline (id, seq, bar code, umi, gene)
s1_seq_i = 9

s2_id_i = 0  # indices for custom pipeline (id, seq, bar code, umi, gene)
s2_seq_i = 9


def is_read(s):
    return not s.startswith('@')


def read_file(file):
    try:
        filehandle = open(file, 'r')
    except IOError:
        print("Could not open file for reading. Ending program...")
        sys.exit()
    return filehandle


def main():
    # Main function
    # grab SAM filename/path from command line arguments
    parser = argparse.ArgumentParser(description='')
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument("-illumina", help='Illumina SAM file', required=True, metavar='')
    required_group.add_argument("-custom", help='Custom Pipeline SAM file', required=True, metavar='')
    required_group.add_argument("-read1", help ='Original read1 fastq file', required=True, metavar='')
    required_group.add_argument("-output", help='.sam output file', required=True, metavar='')
    args = parser.parse_args()

    illumina_dict = {}
    compare_dict = {}

    print('Opening Illumina file for reading...')

    # Open files early to detect errors immediately
    illumina = read_file(args.illumina)
    custom = read_file(args.custom)
    read1 = read_file(args.read1)

    for line in filter(is_read, illumina):
        elements = line.rstrip().split('\t')

        # Use regex to find tags since tag order may vary from read to read
        bc_r = re.compile("XB:Z:(\w+)")  # search for bar code
        bc_f = list(filter(bc_r.match, elements))
        umi_r = re.compile("XU:Z:(\w+)")  # search for umi
        umi_f = list(filter(umi_r.match, elements))

        if umi_f and bc_f:  # check for bar code
            s1_bc = bc_f[0][3:24].upper()
            s1_umi = umi_f[0][3:14].upper()
            illumina_dict[elements[0]] = s1_bc, s1_umi
        else:
            # when bar code is missing, keep record
            illumina_dict[elements[0]] = ' '

    illumina.close()

    print('Illumina file successfully closed. Reading custom output...')

    extra_records_custom = 0
    exact_matches = 0

    for line in filter(is_read, custom):
        elements = line.rstrip().split('\t')

        # Use regex to find tags since tag order may vary from read to read
        bc_r = re.compile("XC:Z:(\w+)")  # search for bar code
        bc_f = list(filter(bc_r.match, elements))
        umi_r = re.compile("XM:Z:(\w+)")  # search for umi
        umi_f = list(filter(umi_r.match, elements))
        gene_r = re.compile("GE:Z:(\w+)")  # search for gene
        gene_f = list(filter(gene_r.match, elements))

        key_c = elements[0]

        if umi_f and bc_f and key_c in illumina_dict:  # check for bar code and if record is in Illumina
            s2_bc = bc_f[0][3:24].upper()
            s2_umi = umi_f[0][3:14].upper()

            if len(illumina_dict[key_c]) == 3:  # confirm that illumina record has bar code. Need to compare.
                if illumina_dict[key_c][1] != s2_bc or illumina_dict[key_c][2] != s2_umi:  # keep record if bar codes mismatch
                    compare_dict[key_c] = illumina_dict[key_c]  # copy over matching records in illumina dictionary
                    # pull from original dictionary (remaining entries are of interest: never appear in custom)
                    del illumina_dict[key_c]
                    key_c = key_c + '-C'  # make a new key with similar name (for sorting)
                    compare_dict[key_c] = s2_bc, s2_umi
                    if gene_f:
                        compare_dict[key_c].append(gene_f[0])
                else:  # bar codes match up in both pipelines. Custom pipeline worked correctly.
                    del illumina_dict[key_c]
                    exact_matches += 1

            else:
                # illumina has no bar code on this key. Custom does have a bar code. Keep record.
                compare_dict[key_c] = illumina_dict[key_c]
                key_c = key_c + '-C'
                compare_dict[key_c] = s2_bc, s2_umi

        elif key_c in illumina_dict and not umi_f:  # custom is missing a bar code and umi
            # check if Illumina does have a bar code. If so, a missing bar code in custom is significant
            if len(illumina_dict[key_c]) == 3:
                compare_dict[key_c] = illumina_dict[key_c]
                key_c = key_c + '-C'
                compare_dict[key_c] = ' '
            else:  # both records lack bar codes. Custom pipeline worked correctly. Remove record.
                del illumina_dict[key_c]
                exact_matches += 1

        elif umi_f and bc_f and key_c not in illumina_dict:  # there is a bar code, but no record in Illumina
            # unhelpful for comparison purposes. Nothing to compare to. Likely mapping issue if it occurs.
            extra_records_custom += 1

        else:  # custom has no bar codes and illumina doesnt even have a record
            # unhelpful for comparison purposes. Nothing to compare to. Likely mapping issue if it occurs.
            extra_records_custom += 1

    extra_records_illumina = len(illumina_dict)
    wrong_matches = int(len(compare_dict) / 2)  # account for how we add custom keys to same dictionary
    custom.close()

    print('Custom SAM file successfully closed. Extracting read1 fastQ sequences...')

    read1_dict = {}
    count = 1
    for line in read1:
        if count == 1:
            identifier = line[1:].split()[0]
        elif count == 2:
            if identifier:
                read1_dict[identifier] = line  # record key(id)-value(read1 seq) pair. Skip record if ID does not exist.
        elif count == 4:
            count = 0
        count += 1

    print('Finished read process and completed dictionaries. Beginning write process...')

    with open(args.output, 'w') as f:
        f.write('Number of records found ONLY in custom pipeline: ' + str(extra_records_custom) + '\n')
        f.write('Number of records found ONLY in Illumina pipeline: ' + str(extra_records_illumina) + '\n')
        f.write('Number of records with exact matches in bar codes: ' + str(exact_matches) + '\n')
        f.write('Number of records with wrong matches in bar codes: ' + str(wrong_matches) + '\n\n')
        sorted_sam_keys = sorted(compare_dict.keys())
        for count,key in enumerate(sorted_sam_keys, 1):
            f.write(key + '\t' + str(compare_dict[key]) + '\n')
            # print(key + '\t' + str(illumina_dict[key]))
            if count % 2 == 0:
                # add original read 1 sequence for comparison
                key = key[:-2]  # remove '-C' since read 1 won't have it
                if key in read1_dict:
                    f.write('Read 1 sequence was: ' + read1_dict[key])
                else:  # this shouldn't happen
                    print('ERROR: Identifier from input not present in original read1 fastq file. Try another file.')
                # print('\n')
                f.write('\n')

    return


if __name__ == "__main__":
    main()




