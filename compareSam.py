#!/usr/bin/env python3
import argparse  # command line options
import sys
import re
#import profile

###
# compareSam.py reads in two aligned+merged+tagged SAM files, one from Illumina and one from a custom
# pipeline, and locates differences in bar codes. The experiment goal is to determine where the discrepancies
# in gene and transcript counts are coming from when comparing these two pipelines. Consistent differences
# in the bar codes and umis' may suggest key differences in the pipelines' demultiplexing algorithms.


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
        bc_tag = 'XB'
        umi_tag = 'XU'
        # Use regex to find tags since tag order may vary from read to read
        bc_s = re.compile('%s:Z:(\w+)' % bc_tag)  # search for bar code
        bc_f = list(filter(bc_s.match, elements))
        umi_s = re.compile('%s:Z:(\w+)' % umi_tag)  # search for UMI
        umi_f = list(filter(umi_s.match, elements))

        if umi_f and bc_f:  # check for bar code
            s1_bc = bc_f[0][3:24].upper()
            s1_umi = umi_f[0][3:14].upper()
            illumina_dict[elements[0]] = [s1_bc, s1_umi]
        else:
            # when bar code is missing, keep record
            illumina_dict[elements[0]] = []

    illumina.close()

    print('Illumina file successfully closed. Reading custom output...')

    extra_records_custom = 0
    exact_matches = 0

    for line in filter(is_read, custom):
        elements = line.rstrip().split('\t')
        bc_tag = 'XC'
        umi_tag = 'XM'
        gene_tag = 'GE'
        # Use regex to find tags since tag order may vary from read to read
        bc_s = re.compile('%s:Z:(\w+)' % bc_tag) # search for bar code
        bc_f = list(filter(bc_s.match, elements))
        umi_s = re.compile('%s:Z:(\w+)' % umi_tag) # search for UMI
        umi_f = list(filter(umi_s.match, elements))
        gene_s = re.compile('%s:Z:(\w+)' % gene_tag)  # search for gene
        gene_f = list(filter(gene_s.match, elements))

        key = elements[0]

        # check if Illumina has the read id
        if key in illumina_dict:
            # check to see if both umi and bar code are present in custom SAM
            if bc_f and umi_f:
                s2_bc = bc_f[0][3:24].upper()
                s2_umi = umi_f[0][3:14].upper()
                # confirm that illumina record has bar code. Need to compare.
                if len(illumina_dict[key]) == 2:
                    # keep record if bar codes mismatch
                    if illumina_dict[key][0] != s2_bc or illumina_dict[key][1] != s2_umi:
                        # copy over matching records in illumina dictionary
                        compare_dict[key] = illumina_dict[key]
                        # remove illumina dictionary entry (remaining entries are of interest: never appear in custom)
                        del illumina_dict[key]
                        key_c = key + '-C'  # make a new key with similar name (for sorting)
                        compare_dict[key_c] = [s2_bc, s2_umi]

                        if gene_f:  # throw in gene info
                            compare_dict[key_c].append(gene_f[0])
                    # bar codes match up in both pipelines. Custom pipeline worked correctly. Do not keep.
                    else:
                        del illumina_dict[key]
                        exact_matches += 1

                else:
                    # illumina has no bar code on this key. Custom does have a bar code. Keep record.
                    compare_dict[key] = illumina_dict[key]
                    key_c = key + '-C'
                    compare_dict[key_c] = [s2_bc, s2_umi]
                    if gene_f:  # throw in gene info
                        compare_dict[key_c].append(gene_f[0])

            # custom is missing a bar code or umi
            else:
                # check if Illumina does have a bar code + umi. If so, a missing bar code in custom is significant
                if len(illumina_dict[key]) == 2:
                    compare_dict[key] = illumina_dict[key]
                    key_c = key + '-C'
                    compare_dict[key_c] = []
                else:  # both records lack bar codes. Custom pipeline worked correctly. Remove record.
                    del illumina_dict[key]
                    exact_matches += 1

        # Illumina does not have this read.
        else:
            # Custom has the bar codes
            #if umi_f and bc_f:
            # unhelpful for comparison purposes. Nothing to compare to. Likely mapping issue if it occurs.
            extra_records_custom += 1
            # custom has no bar codes and illumina doesnt even have a record
            # unhelpful for comparison purposes. Nothing to compare to. Likely mapping issue if it occurs.
            #extra_records_custom += 1

    extra_records_illumina = len(illumina_dict)
    wrong_matches = int(len(compare_dict) / 2)  # account for how we add custom keys to same dictionary
    custom.close()

    print('Custom SAM file successfully closed. Extracting read1 fastQ sequences...')

    count = 1
    for line in read1:

        if count == 1:
            identifier = line[1:].split()[0]

        elif count == 2:
            if identifier:
                # record key(id)-value(read1 seq) pair. Skip record if ID does not exist.
                if identifier in compare_dict.keys():
                    compare_dict[identifier].append('Read 1: ' + line.rstrip())

        elif count == 4:
            count = 0

        count += 1

    read1.close()

    print('Finished read process and completed dictionaries. Beginning write process...')

    with open(args.output, 'w') as f:
        f.write('Number of records found ONLY in custom pipeline: ' + str(extra_records_custom) + '\n')
        f.write('Number of records found ONLY in Illumina pipeline: ' + str(extra_records_illumina) + '\n')
        f.write('Number of records with exact matches in bar codes: ' + str(exact_matches) + '\n')
        f.write('Number of records with wrong matches in bar codes: ' + str(wrong_matches) + '\n\n')
        sorted_sam_keys = sorted(compare_dict.keys())
        for count,key in enumerate(sorted_sam_keys, 1):
            if len(compare_dict[key]) == 1:  # this accounts for illumina entry only having read1 seq
                f.write(key + '\n')
                print(key)
            else:
                f.write(key + '\t' + str(compare_dict[key][:3]) + '\n')
                print(key + '\t' + str(compare_dict[key][:3]))

            if count % 2 == 0:  # after writing both versions of an entry, get the corresponding read1
                # add original read 1 sequence for comparison
                key = key[:-2]  # remove '-C' since read 1 won't have it
                if len(compare_dict[key]) == 2:  # this shouldn't happen. Illumina missing one of umi/bc/read1.
                    print('ERROR: Illumina read missing a umi, bar code, or read1. Try a different fastq read1.')

                f.write('\n')
    return


if __name__ == "__main__":
    main()
