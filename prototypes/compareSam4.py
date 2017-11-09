#!/usr/bin/env python3
import argparse  # command line options
import sys
import re
#import profile

###
# Updates
# Start throwing out entire Illumina records that lack any barcodes.

###
# compareSam.py reads in two aligned+merged+tagged SAM files, one from Illumina and one from a custom
# pipeline, and locates differences in bar codes. The experiment goal is to determine where the discrepancies
# in gene and transcript counts are coming from when comparing these two pipelines. Consistent differences
# in the bar codes and umis' may suggest key differences in the pipelines' demultiplexing algorithms.
extra_records_custom = 0
exact_matches = 0

def is_read(s):
    return not s.startswith('@')


def create_illumina_dictionary(illumina):
    illumina_dict = {}
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
        else:  # throw out all barcode-lacking records
            pass

    illumina.close()
    return illumina_dict


def create_compare_dictionary(illumina_dict, custom):
    compare_dict = {}
    custom_dict = {}
    global extra_records_custom
    global exact_matches

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

        # check if Illumina has the read id. If it does, it must have barcodes for that id
        if key in illumina_dict:
            # check to see if both umi and bar code are present in custom SAM. They should be.
            if bc_f and umi_f:
                s2_bc = bc_f[0][3:24].upper()
                s2_umi = umi_f[0][3:14].upper()
                # We know both custom and Illumina have barcodes for this record. Keep record if bar codes mismatch
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
            # custom is missing a bar code or umi #### THIS SHOULD NEVER HAPPEN
            else:
                print('Custom is missing barcodes/umi at ' + key)

        # Illumina dropped this read or does not have barcodes for this read. Custom however does.
        else:
            custom_dict[key] = elements
            # Custom has the bar codes
            #if umi_f and bc_f:
            # unhelpful for comparison purposes. Nothing to compare to. Likely mapping issue if it occurs.
            extra_records_custom += 1
            # custom has no bar codes and illumina doesnt even have a record

    return compare_dict, custom_dict


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

    print('Opening Illumina file for reading...')

    # Open files early to detect errors immediately
    illumina = read_file(args.illumina)
    custom = read_file(args.custom)
    read1 = read_file(args.read1)

    illumina_dict = create_illumina_dictionary(illumina)

    print('Illumina file successfully closed. Reading custom output...')

    compare_dict, custom_dict = create_compare_dictionary(illumina_dict, custom)
    extra_records_illumina = len(illumina_dict)
    wrong_matches = int(len(compare_dict) / 2)  # account for how we add custom keys to same dictionary
    custom.close()

    print('Custom SAM file successfully closed. Extracting read1 fastQ sequences...')

    illumina_nonunique_id_count = 0
    count = 1
    for line in read1:

        if count == 1:
            identifier = line[1:].split()[0]  # clean first line

        elif count == 2:
            if identifier:
                # record key(id)-value(read1 seq) pair. Skip record if ID does not exist.
                # add read1 sequence info to dictionaries of interest
                if identifier in compare_dict.keys():
                    compare_dict[identifier].append('Read 1: ' + line.rstrip())
                if identifier in custom_dict.keys():
                    custom_dict[identifier].append('Read 1: ' + line.rstrip())
                if identifier in illumina_dict.keys():
                    illumina_dict[identifier].append('Read 1: ' + line.rstrip())
                    illumina_nonunique_id_count += 1
                    # note: entries missing a read1 are oddities. The records don't show up in the fastq

        elif count == 4:
            count = 0  # reset counter

        count += 1

    read1.close()

    print('Finished read process and completed dictionaries. Beginning write process...')

    with open(args.output, 'w') as f:
        f.write('Number of records found ONLY in custom pipeline: ' + str(extra_records_custom) + '\n')
        f.write('Number of records found ONLY in Illumina pipeline: ' + str(extra_records_illumina) + '\n')
        f.write('Number of records with exact matches in bar codes: ' + str(exact_matches) + '\n')
        f.write('Number of records with wrong matches in bar codes: ' + str(wrong_matches) + '\n')
        f.write('Number of illumina records with IDs matching fastq: ' + str(illumina_nonunique_id_count) + '\n\n')
        sorted_sam_keys = sorted(compare_dict.keys())
        for count, key in enumerate(sorted_sam_keys, 2):
            if count % 2 == 0:
                f.write(key + '\t' + str(compare_dict[key][:3]) + '\n')
            else:
                f.write(key + '\t' + str(compare_dict[key][:3]) + '\n\n')
        f.write("\n\nEnd of one to one comparisons\n\n")
        f.write("Illumina-unique reads------------------------------------------------\n\n")
        for key in illumina_dict:
            f.write(key + '\t' + str(illumina_dict[key]) + '\n')
        f.write("\nCustom-unique reads--------------------------------------------------\n\n")
        for key in custom_dict:
            f.write(key + '\t' + str(custom_dict[key][1:]) + '\n')

    return


if __name__ == "__main__":
    main()
