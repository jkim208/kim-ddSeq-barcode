#!/usr/bin/env python3
import argparse  # command line options
import sys


def is_read(s):
   return not s.startswith('@')


def main():
    # Main function
    # grab SAM filename/path from command line arguments
    parser = argparse.ArgumentParser(description='')
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument("-illumina", help='Illumina SAM file', required=True, metavar='')
    required_group.add_argument("-custom", help='Custom Pipeline SAM file', required=True, metavar='')
    required_group.add_argument("-output", help='.sam output file', required=True, metavar='')
    args = parser.parse_args()

    sam_dict = {}

    try:
        illumina = open(args.illumina, 'r')
    except IOError:
        print("Could not open SAM file for reading. Ending program...")
        sys.exit()

    for line in filter(is_read, illumina):
        elements = line.rstrip().split('\t')
        #sam_dict[elements.pop(0)] = elements
        sam_dict[elements[0]] = elements
    illumina.close()

    try:
        custom = open(args.custom, 'r')
    except IOError:
        print("Could not open SAM file for reading. Ending program...")
        sys.exit()

    for line in filter(is_read, custom):
        elements = line.rstrip().split('\t')
        #key = elements.pop(0)
        key = elements[0]
        if key in sam_dict:
            key = key + ' custom'
            sam_dict[key] = elements

    print('Beginning write process...')

    with open(args.output, 'w') as f:
        for stuff in sorted(sam_dict.keys()):
            f.write(str(sam_dict[stuff]) + '\n')

    return


if __name__ == "__main__":
    main()




