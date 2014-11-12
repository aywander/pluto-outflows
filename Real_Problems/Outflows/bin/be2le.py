
## This script reads in a binary file of big endianness and outputs the same
## data in small endian binary format. If the -r flag is given, the input file is
## small endian and output is big endian. -f and -d flags stand for float and
## double precision data.
## 
## Usage is:
## be2le.py [-h] [--fin] [--fout] BE-file LE-file

import numpy as np
import argparse


## Argument parser
parser = argparse.ArgumentParser(description=
    'Convert a homogeneous binary file from big endian to little endian.')

parser.add_argument('fname_be', metavar='BE-file', type=str,
                    help='Big endian input file.')
parser.add_argument('fname_le', metavar='LE-file', type=str,
                    help='Little endian output file.')

parser.add_argument('--fin', dest='fin', type=str, action='store', choices=["<f4", ">f4", "<f8", ">f8"],
                    help='Input data type. Choose from "<f4", ">f4", "<f8", ">f8" (where > is big-endian).')
parser.add_argument('--fout', dest='fout', type=str, action='store', choices=["<f4", ">f4", "<f8", ">f8"],
                    help='Output data type. Choose from "<f4", ">f4", "<f8", ">f8" (where > is big-endian).')

args = parser.parse_args()

## Read the file
a = np.fromfile(args.fname_be, args.fin)

## Write the file in opposite endianness
a.astype(args.fout).tofile(args.fname_le)

