import sys
sys.path.append('camb4py-master/build/lib.linux-x86_64-2.7/camb4py/')
import camb4py
import argparse

# loading camb executable
camb = camb4py.load('CAMB/camb', defaults='CAMB/params.ini')

# reading in parameters.
descr = 'program to write CAMB Pks'
parser = argparse.ArgumentParser(description=descr)
parser.add_argument('--version', action='version', version='%(prog)s v1')
parser.add_argument('--H_0', metavar = 'H_0', type = float, default = 70.0)
parser.add_argument('--ombh2', metavar = 'ombh2', type = float, default = 0.0226)
parser.add_argument('--omch2', metavar = 'omch2', type = float, default = 0.112)
parser.add_argument('--omnuh2', metavar = 'omnuh2', type = float, default = 0.00064)
parser.add_argument('--omk', metavar = 'omk', type = float, default = 0.0)
parser.add_argument('--z', metavar = 'z', type = float, default = 7.0)

args = parser.parse_args()

camb_result_dict = camb(**{'ombh2':args.ombh2, 'omch2':args.omch2, 'omnuh2':args.omnuh2, 'omk':args.omk, 'hubble':args.H_0, 'transfer_redshift(1)':args.z})

table = camb_result_dict["transfer_matterpower"]
nvals = len(table)
f = open("Pks.dat", "w")
for n in range(0, nvals):
    f.write(str(table[n][0])+" "+str(table[n][1])+"\n")
