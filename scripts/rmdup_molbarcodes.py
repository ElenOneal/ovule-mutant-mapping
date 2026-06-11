#!/usr/bin/python
#
# Split raw fastqs by i7 index
#  identify/remove PCR duplicates using i5 molecular barcode
#
#
# Script written by Thom Nelson

import argparse
import gzip
import sys
import re

# parse command line options
parser = argparse.ArgumentParser(description='Use molecular barcodes on i5 index to remove PCR duplicates. Fastq files should be gzipped and of the form:\n  prefix.[1,2,i5,i7].suffix')
parser.add_argument('-p','--prefix', required=True, help='prefix for gzipped fastq, up until, e.g., .1.')
parser.add_argument('-s','--suffix', required=True, help='fastq suffix, e.g.: fq.gz or fastq.gz')
args=parser.parse_args()

prefix = args.prefix
suffix = args.suffix
fwdin  = gzip.open(prefix+".1."+suffix,'r')
revin  = gzip.open(prefix+".2."+suffix,'r')
i5in   = gzip.open(prefix+".i5."+suffix,'r')
i7in   = gzip.open(prefix+".i7."+suffix,'r')


# CREATE OUTPUT FILES FOR EACH INPUT FILE
fwdout = prefix+'.rmdup.1.'+suffix
revout = prefix+'.rmdup.2.'+suffix
i5out  = prefix+'.rmdup.i5.'+suffix
i7out  = prefix+'.rmdup.i7.'+suffix
fwdout = gzip.open(fwdout,'w')
revout = gzip.open(revout,'w')
i5out  = gzip.open(i5out,'w')
i7out  = gzip.open(i7out,'w')

# MOVE THROUGH FASTQ FILES

total = 0 # total reads
dupes  = 0 # reads representing PCR duplicates
passed = 0 # reads with valid barcodes that aren't PCR dupes

# INITIALIZE DICTIONARY OF LISTS TO STORE
#  POTENTIAL PCR DUPES. KEYS WILL BE MOLECULAR BARCODES,
#  VALUES WILL BE LISTS OF CONCATENATED SEQUENCES:
#  FIRST 10 BASES FROM FWD READ, FIRST 10 FROM REV.
uniqseqs = {}

# OPEN A LOG FILE
logout = open(prefix+".rmdup.log",'w')
sys.stderr.write("\n")

while True:
    if total % 1000 == 0:
        sys.stderr.write("%s -- clusters examined: %s; PCR dupes ID'd: %s\r"%(prefix,total,dupes))
        logout.write("%s -- clusters examined: %s; PCR dupes ID'd: %s\n"%(prefix,total,dupes))
    i7header = i7in.readline()
    # check if we've arrived at the end of the file
    if not i7header:
        sys.stderr.write("%s -- clusters examined: %s; PCR dupes ID'd: %s\r"%(prefix,total,dupes))
        logout.write("%s -- clusters examined: %s; PCR dupes ID'd: %s\n"%(prefix,total,dupes))        
        break
    i7seq   = i7in.readline().strip()
    i7plus  = i7in.readline()
    i7phred = i7in.readline()
    total   += 1
    goahead = False
    # use molecular barcode & first 10 bp of both reads to check if PCR duplicate
    i5header = i5in.readline()
    i5seq    = i5in.readline()
    i5plus   = i5in.readline()
    i5phred  = i5in.readline()
    r1header = fwdin.readline()
    r1seq    = fwdin.readline()
    r1plus   = fwdin.readline()
    r1phred  = fwdin.readline()
    r2header = revin.readline()
    r2seq    = revin.readline()
    r2plus   = revin.readline()
    r2phred  = revin.readline()
    r1test   = list(r1seq.strip()) ; r1test = ''.join(r1test[0:10])
    r2test   = list(r2seq.strip()) ; r2test = ''.join(r2test[0:10])
    readtest = r1test+r2test
    if i5seq not in uniqseqs:
        fwdout.write("%s%s%s%s"%(r1header,r1seq,r1plus,r1phred))
        revout.write("%s%s%s%s"%(r2header,r2seq,r2plus,r2phred))
        i5out.write("%s%s%s%s"%(i5header,i5seq,i5plus,i5phred))
        i7out.write("%s%s%s%s"%(i7header,i7seq,i7plus,i7phred))
        uniqseqs[i5seq] = [readtest]
    else:
        if readtest not in uniqseqs[i5seq]:
            fwdout.write("%s%s%s%s"%(r1header,r1seq,r1plus,r1phred))
            revout.write("%s%s%s%s"%(r2header,r2seq,r2plus,r2phred))
            i5out.write("%s%s%s%s"%(i5header,i5seq,i5plus,i5phred))
            i7out.write("%s%s%s%s"%(i7header,i7seq,i7plus,i7phred))
            uniqseqs[i5seq].append(readtest)
        else:
            dupes += 1
fwdin.close()
fwdout.close()
revin.close()
revout.close()
i5in.close()
i5out.close()
i7in.close()
i7out.close()
logout.write("\n\nClusters examined for sample %s\n"%(prefix))
logout.write("Total clusters examined: %s\n"%(total))
logout.write("PCR duplicates ID'd:     %s\n"%(dupes))
sys.stderr.write("\n\nClusters examined for sample %s\n"%(prefix))
sys.stderr.write("Total clusters examined: %s\n"%(total))
sys.stderr.write("PCR duplicates ID'd:     %s\n"%(dupes))
logout.close()