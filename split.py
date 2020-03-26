#!/usr/bin/env python
import os
import sys
import glob
import string

path2in = '/u/scratch/t/tianhao/NovaSeq022020/data/'
path2out = '/u/scratch/t/tianhao/NovaSeq022020/split/'

def split(infile1):
  outfile = infile1.rsplit('_001')[0]
  chunksize = 5000000
  fid = 1
  with open(path2in+infile1) as infile:
    f = open(path2out+outfile+'_%03d.fastq' %fid, 'w')
    for i,line in enumerate(infile):
      f.write(line)
      if not (i+1)%chunksize:
        print "file%d" %fid
        f.close()
        fid += 1
        f = open(path2out+outfile+'_%03d.fastq' %fid, 'w')
    f.close()

def main():
  infiles = sorted(glob.glob(path2in+'preclear_Fn22_2loop*.fastq'))
  for infile in infiles:
    infile = infile.rsplit('/')[-1]
    split(infile)

if __name__ == '__main__':
  main()


