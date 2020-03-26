#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
import sys
from distance import hamming
from time import time

def main():
  workpath  = '/u/scratch/t/tianhao/NovaSeq022020/'
  infile1   = workpath+'split/'+sys.argv[1]
  infile2   = infile1.replace('R1','R2')
  sample    = infile1.rsplit('.fastq')[0].rsplit('/')[-1]
  binder    = sample.rsplit('_')[1] 
  refdict   = readref('ref/ref.txt',binder)
  pridict   = readref('ref/pri.txt',binder)
  
  inhandle1 = open(infile1)
  inhandle2 = open(infile2)
  outfile   = open(workpath+'display/'+sample+'.txt','w')
  handle1   = SeqIO.parse(inhandle1,'fastq')
  handle2   = SeqIO.parse(inhandle2,'fastq')
  counts    = [0,0,0,0,0]
  #0: not mapped
  #1: partially mapped
  #2: not paired
  #3: low quality
  #4: indel
  readcount = 0; start_time = time()
  for record1 in handle1:
    readcount += 1
    if readcount % 100000 == 0:
       print(readcount,time()-start_time,counts)
       #30s for 10000 records.
       #5M lines, 1.25M records, will cost 1 hour
    record2  = handle2.next()
    seq1     = str(record1.seq)
    seq2     = str(record2.seq)

    #mapping
    offsets1,dir1 = map_primer(seq1,pridict,binder)
    offsets2,dir2 = map_primer(seq2,pridict,binder)
    if len(offsets1) == 0 or len(offsets2) == 0: 
      counts[0] += 1
      continue
    if len(offsets1) < 4 or len(offsets2) < 4:
      counts[1] += 1
      continue
    if dir1 == dir2:
      counts[2] += 1
      continue
    
    #extracting loops
    if dir1 == 'F':
      loop1f = record1[offsets1[0]+len(pridict[binder+'_Loop1F']):offsets1[0]+offsets1[1]]
      loop2f = record1[offsets1[0]+offsets1[1]+offsets1[2]+len(pridict[binder+'_Loop2F']):sum(offsets1)]
      loop1r = record2[offsets2[0]+offsets2[1]+offsets2[2]+len(pridict[binder+'_Loop1R']):sum(offsets2)]
      loop2r = record2[offsets2[0]+len(pridict[binder+'_Loop2R']):offsets2[0]+offsets2[1]]
    if dir2 == 'F':
      loop1f = record2[offsets2[0]+len(pridict[binder+'_Loop1F']):offsets2[0]+offsets2[1]]
      loop2f = record2[offsets2[0]+offsets2[1]+offsets2[2]+len(pridict[binder+'_Loop2F']):sum(offsets2)]
      loop1r = record1[offsets1[0]+offsets1[1]+offsets1[2]+len(pridict[binder+'_Loop1R']):sum(offsets1)]
      loop2r = record1[offsets1[0]+len(pridict[binder+'_Loop2R']):offsets1[0]+offsets1[1]] 

    #compare
    if len(loop1f) != len(loop1r) or len(loop2f) != len(loop2r):
      counts[2] += 1
      continue
    ref1  = refdict[binder+'_Loop1']
    ref2  = refdict[binder+'_Loop2']
    if len(loop1f) != len(ref1) or len(loop2f) != len(ref2):
      counts[4] += 1
      continue
    mutlist1 = call_mut(loop1f,loop1r,ref1)
    mutlist2 = call_mut(loop2f,loop2r,ref2)
    if 'lowQ' in mutlist1 or 'lowQ' in mutlist2: 
      counts[3] += 1
      continue
    
    #print result
    if len(mutlist1) == 0:
      outfile.write('WT')
    else:
      outfile.write(':'.join(mutlist1))
    outfile.write('\t')
    if len(mutlist2) == 0:
      outfile.write('WT')
    else:
      outfile.write(':'.join(mutlist2))
    outfile.write('\n')
  
  inhandle1.close()
  inhandle2.close()
  outfile.close()

def call_mut(loop1f,loop1r,ref1):
    mutlist1 = []
    r1seq  = str(loop1r.seq.reverse_complement())
    f1seq  = str(loop1f.seq)
    r1qual = loop1r.letter_annotations["phred_quality"]
    f1qual = loop1f.letter_annotations["phred_quality"]
    _low_qual = 0
    for i in range(len(ref1)):
      if f1seq[i] == r1seq[i] and f1seq[i] != ref1[i]:
        mutlist1.append(ref1[i]+str(i+1)+f1seq[i])
      elif f1seq[i] != ref1[i] and f1qual[i] > 30 and r1qual[i] < 30:
        mutlist1.append(ref1[i]+str(i+1)+f1seq[i])
      elif r1seq[i] != ref1[i] and r1qual[i] > 30 and f1qual[i] < 30:
        mutlist1.append(ref1[i]+str(i+1)+r1seq[i])
      elif f1seq[i] != r1seq[i] or f1seq[i] != ref1[i] or r1seq[i] != ref1[i]:
        _low_qual = 1
    if _low_qual == 1: mutlist1.append('lowQ')
    return mutlist1     

def map_primer(seq,pridict,binder):
  pseq    = pridict[binder+'_Loop1F'] 
  offset1 = map_read(pseq,seq)
  pseq    = pridict[binder+'_Loop2R']
  offset4 = map_read(pseq,seq)
  if offset1 == 'NA' and offset4 == 'NA': return [],'NA'
  if offset1 != 'NA': 
    pseq    = str(Seq(pridict[binder+'_Loop1R']).reverse_complement())
    offset2 = map_read(pseq,seq[offset1:])
    if offset2 == 'NA': return [offset1], 'F'
    pseq    = pridict[binder+'_Loop2F']
    offset3 = map_read(pseq,seq[offset1+offset2:])
    if offset3 == 'NA': return [offset1,offset2], 'F'
    pseq    = str(Seq(pridict[binder+'_Loop2R']).reverse_complement())
    offset4 = map_read(pseq,seq[offset1+offset2+offset3:])
    if offset4 == 'NA': return [offset1,offset2,offset3], 'F'
    return [offset1, offset2, offset3, offset4], 'F'
  if offset4 != 'NA': 
    pseq    = str(Seq(pridict[binder+'_Loop2F']).reverse_complement())
    offset3 = map_read(pseq,seq[offset4:])
    if offset3 == 'NA': return [offset4], 'R'
    pseq    = pridict[binder+'_Loop1R'] 
    offset2 = map_read(pseq,seq[offset4+offset3:])
    if offset2 == 'NA': return [offset4,offset3], 'R'
    pseq    = str(Seq(pridict[binder+'_Loop1F']).reverse_complement())
    offset1 = map_read(pseq,seq[offset4+offset3+offset2:])
    if offset1 == 'NA': return [offset4,offset3,offset2], 'R'
    return [offset4,offset3,offset2,offset1], 'R'
  return [],'NA'
  
def map_read(pseq,seq):
  for i in range(0,len(seq)-len(pseq)):
    qseq = seq[i:i+len(pseq)]
    if hamming(qseq,pseq) < 3:
      return i
  return 'NA'

def readref(infile,binder):
  inhandle = open(infile)
  refdict  = {}
  for record in SeqIO.parse(inhandle,'fasta'):
    if binder not in str(record.id): continue
    refdict[str(record.id)] = str(record.seq).upper()
  inhandle.close()
  return refdict


if __name__ == '__main__':
  main()
