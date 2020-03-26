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
  refdict   = readref('ref/ref_s.txt',binder)
  pridict   = readref('ref/pri_s.txt',binder)
  rpridict  = rc_all(pridict)
  refseq    = refdict[binder]
  ref1      = refseq[7:28]
  ref2      = refseq[169:199]
  
  inhandle1 = open(infile1)
  inhandle2 = open(infile2)
  outfile   = open(workpath+'display/'+sample+'.txt','w')
  logfile   = open(workpath+'display_log/'+sample+'.txt','w')
  handle1   = SeqIO.parse(inhandle1,'fastq')
  handle2   = SeqIO.parse(inhandle2,'fastq')
  counts    = [0,0,0,0,0,0]
  #0: not mapped
  #1: constant region mutation
  #2: not paired
  #3: low quality
  #4: indel
  #5: good read
  readcount = 0; start_time = time()
  for record1 in handle1:
    readcount += 1
    if readcount % 100000 == 0:
       print(readcount,time()-start_time,counts)
       #3s for 10000 records.
       #5M lines, 1.25M records, will cost 6 min
    record2  = handle2.next()
    seq1     = str(record1.seq)
    seq2     = str(record2.seq)

    #mapping
    offset1,offset3,dir1 = map_primer(seq1,pridict,binder,rpridict)
    offset2,offset4,dir2 = map_primer(seq2,pridict,binder,rpridict)
    if offset1 == 'NA' or offset2 == 'NA' or offset3 == 'NA' or offset4 == 'NA': 
      counts[0] += 1
      continue
    if dir1 == dir2:
      counts[2] += 1
      continue
    if dir1 == 'F':
      fseq = seq1[offset1:offset3]
      rseq = str(Seq(seq2[offset2:offset4]).reverse_complement())
    if dir1 == 'R':
      fseq = seq2[offset2:offset4]
      rseq = str(Seq(seq1[offset1:offset3]).reverse_complement())
    if len(fseq) != len(rseq) or len(fseq) != len(refseq):
      counts[4] += 1
      continue
   
    #identify backbone mutations
    fc   = fseq[:7]+fseq[28:169]+fseq[199:]
    rc   = rseq[:7]+rseq[28:169]+rseq[199:]
    refc = refseq[:7]+refseq[28:169]+refseq[199:]
    _bg = []
    if fc != rc or refc != fc:
      counts[1] += 1
      for i in range(len(refc)):
        if refc[i] != fc[i] and fc[i] == rc[i]:
          _bg.append(refc[i]+str(i+1)+fc[i])
    
    #extracting loops
    if dir1 == 'F':
      loop1f = record1[offset1+7:offset1+28]
      loop1r = record2[offset2+174:offset2+195]
      loop2f = record1[offset1+169:offset1+199]
      loop2r = record2[offset2+3:offset2+33]    
    else:
      loop1f = record2[offset2+7:offset2+28]
      loop1r = record1[offset1+174:offset1+195]
      loop2f = record2[offset2+169:offset2+199]
      loop2r = record1[offset1+3:offset1+33]    
       
    #compare
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
    outfile.write('\t')
    if len(_bg) == 0:
      outfile.write('WT')
    else:
      outfile.write(':'.join(_bg))
    outfile.write('\n')
    counts[5] += 1
  
  inhandle1.close()
  inhandle2.close()
  outfile.close()
  logfile.write('not mapped: '+str(counts[0])+'\n')
  logfile.write('constant region mutation: '+str(counts[1])+'\n')
  logfile.write('not paired: '+str(counts[2])+'\n')
  logfile.write('low quality: '+str(counts[3])+'\n')
  logfile.write('indel: '+str(counts[4])+'\n')
  logfile.write('good read: '+str(counts[5])+'\n')
  logfile.write('readcount: '+str(readcount)+'\n')
  logfile.close()

def call_mut(loop1f,loop1r,ref1):
    mutlist1 = []
    r1seq  = str(loop1r.seq.reverse_complement())
    f1seq  = str(loop1f.seq)
    r1qual = loop1r.letter_annotations["phred_quality"][::-1]
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

def map_primer(seq,pridict,binder,rpridict):
  pseq    = pridict[binder+'_F'] 
  offset1 = map_read(pseq,seq)
  pseq    = pridict[binder+'_R']
  offset2 = map_read(pseq,seq)
  if offset1 == 'NA' and offset2 == 'NA': return 'NA','NA','NA'
  if offset1 != 'NA': 
    pseq = rpridict[binder+'_R']
    offset3 = map_read(pseq,seq[offset1+197:])
    if offset3 == 'NA': return offset1,'NA','F'
    return offset1, offset3+offset1+197-len(pseq), 'F'
  if offset2 != 'NA': 
    pseq = rpridict[binder+'_F']
    offset4 = map_read(pseq,seq[offset2+197:])
    if offset4 == 'NA': return offset2,'NA','R'
    return offset2, offset4+offset2+197-len(pseq), 'R'
  return 'NA','NA','NA'
  
def map_read(pseq,seq):
  for i in range(0,10):
    qseq = seq[i:i+len(pseq)]
    if len(qseq) < len(pseq):
      return 'NA'
    if hamming(qseq,pseq) < 3:
      return i+len(pseq)
  return 'NA'

def rc_all(pridict):
  newdict = {}
  for pri in pridict:
    newdict[pri] = str(Seq(pridict[pri]).reverse_complement())
  return newdict

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
