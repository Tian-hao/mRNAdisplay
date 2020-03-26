#!/usr/bin/env python
import os
import sys
import glob
import string
from Bio.Seq import Seq
from Bio import SeqIO
from time import time


def parse_mut(mut,ref):
  if mut == 'WT':
    return ref
  mut = mut.rsplit(':')
  newseq = ref
  for m in mut:
    newseq = newseq[:int(m[1:-1])-1]+m[-1]+newseq[int(m[1:-1]):]
  return newseq

def callmut(aa,ref):
  refaa = str(Seq(ref).translate())
  mutlist = []
  for i in range(len(aa)):
    if aa[i] != refaa[i]:
      mutlist.append(refaa[i]+str(i+1)+aa[i])
  return ':'.join(mutlist)

def makefull(nucs,refseq):
  newseq = refseq
  nucs = nucs.rsplit('_')
  if nucs[0] != 'WT':
    muts = nucs[0].rsplit(':')
    for mut in muts:
      newseq = newseq[:7+int(mut[1:-1])]+mut[-1]+newseq[7+int(mut[1:-1])+1:]
  if nucs[1] != 'WT':
    muts = nucs[1].rsplit(':')
    for mut in muts:
      newseq = newseq[:169+int(mut[1:-1])]+mut[-1]+newseq[169+int(mut[1:-1])+1:]
  if nucs[2] != 'WT':
    muts = nucs[2].rsplit(':')
    for mut in muts:
      pos = int(mut[1:-1])
      if pos < 7:
        newseq = newseq[:pos]+mut[-1]+newseq[pos+1:]
      elif pos < 148:
        newseq = newseq[:pos+21]+mut[-1]+newseq[pos+22:]
      else:
        newseq = newseq[:pos+51]+mut[-1]+newseq[pos+52:]
  return newseq

def main():
  workpath = '/u/scratch/t/tianhao/NovaSeq022020/display/'
  #binders  = ['Fn17','Fn22']
  loops    = ['1loop','2loop']
  binders   = [sys.argv[1]]
  #loops     = ['2loop']
  for binder in binders:
    for loop in loops:
      main_mapper(binder,loop,workpath)

def main_mapper(binder,loop,workpath):
  print(binder,loop)
  libs     = ['Input','N-lysate','NP-lysate','N-purified','preclear']
  genodict = {} #genotype: lib_name: count
  aadict   = {} #aaseq: lib_name: count
  depdict  = {} #lib_name: count 
  outfile  = open('result/display_readcount_'+binder+'_'+loop+'.txt','w') #a small result file grouped by aa sequence
  genofile = open(workpath+'/result_nuc_display_readcount_'+binder+'_'+loop+'.txt','w') #full result with each nucleotide info
  depfile  = open('result/display_depth_'+binder+'_'+loop+'.txt','w')
  logresult= open('result/display_log_'+binder+'_'+loop+'.txt','w')
  for lib in libs:
    depdict[lib] = 0
      
  reffile = open('ref/ref.txt')
  refdict = {}
  for record in SeqIO.parse(reffile,'fasta'):
    refdict[str(record.id)] = str(record.seq).upper()
  reffile.close()  
  refseq = refdict[binder]

  infiles = sorted(glob.glob(workpath+'*.txt'))
  start_time = time()
  for infile in infiles:
    if binder != infile.rsplit('_')[1]: continue
    if loop   != infile.rsplit('_')[2]: continue
    lib    = infile.rsplit('_')[0].rsplit('/')[-1]
    print(infile, 'using time', time()-start_time)
    handle = open(infile)
    for line in handle:
      if line[-1] != '\n': break
      line = line.rstrip().rsplit('\t')
      depdict[lib] += 1
      nucs = line[0]+'_'+line[1]+'_'+line[2]
      if nucs not in genodict: 
        genodict[nucs] = {}
        for libi in libs:
          genodict[nucs][libi] = 0
      genodict[nucs][lib] += 1
    handle.close()
  
  #writefile
  depfile.write('Library\tTotal depth\n')
  for lib in libs:
    depfile.write(lib+'\t'+str(depdict[lib])+'\n')
  depfile.close()
   
  genofile.write('Loop1Seq\tLoop2Seq\tLoop1aa\tLoop2aa\tLoop1Mut\tLoop2Mut\tfullSeq\tbackboneMut\tbackboneAA')
  for lib in libs:
    genofile.write('\t'+lib)
  genofile.write('\n')
  print('total count: '+str(len(genodict)))
  nuccount = 0
  start_time = time()
  for nucs in genodict:
    nuccount += 1
    if nuccount % 1000000 == 0: 
      print('processing mutation count '+str(nuccount)+' using time '+str(time()-start_time))
    nucss = nucs.rsplit('_')
    nuc1 = parse_mut(nucss[0],refdict[binder+'_Loop1'])
    nuc2 = parse_mut(nucss[1],refdict[binder+'_Loop2'])
    bgmut = nucss[2]
    fullseq = makefull(nucs,refseq)
    bgaa = callmut(str(Seq('A'+fullseq+'C').translate()),'A'+refseq+'C')
    aa1 = str(Seq(nuc1).translate())
    aa2 = str(Seq(nuc2).translate())
    mut1 = callmut(aa1,refdict[binder+'_Loop1'])
    mut2 = callmut(aa2,refdict[binder+'_Loop2'])
    genofile.write(nuc1+'\t'+nuc2+'\t'+aa1+'\t'+aa2+'\t'+mut1+'\t'+mut2+'\t'+fullseq+'\t'+bgmut+'\t'+bgaa)
    for lib in libs:
      genofile.write('\t'+str(genodict[nucs][lib]))
    genofile.write('\n')

    aas = aa1+'_'+aa2
    if aas not in aadict:
      aadict[aas] = {}
      for lib in libs:
        aadict[aas][lib] = 0
    for lib in libs:
      aadict[aas][lib] += genodict[nucs][lib]
  genofile.close()

  outfile.write('Loop1\tLoop2')
  for lib in libs:
    outfile.write('\t'+lib)
  outfile.write('\n')
  for muts in aadict:
    mutss = muts.rsplit('_')
    outfile.write(mutss[0]+'\t'+mutss[1])
    for lib in libs:
      outfile.write('\t'+str(aadict[muts][lib]))
    outfile.write('\n')
  outfile.close()

  logfiles = sorted(glob.glob('/u/scratch/t/tianhao/NovaSeq022020/display_log/*.txt'))
  logdict = {}
  for logfile in logfiles:
    if binder != logfile.rsplit('_')[2]: continue
    if loop   != logfile.rsplit('_')[3]: continue
    lib    = logfile.rsplit('_')[1].rsplit('/')[-1]
    handle = open(logfile)
    for line in handle:
      line = line.rstrip().rsplit(': ')
      record = line[0]
      count = float(line[1])
      if record not in logdict: logdict[record] = {}
      if lib not in logdict[record]: logdict[record][lib] = 0
      logdict[record][lib] += count
    handle.close()
  
  logresult.write('lib')
  for record in logdict:
    logresult.write('\t'+record)
  logresult.write('\n')
  for lib in libs:
    logresult.write(lib)
    for record in logdict:
      logresult.write('\t'+str(logdict[record][lib]))
    logresult.write('\n')
  logresult.close()

      
if __name__ == '__main__':
  main()
