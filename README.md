# mRNAdisplay

codes to map mRNA display data and find hits that strongly bind to the bait. 

also a tensor flow implementation to cross validate screen data. 

## 1. System requirements

All software dependencies and operating systems (including version numbers).  
python > 3.6.  
R > 3.2.  
Matplotlib 3.2.1    
Tensorflow 2.0    
Biopython 1.76.  
Mapper1 and Mapper2 was run on UCLA's cluster server, Hoffman2. It is powered by Sun Grid Engine.  
http://www.hoffman2.idre.ucla.edu.  
analysis.ipynb was run on a CentOS 6.7 desktop with Anaconda.   
 
Versions the software has been tested on:  
I only tested the sottware on Hoffman2 and CentOS 6.7.   

Any required non-standard hardware:  
No.  

## 2. Installation guide  
Instructions:   
1. On a Sun Grid Engine, load biopython first.   
Use mywrapper.sh to run Mapper1.py in batch.  
Then run Mapper2.py.  

2. On a linus desktop, type.  
conda install -c conda-forge tensorflow  
Then open analysis.ipnb to analyze the mRNA display data.   

Typical install time on a "normal" desktop computer:    
<30 minutes.   

## 3. Demo  
Instructions to run on data.   
On a Sun Grid Engine, split sequencing file into chunks, then use mywrapper.sh to run Mapper1.py in batch.    
After Mapper1.py is done, run Mapper2.py to summarize reads.    
Then download readcount to local desktop.    
Use Jupyter notebook to open analysis.ipnb.    

Expected output:    
Mapper1.py maps mRNA display reads to reference sequences and save mutagenized region to temperary files.    
Mapper2.py summarizes mapped reads and count the occurance of mutations.    
analysis.ipnb analyzes mutation counts and make figures.    

Expected run time for demo on a "normal" desktop computer:    
Mapper1.py takes ~30min for each million reads.    
Mapper2.py takes ~10 hours to summarize ~600 million reads.   
analysis.ipnb takes ~2 hours to run all the analysis.   

## 4. Instructions for use
How to run the software on your data.   
Download data from:   
https://www.ncbi.nlm.nih.gov/bioproject/615649   
Use split.py to split the data into small chunks.   
Run codes as instructed above.   
Send emails to tianhao@ucla.edu for help.   
