#!/usr/bin/env python
# coding: utf-8

# Command line tools for creating .fastq files with varrying amounts of convolution,
# to test scRNA-seq based deconvolution methods

# In[1]:


import os
import sys
import subprocess
from sys import argv


# In[2]:


def count_fastq_reads(fastq_file_path: str):
    '''
    Returns the number of reads in a fastq or fastq.qz file
    given by fastq_file_path
    '''
    # subprocess.Popen(["cat", "/etc/services"], stdout=subprocess.PIPE, shell=True)
    if fastq_file_path[-2::] == 'gz':
        qs = "zcat {} | wc -l".format(fastq_file_path)
        n_lines = subprocess.run([qs],stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 shell=True)
    else:
        qs = "cat {} | wc -l".format(fastq_file_path)
        n_lines = subprocess.run([qs],stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 shell=True)
    
    n_lines = int(n_lines.stdout)
    n_reads = int(n_lines) / 4
    return n_reads


# In[3]:


def convolute_raw_data(sampleA: str, sampleB: str,
                       out_path: str, fractionContam: float,
                       seed=42,
                       **kwargs):
    """
    A function to create a convoluted mixture of reads from sampleA
    and SampleB. sampleA and sampleB are file paths to forward or 
    reverse reads for each sample.
    The output file will contain the same number of reads as sampleA. 
    fractionContam indicates the proportion of reads from sampleA 
    that will be replaced with reads from sampleB.
    
    Also, one can pass the number of reads in sampleA and sampleB
    as kwargs sampleA_reads and sampleB_reads if those quantities 
    are known, which will improve performance somewhat. 
    
    out_path indicates where the convoluted file will be written. 
    """
    # get the number of reads in sampleA
    if 'sampleA_reads' not in kwargs and 'sampleB_reads' not in kwargs:
        sample_A_reads = count_fastq_reads(sampleA)
        sample_B_reads = count_fastq_reads(sampleB)
    else:
        sample_A_reads = kwargs['sampleA_reads']
        sample_B_reads = kwargs['sampleB_reads']
        
    # calc the number of reads needed for each new file
    n_contam_reads = sample_A_reads * fractionContam
    n_sampleA_reads = sample_A_reads - n_contam_reads
    n_contam_reads = int(n_contam_reads)
    n_sampleA_reads = int(n_sampleA_reads)
    # rounding can make end file not exactly same as 
    # original sampleA
    # write the convoluted files
    cmd_1 = "seqtk sample -s {} {} {} > {}"
    cmd_1 = cmd_1.format(seed, sampleA, n_sampleA_reads, out_path)
    # command to radd contaminating reads
    cmd_2 = "seqtk sample -s {} {} {} >> {}"
    cmd_2 = cmd_2.format(seed, sampleB, n_contam_reads, out_path)
    # execute the commands
    subprocess.run([cmd_1], shell=True)
    subprocess.run([cmd_2], shell=True)


# In[4]:


# seqtk sample -s 123 read1.fq 100000 > sub_read1.fq


# In[5]:


#convolute_raw_data(sampleA='spam', sampleB='eggs', out_path='/', fractionContam=42)


# #fractionContam indicates the proportion of reads 
