import argparse
import gzip
import csv
import random
import numpy as np



def read_seq(seq_file):
    seq_list = []
    seq = ''
    with gzip.open(seq_file, 'rt') as fp:
        for line in fp:
            
            if line[0] == '>':
                name = line[1:-1]
                if len(seq):
                    seq_list.append(seq)                    
                seq = ''
            else:
                seq = seq + line[:-1]
        if len(seq):
            seq_list.append(seq) 
    
    return np.array(seq_list)
def load_label_seq(seq_file):
    label_list = []
    seq = ''
    with gzip.open(seq_file, 'rt') as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:-1]
                posi_label = name.split(';')[-1]
                label = posi_label.split(':')[-1]
                label_list.append(int(label))
    return np.array(label_list)
    
def convert_format(parser):
	
	CLIP_data = parser.CLIP_data
	output = parser.output
	s=read_seq(CLIP_data)
	l=load_label_seq(CLIP_data)  
	with gzip.open(output, 'wb') as fw:
		fw.write('sequence	label\n'.encode('utf-8'))
		
		for i in range(len(s)):
			
			fw.write(s[i].encode('utf-8'))
			fw.write('	'.encode('utf-8'))
			fw.write(str(l[i]).encode('utf-8'))
			fw.write('\n'.encode('utf-8'))
			
                      
def parse_arguments(parser):
## data
    parser.add_argument('--CLIP_data', type=str, default='/s/chopin/a/grad/amenit/Downloads/DeepBind/clip/2_PARCLIP_AGO2MNASE_hg19/30000/test_sample_0/sequences.fa.gz',
                        help='path for CLIP-seq data ')
    
    parser.add_argument('--output', type=str, default='train_CLIP.gz',
                        help='path for output file that matches format needed for deepRAM ')

    args = parser.parse_args()

    return args
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='preprocessing CLIP-seq file')
    args = parse_arguments(parser)
    convert_format(args)
