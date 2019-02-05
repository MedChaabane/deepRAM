import argparse
import gzip
import csv
import random



def dinucshuffle(sequence):
    b=[sequence[i:i+2] for i in range(0, len(sequence), 2)]
    random.shuffle(b)
    d=''.join([str(x) for x in b])
    return d
    
def convert_format(parser):
	
	ENCODE_data = parser.ENCODE_data
	output = parser.output
	with gzip.open(output, 'wb') as fw:
		fw.write('sequence	label\n'.encode('utf-8'))
		with gzip.open(ENCODE_data, 'rt') as data:
			 next(data)
			 reader = csv.reader(data,delimiter='\t')
			 for row in reader:
				 
				 fw.write(row[2].encode('utf-8'))
				 fw.write('	'.encode('utf-8'))
				 fw.write('1'.encode('utf-8'))
				 fw.write('\n'.encode('utf-8'))
				 fw.write(dinucshuffle(row[2]).encode('utf-8'))
				 fw.write('	'.encode('utf-8'))
				 fw.write('0'.encode('utf-8'))
				 fw.write('\n'.encode('utf-8'))
	
	              
                      
def parse_arguments(parser):
## data
    parser.add_argument('--ENCODE_data', type=str, default='/s/chopin/a/grad/amenit/Desktop/ENCODE-models/*SRF_H1-hESC_SRF_HudsonAlpha_AC.seq.gz/SRF_H1-hESC_SRF_HudsonAlpha_AC.seq.gz',
                        help='path for ENCODE chip-seq data ')
    
    parser.add_argument('--output', type=str, default='train_ChIP.gz',
                        help='path for output file that matches format needed for deepRAM ')

    args = parser.parse_args()

    return args
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='preprocessing ChIP-seq file')
    args = parse_arguments(parser)
    convert_format(args)
