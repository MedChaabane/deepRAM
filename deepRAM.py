

import csv
import math 
import random
import gzip
import torch
from sklearn import metrics
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence
import gensim
import multiprocessing
import numpy as np
from torch.utils.data import Dataset, DataLoader
from extract_motifs import get_motif
import torch.nn as nn
import os
import argparse
import warnings


warnings.filterwarnings("ignore")
# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(device)
import torch.nn.functional as F


def seqtopad(sequence,motlen):
    rows=len(sequence)+2*motlen-2
    S=np.empty([rows,4])
    base= bases if data_type=='DNA' else basesRNA
    for i in range(rows):
        for j in range(4):
            if i-motlen+1<len(sequence) and sequence[i-motlen+1]=='N' or i<motlen-1 or i>len(sequence)+motlen-2:
                S[i,j]=np.float32(0.25)
            elif sequence[i-motlen+1]==base[j]:
                S[i,j]=np.float32(1)
            else:
                S[i,j]=np.float32(0)
    return np.transpose(S)
def dinucshuffle(sequence):
    b=[sequence[i:i+2] for i in range(0, len(sequence), 2)]
    random.shuffle(b)
    d=''.join([str(x) for x in b])
    return d

def logsampler(a,b):
        x=np.random.uniform(low=0,high=1)
        y=10**((math.log10(b)-math.log10(a))*x + math.log10(a))
        return y
    
def sqrtsampler(a,b):
        
        x=np.random.uniform(low=0,high=1)
        y=(b-a)*math.sqrt(x)+a
        return y
      
# input of shape(batch_size,inp_chan,iW)

class Network(nn.Module):
    def __init__(self, nummotif,motiflen,RNN_hidden_size,hidden_size,hidden,dropprob,sigmaConv,sigmaNeu,sigmaRNN,xavier_init):
      
        super(Network, self).__init__()
        
        self.hidden=hidden
        self.RNN_hidden_size=RNN_hidden_size
        
        self.dropprob=dropprob
        self.sigmaConv=sigmaConv
        self.sigmaNeu=sigmaNeu
        self.hidden_size=hidden_size  
        self.input_channels=4
        
        # Embedding
        if embedding:
              model1 = gensim.models.Word2Vec.load(word2vec_model)
              weights = torch.FloatTensor(model1.wv.vectors)
              self.embedding = nn.Embedding.from_pretrained(weights, freeze=False)
              self.input_channels=Embsize
              
        # Convnet
        self.ConvWeights=[]
        self.ConvBias=[]
        self.wConv=torch.randn(nummotif,self.input_channels,motiflen).to(device)
        self.wRect=torch.randn(nummotif).to(device)
        self.ConvWeights.append(self.wConv)
        self.ConvBias.append(self.wRect)
        conv_channels=nummotif
        
        if conv:
              self.FC_size= nummotif
              self.input_channels=nummotif
              for c in range(1,conv_layers):
                  Wconv=torch.randn(int(1.5*c*nummotif),conv_channels,motiflen).to(device)
                  Bconv=torch.randn(int(1.5*c*nummotif)).to(device)
                  self.ConvWeights.append(Wconv)
                  self.ConvBias.append(Bconv)
                  conv_channels=int(1.5*c*nummotif)
                  self.FC_size= int(1.5*c*nummotif)
                  self.input_channels=int(1.5*c*nummotif)
              torch.nn.init.normal_(self.ConvWeights[0],mean=0,std=sigmaConv)
              
              ind=0
              for weights in self.ConvWeights:
                  weights.requires_grad=True
                  if ind>0:
                    if dilation>1:
                      torch.nn.init.normal_(weights,mean=0,std=0.1)
                      print('ffffffff')
                    else:
                      torch.nn.init.xavier_uniform(weights)
                  ind=ind+1
              for weights in self.ConvBias:
                  weights.requires_grad=True
                  torch.nn.init.normal_(weights)
             

              
        
        # RNN
        self.rnn = nn.GRU(self.input_channels, RNN_hidden_size, num_layers=1, bidirectional=False).to(device)
        if RNN:
              if RNN_type=='GRU':
                  self.rnn = nn.GRU(self.input_channels, RNN_hidden_size, num_layers=RNN_layers, bidirectional=False).to(device)
                  self.FC_size= RNN_hidden_size
              elif RNN_type=='BiGRU':
                  self.rnn = nn.GRU(self.input_channels, RNN_hidden_size, num_layers=RNN_layers, bidirectional=True).to(device)
                  self.FC_size= 2*RNN_hidden_size
              elif RNN_type=='LSTM':
                  self.rnn = nn.LSTM(self.input_channels, RNN_hidden_size, num_layers=RNN_layers, bidirectional=False).to(device)
                  self.FC_size= RNN_hidden_size
              elif RNN_type=='BiLSTM':
                  self.rnn = nn.LSTM(self.input_channels, RNN_hidden_size, num_layers=RNN_layers, bidirectional=True).to(device)
                  self.FC_size= 2*RNN_hidden_size
              if not xavier_init:
                for layer_p in self.rnn._all_weights:
                  for p in layer_p:
                    if 'weight' in p:
                      torch.nn.init.normal_(self.rnn.__getattr__(p),mean=0,std=sigmaRNN)
              else:
                for layer_p in self.rnn._all_weights:
                  for p in layer_p:
                    if 'weight' in p:
                      torch.nn.init.xavier_uniform(self.rnn.__getattr__(p))

        # FC    
        self.wHidden=torch.randn(self.FC_size,self.hidden_size).to(device)
        self.wHiddenBias=torch.randn(self.hidden_size).to(device)
        self.wHidden.requires_grad=True
        self.wHiddenBias.requires_grad=True
        
        if not self.hidden:
            self.wNeu=torch.randn(self.FC_size,1).to(device)
            self.wNeuBias=torch.randn(1).to(device)
            if not xavier_init:
              torch.nn.init.normal_(self.wNeu,mean=0,std=self.sigmaNeu)
              torch.nn.init.normal_(self.wNeuBias,mean=0,std=self.sigmaNeu)
            else:
              torch.nn.init.xavier_uniform(self.wNeu)
            
            

        else:
           
            self.wNeu=torch.randn(self.hidden_size,1).to(device)
            self.wNeuBias=torch.randn(1).to(device)     
            if not xavier_init:  
              torch.nn.init.normal_(self.wNeu,mean=0,std=self.sigmaNeu)
              torch.nn.init.normal_(self.wNeuBias,mean=0,std=self.sigmaNeu)
              torch.nn.init.normal_(self.wHidden,mean=0,std=self.sigmaNeu)
              torch.nn.init.normal_(self.wHiddenBias,mean=0,std=self.sigmaNeu)
            else:
              torch.nn.init.xavier_uniform(self.wNeu)
              torch.nn.init.xavier_uniform(self.wHidden)
            
  
           

        self.wNeu.requires_grad=True
        self.wNeuBias.requires_grad=True
        
        self.dropout=torch.nn.Dropout(p=dropprob, inplace=False)
        self.max=torch.nn.MaxPool1d(3, stride=1)
        
        
    def get_weights(self):
        ll=[]
        for layer_p in self.rnn._all_weights:
          for p in layer_p:
            if 'weight' in p:
               ll.append(self.rnn.__getattr__(p))
        return ll+self.ConvWeights+self.ConvBias
        
    def layer1out(self,x):
		
        if type(x) is np.ndarray:
          x = torch.from_numpy(x.astype(np.float32))
        #x = Variable(x, volatile=True)
        if torch.cuda.is_available():
          x = x.to(device)
        if embedding:
          print(x.shape)
          x= self.embedding(x)
          x=x.permute(0,2,1)
          print(x.shape)
        if conv:
          x=F.conv1d(x, self.wConv, bias=self.wRect, stride=1, padding=0)
          out=x.clamp(min=0)
          print(out.shape)
          temp = out.data.cpu().numpy()
        else:
          print('you need to have CNN to visualize motifs')
        return temp
        
    def forward(self, x):
      
        if embedding:
          # shape of x : batch_size x seq_len
          x= self.embedding(x)
          x=x.permute(0,2,1)
          # shape of x_emb : batch_size x embd_size x seq_len
          
#         else:
#           # shape of x : batch_size x 4 x seq_len
         
        if conv:
          x=F.conv1d(x, self.ConvWeights[0], bias=self.ConvBias[0], stride=1, padding=0)
          x=x.clamp(min=0)
          x=self.max(x)
          
          for c in range(1,len(self.ConvWeights)):
            x=F.conv1d(x, self.ConvWeights[c], bias=self.ConvBias[c], stride=1, padding=0,dilation=dilation)
            x=x.clamp(min=0)
            x=self.max(x)
            
            
             
          
			
        if RNN:
          if conv:
            
            x=self.dropout(x)
          x=x.permute(2,0,1)
          # shape of x :  seq_len x batch_size  x features
          output, _ = self.rnn(x)
          # shape of output :  seq_len x batch_size  x num_directions * features
          
          if RNN_type== 'BiLSTM' or RNN_type=='BiGRU':

              Normal_RNN=output[-1, :, :self.RNN_hidden_size]
              Rev_RNN=output[0, :, self.RNN_hidden_size:]
              x = torch.cat((Normal_RNN, Rev_RNN), 1)
              x=self.dropout(x)
              #shape of x: batch_size x 2*hidden_size
#               print(x.shape)
              
          else:
                      ## from (1, N, hidden) to (N, hidden)
              x = output[-1, :, :]
              x=self.dropout(x)
#               print(hn.shape)
#               x = hn.view(hn.size()[1], hn.size(2))
              # shape of x: batch_size x hidden_size
              #print(x.shape)
          
          
        else:
          x, _ = torch.max(x, dim=2)
          #print(x.shape)

          # shape of x : batch_size x numb_filters
          x=self.dropout(x)
        if self.hidden:
          x=x @ self.wHidden
          x.add_(self.wHiddenBias)
          x=x.clamp(min=0)
          x=self.dropout(x)
        x=x @ self.wNeu
        x.add_(self.wNeuBias)
          
        
        
        return torch.sigmoid(x)



            




class Chip():
    def __init__(self,filename,motiflen=24):
        self.file = filename
        self.motiflen = motiflen
  
    def openFile(self):
        train_dataset=[]
        sequences=[]
        with gzip.open(self.file, 'rt') as data:
                next(data)
                reader = csv.reader(data,delimiter='\t')

                if embedding:

                        for row in reader:

                            ## When using Embedding
                            sequences.append(row[0])
                            
                            train_dataset.append([row[0],[int(row[1])]])
                else:
                        for row in reader:
                                
                                train_dataset.append([seqtopad(row[0],self.motiflen),[int(row[1])]])
  
        random.shuffle(train_dataset)
        size=int(len(train_dataset)/3)
        firstvalid=train_dataset[:size]
        secondvalid=train_dataset[size:size+size]
        thirdvalid=train_dataset[size+size:]
        firsttrain=secondvalid+thirdvalid
        secondtrain=firstvalid+thirdvalid
        thirdtrain=firstvalid+secondvalid
        return firsttrain,firstvalid,secondtrain,secondvalid,thirdtrain,thirdvalid,train_dataset,sequences
        
def Gen_Words(sequences,kmer_len,s):
    out=[]

    for i in sequences:

        kmer_list=[]
        for j in range(0,(len(i)-kmer_len)+1,s):

              kmer_list.append(i[j:j+kmer_len])

        out.append(kmer_list)

    return out



class chipseq_dataset(Dataset):
    """ Diabetes dataset."""

    def __init__(self,xy=None):
        self.x_data=np.asarray([el[0] for el in xy],dtype=np.float32)
        self.y_data =np.asarray([el[1] for el in xy ],dtype=np.float32)
        self.x_data = torch.from_numpy(self.x_data)
        self.y_data = torch.from_numpy(self.y_data)
        self.len=len(self.x_data)
      

    def __getitem__(self, index):
        return self.x_data[index], self.y_data[index]

    def __len__(self):
        return self.len

class chipseq_dataset_embd(Dataset):
    """ Diabetes dataset."""

    def __init__(self,xy=None,model=None,kmer_len=5,stride=2):
      
        self.kmer_len= kmer_len
        self.stride= stride
        data=[el[0] for el in xy]
        words_doc= self.Gen_Words(data,self.kmer_len,self.stride)
#         print(words_doc[0])
        x_data=[self.convert_data_to_index(el,model.wv) for el in words_doc]
#         print(x_data[0])
       
        
        self.x_data=np.asarray(x_data,dtype=np.float32)
        self.y_data =np.asarray([el[1] for el in xy ],dtype=np.float32)
        self.x_data = torch.LongTensor(self.x_data)
        self.y_data = torch.from_numpy(self.y_data)
        self.len=len(self.x_data)
      

    def __getitem__(self, index):
        return self.x_data[index], self.y_data[index]

    def __len__(self):
        return self.len
      
    def Gen_Words(self,pos_data,kmer_len,s):
        out=[]
        
        for i in pos_data:

            kmer_list=[]
            for j in range(0,(len(i)-kmer_len)+1,s):

                  kmer_list.append(i[j:j+kmer_len])
                
            out.append(kmer_list)
            
        return out


    def convert_data_to_index(self, string_data, wv):
      index_data = []
      for word in string_data:
          if word in wv:
              index_data.append(wv.vocab[word].index)
      return index_data
      
      
class Chip_test():
    def __init__(self,filename,motiflen=24):
        self.file = filename
        self.motiflen = motiflen
    def openFile(self):
        test_dataset=[]
        seq=[]


        with gzip.open(self.file, 'rt') as data:
            next(data)
            reader = csv.reader(data,delimiter='\t')
            if embedding:
                if evaluate_performance: 
                    for row in reader:
						## When using Embedding
                        test_dataset.append([row[0],[int(row[1])]])
                        seq.append(row[0])
                else:
                    for row in reader:
						## just adding fake label but it will not be used
                        test_dataset.append([row[0],[1]])
                        seq.append(row[0])
						
            else:
                if evaluate_performance: 
                    for row in reader:
                        test_dataset.append([seqtopad(row[0],self.motiflen),[int(row[1])]])
                        seq.append(row[0])
                else:
                    for row in reader:
						## just adding fake label but it will not be used
                        test_dataset.append([seqtopad(row[0],self.motiflen),[1]])
                        seq.append(row[0])
        return test_dataset,seq    



train_dataloader=[]
valid_dataloader=[]
test_loader=[]
sequences=[]
seq=[]
def Load_Data(train_file,test_file):
	global nummotif
	global motiflen
	global train_dataloader
	global valid_dataloader
	global test_loader
	global alldataset_loader
	global sequences
	global seq
	global motif_loader
	global seq_motif
	
	if embedding:
		nummotif=32 #number of motifs to discover
		motiflen=10
	if train:
		chipseq=Chip(train_file)
		train1,valid1,train2,valid2,train3,valid3,alldataset,sequences=chipseq.openFile()

		#### word2vect model training


		print(embedding)
		if embedding and word2vect_train:
			print('training word2vec model')
			document= Gen_Words(sequences,kmer_len,stride)
			model = gensim.models.Word2Vec (document, window=int(12 / stride), min_count=0, size=Embsize,workers=multiprocessing.cpu_count())
			model.train(document,total_examples=len(document),epochs=Embepochs)
			model.save(word2vec_model)





		if embedding:
			model1 = gensim.models.Word2Vec.load(word2vec_model)
			train1_dataset=chipseq_dataset_embd(train1,model1,kmer_len,stride)
			train2_dataset=chipseq_dataset_embd(train2,model1,kmer_len,stride)
			train3_dataset=chipseq_dataset_embd(train3,model1,kmer_len,stride)
			valid1_dataset=chipseq_dataset_embd(valid1,model1,kmer_len,stride)
			valid2_dataset=chipseq_dataset_embd(valid2,model1,kmer_len,stride)
			valid3_dataset=chipseq_dataset_embd(valid3,model1,kmer_len,stride)
			alldataset_dataset=chipseq_dataset_embd(alldataset,model1,kmer_len,stride)
		else:
			train1_dataset=chipseq_dataset(train1)
			train2_dataset=chipseq_dataset(train2)
			train3_dataset=chipseq_dataset(train3)
			valid1_dataset=chipseq_dataset(valid1)
			valid2_dataset=chipseq_dataset(valid2)
			valid3_dataset=chipseq_dataset(valid3)
			alldataset_dataset=chipseq_dataset(alldataset)

		train_loader1 = DataLoader(dataset=train1_dataset,batch_size=batch_size,shuffle=True)
		train_loader2 = DataLoader(dataset=train2_dataset,batch_size=batch_size,shuffle=True)
		train_loader3 = DataLoader(dataset=train3_dataset,batch_size=batch_size,shuffle=True)
		valid1_loader = DataLoader(dataset=valid1_dataset,batch_size=batch_size,shuffle=True)
		valid2_loader = DataLoader(dataset=valid2_dataset,batch_size=batch_size,shuffle=True)
		valid3_loader = DataLoader(dataset=valid3_dataset,batch_size=batch_size,shuffle=True)
		alldataset_loader=DataLoader(dataset=alldataset_dataset,batch_size=batch_size,shuffle=True)

		train_dataloader=[train_loader1,train_loader2,train_loader3]
		valid_dataloader=[valid1_loader,valid2_loader,valid3_loader]

	#### test dataset

	if embedding:
		model1 = gensim.models.Word2Vec.load(word2vec_model)
	chipseq_test=Chip_test(test_file)   
	motif_test=Chip_test(test_file,1)  
	
	test_data, seq=chipseq_test.openFile()
	motif_data, seq_motif=motif_test.openFile()
	if embedding:
		test_dataset=chipseq_dataset_embd(test_data,model1,kmer_len,stride)
		motif_dataset=chipseq_dataset_embd(motif_data,model1,kmer_len,stride)
	else:
		test_dataset=chipseq_dataset(test_data)
		motif_dataset=chipseq_dataset(motif_data)
	
	test_loader = DataLoader(dataset=test_dataset,batch_size=128,shuffle=True)
	motif_loader = DataLoader(dataset=motif_dataset,batch_size=10000000,shuffle=False)

	

def Calibration():
	print('start')
	best_AUC=0
	device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
	print(device)
	# device='cpu'
	learning_steps_list=[5000,10000,15000,20000,25000,30000,35000,40000]
	for number in range(40):
		# hyper-parameters
		RNN_hidden_size_list=[20,50,80,100]
		RNN_hidden_size=random.choice(RNN_hidden_size_list)
		dropoutList=[0,0.15,0.3,0.45,0.6] 
		dropprob=random.choice(dropoutList)
		hidden_list=[True,False]
		hidden=random.choice(hidden_list)
		xavier_List=[True,True,False] 
		xavier=random.choice(xavier_List)
		hidden_size_list=[32,64]
		hidden_size=random.choice(hidden_size_list)
		optim_list=['SGD','Adagrad','Adagrad']
		optim=random.choice(optim_list)
		learning_rate=logsampler(0.005,0.5) 
		momentum_rate=sqrtsampler(0.95,0.99)  
		sigmaConv=logsampler(10**-6,10**-2)   
		sigmaNeu=logsampler(10**-3,10**-1) 
		sigmaRNN=logsampler(10**-4,10**-1) 
		weightDecay=logsampler(10**-10,10**-1) 
		nummotif_list=[16]
		nummotif1=random.choice(nummotif_list)
		
		
		
	
		
		model_auc=[[],[],[]]
		for kk in range(3):
			model = Network(nummotif1,motiflen,RNN_hidden_size,hidden_size,hidden,dropprob,sigmaConv,sigmaNeu,sigmaRNN,xavier).to(device)
			if optim=='SGD':
				optimizer = torch.optim.SGD(model.get_weights()+[model.wNeu,model.wNeuBias,model.wHidden,model.wHiddenBias], lr=learning_rate,momentum=momentum_rate,nesterov=True
											,weight_decay=weightDecay)
			else:
				optimizer = torch.optim.Adagrad(model.get_weights()+[model.wNeu,model.wNeuBias,model.wHidden,model.wHiddenBias], lr=learning_rate,weight_decay=weightDecay)

			train_loader=train_dataloader[kk]
			valid_loader=valid_dataloader[kk]
			
			
			learning_steps=0
			while learning_steps<=40000:
			   
				auc=[]
				model.train()
				for i, (data, target) in enumerate(train_loader):
					
					data = data.to(device)
					target = target.to(device)

					# Forward pass
					output = model(data)          
					loss = F.binary_cross_entropy(output,target)
					optimizer.zero_grad()
					loss.backward()
					optimizer.step()
					learning_steps+=1

					if learning_steps% 5000==0:
						
							
							
							
						with torch.no_grad():
							model.eval()
							auc=[]
							for j, (data1, target1) in enumerate(valid_loader):
								data1 = data1.to(device)
								target1 = target1.to(device)
								
								# Forward pass
								output = model(data1)
								
								pred=output.cpu().detach().numpy().reshape(output.shape[0])
								labels=target1.cpu().numpy().reshape(output.shape[0])
								if output.shape[0]>60:
									auc.append(metrics.roc_auc_score(labels, pred))
							#print(np.mean(auc))
							model_auc[kk].append(np.mean(auc))
							
							model.train()
						
		
		print('                   ##########################################               ')

		for n in range(8):
			AUC=(model_auc[0][n]+model_auc[1][n]+model_auc[2][n])/3
			#print(AUC)
			if AUC>best_AUC:
				best_AUC=AUC
				best_learning_steps=learning_steps_list[n]
				best_LearningRate=learning_rate
				best_LearningMomentum=momentum_rate
				best_sigmaConv=sigmaConv
				best_dropprob=dropprob
				best_sigmaNeu=sigmaNeu
				best_RNN_hidden_size=RNN_hidden_size
				best_weightDecay=weightDecay
				best_hidden=hidden
				best_sigmaRNN=sigmaRNN
				best_xavier=xavier
				best_optim=optim
				best_nummotif=nummotif
				best_hidden_size=hidden_size


	print('best_AUC=',best_AUC)            
	print('best_learning_steps=',best_learning_steps)      
	print('best_LearningRate=',best_LearningRate)
	print('best_LearningMomentum=',best_LearningMomentum)
	print('best_sigmaConv=',best_sigmaConv)
	print('best_dropprob=',best_dropprob)
	print('best_sigmaNeu=',best_sigmaNeu)
	print('best_RNN_hidden_size',best_RNN_hidden_size)
	print('best_weightDecay=',weightDecay)
	print('best_hidden=',best_hidden)
	print('best_sigmaRNN=',best_sigmaRNN)
	print('best_xavier=',best_xavier)
	print('best_optim=',best_optim)
	print('best_nummotif=',best_nummotif)
	print('best_hidden_size=',best_hidden_size)
	
	best_hyperparameters = {'best_learning_steps': best_learning_steps,'best_LearningRate':best_LearningRate,'best_LearningMomentum':best_LearningMomentum,'best_sigmaConv':best_sigmaConv,
	                         'best_dropprob':best_dropprob,'best_sigmaNeu':best_sigmaNeu,'best_RNN_hidden_size':best_RNN_hidden_size,
	                         'best_weightDecay':best_weightDecay,'best_hidden':best_hidden,'best_sigmaRNN':best_sigmaRNN,'best_xavier':best_xavier,'best_optim':best_optim,'best_nummotif':best_nummotif,'best_hidden_size':best_hidden_size}
	torch.save(best_hyperparameters, model_dir+'best_hyperpamarameters.pth')
	return best_hyperparameters
def Train_model():
	best_hyperparameters=torch.load(model_dir+'best_hyperpamarameters.pth')
	best_learning_steps=best_hyperparameters['best_learning_steps']
	best_LearningRate=best_hyperparameters['best_LearningRate']
	best_LearningMomentum=best_hyperparameters['best_LearningMomentum']
	best_sigmaConv=best_hyperparameters['best_sigmaConv']
	best_dropprob=best_hyperparameters['best_dropprob']
	best_sigmaNeu=best_hyperparameters['best_sigmaNeu']
	best_RNN_hidden_size=best_hyperparameters['best_RNN_hidden_size']
	best_weightDecay=best_hyperparameters['best_weightDecay']
	best_hidden=best_hyperparameters['best_hidden']
	best_sigmaRNN=best_hyperparameters['best_sigmaRNN']
	best_xavier=best_hyperparameters['best_xavier']
	best_optim=best_hyperparameters['best_optim']
	best_nummotif=best_hyperparameters['best_nummotif']
	best_hidden_size=best_hyperparameters['best_hidden_size']
	best_AUC=0


	for number_models in range(5):
	  model = Network(best_nummotif,motiflen,best_RNN_hidden_size,best_hidden_size,best_hidden,best_dropprob,best_sigmaConv,best_sigmaNeu,best_sigmaRNN,best_xavier).to(device)
	  if best_optim=='SGD':
		  optimizer = torch.optim.SGD(model.get_weights()+[model.wNeu,model.wNeuBias,model.wHidden,model.wHiddenBias], lr=best_LearningRate,momentum=best_LearningMomentum,nesterov=True,weight_decay=best_weightDecay)
	  else:
		  optimizer = torch.optim.Adagrad(model.get_weights()+[model.wNeu,model.wNeuBias,model.wHidden,model.wHiddenBias], lr=best_LearningRate,weight_decay=best_weightDecay)


	  train_loader=alldataset_loader
	  valid_loader=alldataset_loader
	  learning_steps=0
	  model.train()
	  while learning_steps<=best_learning_steps:
	  
		  for i, (data, target) in enumerate(train_loader):
			  data = data.to(device)
			  target = target.to(device)
			  
				# Forward pass
			  output = model(data)
			  loss = F.binary_cross_entropy(output,target)

			  optimizer.zero_grad()
			  loss.backward()
			  optimizer.step()
			  learning_steps+=1
			  
	  with torch.no_grad():
		  model.eval()
		  auc=[]
		  for i, (data, target) in enumerate(valid_loader):
			  data = data.to(device)
			  target = target.to(device)
			  
			  # Forward pass
			  output = model(data)
			  
			  pred=output.cpu().detach().numpy().reshape(output.shape[0])
			  labels=target.cpu().numpy().reshape(output.shape[0])
			  if output.shape[0]>30:
				  auc.append(metrics.roc_auc_score(labels, pred))
	  #             
		  AUC_training=np.mean(auc)
		  print('AUC on training data for model ',number_models+1,' = ',AUC_training)
		  if AUC_training>best_AUC:
			  best_AUC=AUC_training
			  best_model=model
	torch.save(best_model, model_path)
	#torch.save(best_model.state_dict(), model_dir+'best_model.pkl')
	return best_model
#### save model .pkl
#### load model
def Test_Motifs():
	

	model = torch.load(model_path)
	with torch.no_grad():
		model.eval()
		for i, (data, target) in enumerate(motif_loader):
			seqRNA= seq_motif
			if data_type=='RNA':
				seqRNA=[sequ.replace('T','U') for sequ in seqRNA]
			detect_motifs(model,seqRNA , data, motif_dir)

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'
def detect_motifs(model, test_seqs, X_train, output_dir):


	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
		
	for param in model.parameters():
		layer1_para =  param.data.cpu().numpy()
		break
		
	layer1_para=model.wConv.data.cpu().numpy()

	filter_outs = model.layer1out(X_train)
	get_motif(layer1_para, filter_outs, test_seqs, dir1 = output_dir,embd=embedding,data=data_type,kmer=kmer_len,s=stride,tomtom=tomtom_dir)

def test_predict():
	

	model = torch.load(model_path)
	
	

	with torch.no_grad():
		  model.eval()
		  auc=[]
		 
		  for i, (data, target) in enumerate(test_loader):
			  data = data.to(device)
			  target = target.to(device)
			  
			 
			  # Forward pass
			  output = model(data)
			  pred=output.cpu().detach().numpy().reshape(output.shape[0])
			  
			  fw = open(out_file, 'w')
			  myprob = "\n".join(map(str, pred[:]))
			
			  fw.write(myprob)
			  
			
			  labels=target.cpu().numpy().reshape(output.shape[0])
			  if evaluate_performance:
				  if output.shape[0]>50:
					  auc.append(metrics.roc_auc_score(labels, pred))
		  if evaluate_performance:                     
			  AUC_test=np.mean(auc)
			  print('AUC on test data = ',AUC_test)
			  fw.write('\nAUC on test data = ')
			  fw.write(str(AUC_test))
	fw.close()
	
######### Global variables #########
embedding=False
conv=True
RNN=False
RNN_type='BiLSTM'

bases='ACGT' #DNA bases
basesRNA='ACGU'#RNA bases
batch_size=128
evaluate_performance=False
train=True
model_dir='models/'
# embedding hyper-parameters
Embepochs=100
Embsize=50
kmer_len=3
stride=1
word2vect_train=True
word2vec_model='models/word2vec_model'
# CNN hyper-parameters
nummotif=16 #number of motifs to discover
motiflen=24

################################

	  

def run_deepRAM(parser):
	global embedding
	global conv
	global RNN
	global RNN_type
	global kmer_len
	global stride
	global word2vec_train 
	global word2vec_model 
	global evaluate_performance
	global train
	global model_dir
	global model_path
	global out_file
	global motif
	global motif_dir
	global tomtom_dir
	global data_type
	global conv_layers
	global RNN_layers
	global dilation
	
	train_data = parser.train_data
	test_data = parser.test_data
	data_type=parser.data_type
	train = parser.train
	predict=parser.predict_only
	model_dir = parser.models_dir
	model_path=parser.model_path
	out_file = parser.out_file
	motif=parser.motif
	motif_dir=parser.motif_dir
	tomtom_dir=parser.tomtom_dir
	evaluate_performance=parser.evaluate_performance
	embedding = parser.Embedding
	conv = parser.Conv
	RNN = parser.RNN
	RNN_type = parser.RNN_type
	dilation=parser.dilation
	kmer_len = parser.kmer_len
	stride = parser.stride
	word2vec_train = parser.word2vec_train
	word2vec_model = parser.word2vec_model
	conv_layers=parser.conv_layers
	RNN_layers=parser.RNN_layers
	print(embedding)
	if not os.path.exists(model_dir):
		  os.makedirs(model_dir)


	if predict:
		  train = False

	if train:
		
		  print('Load Data')
		  Load_Data(train_data,test_data)
		
		  print('Automatic Calibration')
		
		  best_hyperparameters=Calibration()
		  print('Training 6 models using best hyper-parameters set')
		  model=Train_model()
		  print('Predicting sequence specificities')
		  model=test_predict()
		  
		  

	else:
		  print('Load Data')
		  Load_Data(train_data,test_data)
		  print('Predicting sequence specificities')
		  model=test_predict()
	if motif:
		  Test_Motifs()
			  
                      
                      
def parse_arguments(parser):
## data
    parser.add_argument('--train_data', type=str, default='train.fa.gz',
                        help='path for training data with format: sequence 	label')
    
    parser.add_argument('--test_data', type=str, default='seq.fa.gz',
                        help='path for test data containing test sequences with or without label')
    parser.add_argument('--data_type', type=str, default='DNA',
                        help='type of data: DNA or RNA ')

## model
    parser.add_argument('--train', type=boolean_string, default=True, help='use this option for automatic calibration, training model using train_data and predict labels for test_data')
    parser.add_argument('--predict_only', type=boolean_string, default=False, help='use this option to load pretrained model (found in model_path) and use it to predict test sequences (train will be set to False).')
    parser.add_argument('--evaluate_performance', type=boolean_string, default=True, help='use this option to calculate AUC on test_data. If True, test_data should be format: sequence label')
    
    parser.add_argument('--models_dir', type=str, default='models/',
                        help='The directory to save the trained models for future prediction including best hyperparameters and embedding model')
    parser.add_argument('--model_path', type=str, default='DeepBind.pkl',
                        help='If train is set to True, This path will be used to save your best model. If train is set to False, this path should have the model that you want to use for prediction ')
    
    parser.add_argument('--motif', type=boolean_string, default=False, help='use this option to generate motif logos')

    parser.add_argument('--motif_dir', type=str, default='motifs',
                        help='directory to save motifs logos ')
    parser.add_argument('--tomtom_dir', type=str, default='motifs',
                        help='directory of TOMTOM, i.e:meme-5.0.3/src/tomtom')

    parser.add_argument('--out_file', type=str, default='prediction.txt',
                        help='The output file used to store the prediction probability of testing data')


## architecture
    parser.add_argument('--Embedding', type=boolean_string, default=False, help='Use embedding layer: True or False')

    parser.add_argument('--Conv', type=boolean_string, default=True, help='Use conv layer: True or False')


    parser.add_argument('--RNN', type=boolean_string, default=True, help='Use RNN layer: True or False')
    
    parser.add_argument('--RNN_type', type=str, default='BiLSTM', help='RNN type: LSTM or GRU or BiLSTM or BiGRU')
    
  
## Embedding 
    
    parser.add_argument('--kmer_len', type=int, default='3', help='length of kmer used for embedding layer, default=3')
    
    parser.add_argument('--stride', type=int, default='1', help='stride used for embedding layer, default=1')
    
    parser.add_argument('--word2vec_train', type=boolean_string, default=True, help='set it to False if you have already trained word2vec model. If you set it to False, you need to specify the path for word2vec model in word2vec_model argument.')
    
    parser.add_argument('--word2vec_model', type=str, default='word2vec', help='If word2vec_train is set to True, This path will be used to save your word2vec model. If word2vec_train is set to False, this path should have the word2vec model that you want to use for embedding layer')

    parser.add_argument('--conv_layers', type=int, default='1', help='number of convolutional modules')
    parser.add_argument('--dilation', type=int, default='1', help='the spacing between kernel elements for convolutional modules (except the first convolutional module)')
    parser.add_argument('--RNN_layers', type=int, default='1', help='number of RNN layers')
    
    args = parser.parse_args()

    return args
def main():
    parser = argparse.ArgumentParser(description='sequence specificities prediction using deep learning approach')
    args = parse_arguments(parser)
    run_deepRAM(args)    
if __name__ == "__main__":
    main()
