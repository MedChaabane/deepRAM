

   <img src="https://github.com/MedChaabane/deepRAM/blob/master/CSU-Ram.jpg" width="160">
   
  
# deepRAM

deepRAM is an end-to-end deep learning toolkit for predicting protein binding sites and motifs. It helps users run experiments using many state-of-the-art deep learning methods and addresses the challenge of selecting model parameters in deep learning models using a fully automatic model selection strategy. This helps avoid hand-tuning and thus removes any bias in running experiments, making it user friendly without losing its flexibility. While it was designed with ChIP-seq and CLIP-seq data in mind, it can be used for any DNA/RNA sequence binary classification problem.

deepRAM allows users the flexibility to choose a deep learning model by selecting its different components:  input sequence representation (one-hot or k-mer embedding), whether to use a CNN and how many layers, and whether to use an RNN, and the number of layers and their type. For CNNs, the user can choose to use dilated convolution as well.
 <br><br>
## Dependency <br>
We recommend to use [Anaconda 3](https://www.anaconda.com/download/) platform. 
-  <a href=https://www.python.org/downloads/source/>Python 3.6 </a> <br>
-  <a href=https://pytorch.org/>PyTorch 1.0 library </a> (Deep learning library) <br>
-  <a href=https://github.com/scikit-learn/scikit-learn>sklearn</a> (Machine learning library)<br>
-  <a href=https://anaconda.org/anaconda/gensim>gensim</a> (library used to train word2vec algorithm) <br>
-  <a href=https://anaconda.org/anaconda/numpy>numpy</a> <br>

## Usage

```bash 
usage: deepRAM.py [-h] [--train_data TRAIN_DATA] [--test_data TEST_DATA]
                  [--data_type DATA_TYPE] [--train TRAIN]
                  [--predict_only PREDICT_ONLY]
                  [--evaluate_performance EVALUATE_PERFORMANCE]
                  [--models_dir MODELS_DIR] [--model_path MODEL_PATH]
                  [--motif MOTIF] [--motif_dir MOTIF_DIR]
                  [--tomtom_dir TOMTOM_DIR] [--out_file OUT_FILE]
                  [--Embedding EMBEDDING] [--Conv CONV] [--RNN RNN]
                  [--RNN_type RNN_TYPE] [--kmer_len KMER_LEN]
                  [--stride STRIDE] [--word2vec_train WORD2VEC_TRAIN]
                  [--word2vec_model WORD2VEC_MODEL]
                  [--conv_layers CONV_LAYERS] [--dilation DILATION]
                  [--RNN_layers RNN_LAYERS]

sequence specificities prediction using deep learning approach

optional arguments:
  -h, --help            show this help message and exit
  --train_data TRAIN_DATA
                        path for training data with format: sequence label
  --test_data TEST_DATA
                        path for test data containing test sequences with or
                        without label
  --data_type DATA_TYPE
                        type of data: DNA or RNA. default: DNA
  --train TRAIN         use this option for automatic calibration, training
                        model using train_data and predict labels for
                        test_data. default: True
  --predict_only PREDICT_ONLY
                        use this option to load pretrained model (found in
                        model_path) and use it to predict test sequences
                        (train will be set to False). default: False
  --evaluate_performance EVALUATE_PERFORMANCE
                        use this option to calculate AUC on test_data. If
                        True, test_data should be format: sequence label.
                        default: False
  --models_dir MODELS_DIR
                        The directory to save the trained models for future
                        prediction including best hyperparameters and
                        embedding model. default: models/
  --model_path MODEL_PATH
                        If train is set to True, This path will be used to
                        save your best model. If train is set to False, this
                        path should have the model that you want to use for
                        prediction. default: BestModel.pkl
  --motif MOTIF         use this option to generate motif logos. default:
                        False
  --motif_dir MOTIF_DIR
                        directory to save motifs logos. default: motifs
  --tomtom_dir TOMTOM_DIR
                        directory of TOMTOM, i.e:meme-5.0.3/src/tomtom
  --out_file OUT_FILE   The output file used to store the prediction
                        probability of testing data
  --Embedding EMBEDDING
                        Use embedding layer: True or False. default: False
  --Conv CONV           Use conv layer: True or False. default: True
  --RNN RNN             Use RNN layer: True or False. default: False
  --RNN_type RNN_TYPE   RNN type: LSTM or GRU or BiLSTM or BiGRU. default:
                        BiLSTM
  --kmer_len KMER_LEN   length of kmer used for embedding layer, default= 3
  --stride STRIDE       stride used for embedding layer, default= 1
  --word2vec_train WORD2VEC_TRAIN
                        set it to False if you have already trained word2vec
                        model. If you set it to False, you need to specify the
                        path for word2vec model in word2vec_model argument.
                        default: True
  --word2vec_model WORD2VEC_MODEL
                        If word2vec_train is set to True, This path will be
                        used to save your word2vec model. If word2vec_train is
                        set to False, this path should have the word2vec model
                        that you want to use for embedding layer. default:
                        word2vec
  --conv_layers CONV_LAYERS
                        number of convolutional modules. default= 1
  --dilation DILATION   the spacing between kernel elements for convolutional
                        modules (except the first convolutional module).
                        default= 1
  --RNN_layers RNN_LAYERS
                        number of RNN layers. default= 1
```

## Motifs identification and visualization

You need to install <a href=http://weblogo.berkeley.edu/> WebLogo </a> and TOMTOM in <a href=http://meme-suite.org> MEME Suite </a> to match identifyed motifs with known motifs of Transcription Factors and RBPs. Read documentations about installation and usage.

## Installation
1) Download deepRAM
```bash
git clone https://github.com/MedChaabane/deepRAM.git

cd deepRAM
```

2) Install required packages 
```bash
pip3 install -r Prerequisites
```
3) Install deepRAM
```bash
python setup.py install
```
## Datasets
1. ChIP-seq datasets can be downloaded from:
http://tools.genes.toronto.edu/deepbind/nbtcode

2) CLIP-seq datasets can be downloaded from: https://github.com/xypan1232/iDeepS/tree/master/datasets/clip

We have provided two preprocessing scripts to change the format of the used datasets to a format compatible with deepRAM input data format (deepRAM input data format: sequence label. See [Example input data](https://github.com/MedChaabane/deepRAM/blob/master/datasets/example-input-data.gz)): 
- [preprocess_1.py](https://github.com/MedChaabane/deepRAM/blob/master/preprocess_1.py) can be used for [DeepBind](https://www.nature.com/articles/nbt.3300)-ENCODE-ChIP-seq-data-like format and, 
- [preprocess_2.py](https://github.com/MedChaabane/deepRAM/blob/master/preprocess_2.py) can be used for [iONMF](https://www.ncbi.nlm.nih.gov/pubmed/26787667)-CLIP-seq-data-like format.
## Example with CLIP-seq
#### preprocess CLIP-seq files (train and test) to match deepRAM data format: sequence label
```bash
python preprocess_2.py --CLIP_data datasets/CLIP-seq/1_PARCLIP_AGO1234_hg19/30000/training_sample_0/sequences.fa.gz --output CLIP_train.gz
```
```bash
python preprocess_2.py --CLIP_data datasets/CLIP-seq/1_PARCLIP_AGO1234_hg19/30000/test_sample_0/sequences.fa.gz --output CLIP_test.gz
```
#### train DeepBind architecture with CLIP_train.gz and evaluate performance on CLIP_test.gz
```
python deepRAM.py --train_data CLIP_train.gz --test_data CLIP_test.gz --data_type RNA --train True --evaluate_performance True --model_path DeepBind.pkl --out_file prediction.txt --Embedding False --Conv True --RNN False --conv_layers 1 
```
#### visualizating motifs and matching them with known motifs 
```
python deepRAM.py --test_data CLIP_test.gz --data_type RNA --predict_only True --model_path DeepBind.pkl --motif True --motif_dir motifs --tomtom_dir meme-5.0.3/src/tomtom --out_file prediction.txt --Embedding False --Conv True --RNN False --conv_layers 1
```
make sure to specify the directory of TOMTOM in --tomtom_dir argument
