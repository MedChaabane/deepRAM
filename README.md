# deepRAM

   <img src="https://github.com/MedChaabane/deepRAM/blob/master/CSU-Ram.jpg" width="150">

deepRAM is an end-to-end deep learning toolkit for predicting protein binding sites and motifs. It helps users run experiments using many state-of-the-art deep learning methods and addresses the challenge of selecting model parameters in deep learning models using a fully automatic model selection strategy. This helps avoid hand-tuning and thus removes any bias in running experiments, making it user friendly without losing its flexibility. While it was designed with ChIP-seq and CLIP-seq data in mind, it can be used for any DNA/RNA sequence binary classification problem.

deepRAM allows users the flexibility to choose a deep learning model by selecting its different components:  input sequence representation (one-hot or k-mer embedding), whether to use a CNN and how many layers, and whether to use an RNN, and the number of layers and their type. For CNNs, the user can choose to use dilated convolution as well.
 <br><br>
## Dependency <br>
We recommend to use [Anaconda 3](https://www.anaconda.com/download/) platform. 
python 3.6 <br>
<a href=https://pytorch.org/>PyTorch 1.0 library </a> (Deep learning library) <br>
<a href=https://github.com/scikit-learn/scikit-learn>sklearn</a> (Machine learning library)<br>
<a href=https://anaconda.org/anaconda/gensim>gensim</a> (library used to train word2vec algorithm) <br>
<a href=https://anaconda.org/anaconda/numpy>numpy</a> <br>

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

sequence specificities prediction of DNA- and RNA-binding proteins using deep learning approach

optional arguments:
  -h, --help            show this help message and exit
  --train_data TRAIN_DATA
                        path for training data with format: sequence label
  --test_data TEST_DATA
                        path for test data containing test sequences with or
                        without label
  --data_type DATA_TYPE
                        type of data: DNA or RNA
  --train TRAIN         use this option for automatic calibration, training
                        model using train_data and predict labels for
                        test_data
  --predict_only PREDICT_ONLY
                        use this option to load pretrained model (found in
                        model_path) and use it to predict test sequences
                        (train will be set to False).
  --evaluate_performance EVALUATE_PERFORMANCE
                        use this option to calculate AUC on test_data. If
                        True, test_data should be format: sequence label
  --models_dir MODELS_DIR
                        The directory to save the trained models for future
                        prediction including best hyperparameters and
                        embedding model
  --model_path MODEL_PATH
                        If train is set to True, This path will be used to
                        save your best model. If train is set to False, this
                        path should have the model that you want to use for
                        prediction
  --motif MOTIF         use this option to generate motif logos
  --motif_dir MOTIF_DIR
                        directory to save motifs logos
  --tomtom_dir TOMTOM_DIR
                        directory of TOMTOM, i.e:meme-5.0.3/src/tomtom
  --out_file OUT_FILE   The output file used to store the prediction
                        probability of testing data
  --Embedding EMBEDDING
                        Use embedding layer: True or False
  --Conv CONV           Use conv layer: True or False
  --RNN RNN             Use RNN layer: True or False
  --RNN_type RNN_TYPE   RNN type: LSTM or GRU or BiLSTM or BiGRU
  --kmer_len KMER_LEN   length of kmer used for embedding layer, default=3
  --stride STRIDE       stride used for embedding layer, default=1
  --word2vec_train WORD2VEC_TRAIN
                        set it to False if you have already trained word2vec
                        model. If you set it to False, you need to specify the
                        path for word2vec model in word2vec_model argument.
  --word2vec_model WORD2VEC_MODEL
                        If word2vec_train is set to True, This path will be
                        used to save your word2vec model. If word2vec_train is
                        set to False, this path should have the word2vec model
                        that you want to use for embedding layer
  --conv_layers CONV_LAYERS
                        number of convolutional modules
  --dilation DILATION   the spacing between kernel elements for convolutional
                        modules (except the first convolutional module)
  --RNN_layers RNN_LAYERS
                        number of RNN layers

```

## Motifs identification and visualization

You need to install <a href=http://weblogo.berkeley.edu/> WebLogo </a> and TOMTOM in <a href=http://meme-suite.org> MEME Suite </a> to match identifyed motifs with known motifs of Transcription Factors and RBPs. Read documentations about installation and usage.

## Installation
1 Download deepRAM
```bash
git clone https://github.com/MedChaabane/deepRAM.git

cd deepRAM
```

2 Install required packages 
```bash
pip3 install -r Prerequisites.txt
```
3 Install deepRAM
```bash
python setup.py install
```
