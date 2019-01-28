# DeepRAM
![pipeline](https://github.com/MedChaabane/deepRAM/blob/master/CSU-Ram.jpg)

we proposed an end-to-end deep learning toolkit, deepRAM, for predicting protein binding sites and motifs. It helps users run experiments using many state-of-the-art methods and addresses the challenge of selecting model parameters in deep learning models using a fully automatic model selection strategy. This helps avoid hand-tuning and thus removes any bias in running experiments, making it user friendly without losing its flexibility. While it was designed with ChIP-seq and CLIP-seq data in mind, it can be used for any DNA/RNA sequence binary classification problem.

deepRAM allows users the flexibility to choose a deep learning model by selecting its different components:  input sequence representation (one-hot or k-mer embedding), whether to use a CNN and how many layers, and whether to use an RNN, and the number of layers and their type. For CNNs the user can choose to use dilated convolution as well.
 <br><br>
# Dependency <br>
<a href=https://github.com/fchollet/keras/>keras 1.2.0 library</a> <br>
<a href=https://github.com/scikit-learn/scikit-learn>sklearn</a> <br>
<a href=http://www.h5py.org/>h5py</a>, install it using "pip install h5py" <br>
python 2.7 <br>
