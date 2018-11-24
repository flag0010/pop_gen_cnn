# Population genetic inference with Convolutional Neural Networks

This directory contains code for training and testing the neural networks described in this paper: _The Unreasonable Effectiveness of Convolutional Neural Networks in Population Genetic Inference_
by Lex Flagel, Yaniv Brandvain, and Daniel Schrider. https://doi.org/10.1101/336073

Each folder contains the code and README files for a separate problem discussed in the paper.  All code is presented as is and runs in python 2.7 unless stated otherswise. Code is specific to the problems described in the paper, but can be modified to address other problems in population genetics.  The folders called tajimas_D and data_prep_tricks were not used in the paper above, but were used in the development of ideas presented in the paper.  If you are looking for a good starting place to play with convolutional neural networks, check out the tajimas_D directory.  It's a very simple model that will run in a minute or two on a laptop. The rest of the code in this repo will require a computer with significant amounts of RAM (~64 GB) and Tensorflow/Keras combined with a GPU (we used on an Nvidia K80).

