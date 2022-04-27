This is a toy model that is not in the paper. It is the simplest CNN we made (the whole thing fits in ~70 lines of code), and is useful for learning and for tweaking to see what happens. It also runs easily on a laptop in minutes, no GPU needed. So it's a good starting place. The model computes Tajima's D on made up sequence alignments, and then trains a 3 layer convolutional neural net to predict tajima's D.  
to run follow the steps below

1) install anaconda python 2.7 (https://www.anaconda.com/download/)

then install keras and tensorflow:

2) `pip install keras==2.0.0`
3) `conda install tensorflow==1.15`

then run

`python make.tajD.pred.py` 
