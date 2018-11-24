This directory trains a neural network to predict selective sweeps from sequence data.  It detects 5 classes:  a recent hard sweep (i.e. fixation of a de novo beneficial mutation), a recent soft sweep (i.e. fixation of a beneficial but previously neutral segregating polymorphism), a region linked to a nearby hard sweep, a region linked to a nearby soft sweep, and a neutrally evolving region.

To run first download this very large npz file and place it in this directory: https://drive.google.com/open?id=165CTU51MroaCHscZzZBSOPak12kju7Ut
The npz file is nearly 14 Gb. Unlike the other models in this repo, the training data for this one is too large to read into memory at one time. So we read it piecemeal and train the CNN seqeuntially.

To train the CNN first run: `python3 train.model.PYTHON3.py`

This will sequentially train the CNN for 3 epochs and at the end run on the test data producing the final confusion matrix found in the paper. There are 5 classes of sweeps, and in the model they are encoded as follows:
neutral=0, hard=1, hard-near=2, soft=3, soft-near=4

For now this CNN can only be run using our previously produced training file used in the manuscript. Later I will update with discoal code used for producing your own training data.   
