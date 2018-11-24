REIMPLEMENATION OF FILET USING DEEP NEURAL NET

The original FILET approach (https://github.com/kern-lab/FILET) used a decision tree ensemble to detect introgression.  It relied on a user specified feature vector based pop gen summary stats. The advantage of using a neural network is that the network itself extracts a feature vector during training, meaning it's not limited by the current the state of the art in pop gen summary stats.

Here we use a convolutional neural network architecture.  First I encoded the training data (bi-allelic SNPs encoded as binary matrices) as a sort of binary image. An example image can be found at the link below.

https://github.com/flag0010/FILET/blob/master/Figure_1.png

There are 34 lines, 20 from pop1 and 14 from pop2. Normally SNP matrices are displayed with individuals on rows and SNPs (segregating sites) on cols, but I transposed these SNP matrices for reasons I'll explain below.  Thus the image and SNP matrix are 1102x34. 

Next I had to deal with the fact that the training SNP matrices are of variable length.  In the training data supplied with FILET I found the longest had 1102 SNPs. I padded all training data to this length by adding zeros.  This is why the the bottom of the matrix image above is all black.

Because the SNP matrices are transposed a 1D convolution (col-wise) will stride across a single chromosome for all SNPs.    

Finally, as with the original version, the output classification includes 3 possibilities.  Migration from pop1 into pop2, the reverse, and no migration. I encoded these as follows: no-mig. = 0, 1->2 = 1, 2->1 = 2

To run you need python 2.7 and python 3, with numpy, keras and tensorflow. Ideally with the python 3 version set up on a gpu. 

Next download the 2 large gzipped files from the link below and add to this working directory: https://www.dropbox.com/sh/molccbxdm7zasbi/AABsVGfdZJvmF1PLkC3pY8Mva?dl=0

Then extract the data by running the following (coded for python 2.7):
```python extract.big.data.set.py```

This will produce a file called ```big_sim.npz```, which is the training and validation data.
(or alternatively you can skip the steps above and download `big_sim.npz` from this data repository: https://conservancy.umn.edu/handle/11299/198335)

Then in python3 fit a CNN:
```python train.neural.net.introgression.CNN.PYTHON3.py```

After 19 epochs I got an accuracy on validation data of approx 89%.  I stopped and output my trained model as ```big.data.89.2.acc.mod```.

This model can now be used on test data to get a final estimate of model quality (back in python 2.7!):
```extract.test.data.and.get.final.model.confusion.matrix.py```
