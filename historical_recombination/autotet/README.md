Performs the training and testing of a CNN to predict rho on an unphased (and un-genotype-resolved) autotetraploid. 
Requires installation of keras, tensorflow, numpy, scikit-learn, matplotlib, and a set of functions called common (https://github.com/flag0010/python_common/blob/master/common.py).

All python code runs on python2.7, EXCEPT the model training.  That was done in the cloud on python 3.6!  Sorry!

To run, first run the shell script to produce coalescent simulations (requires ms)

`./autotet.run.ms.sh`

then gzip the results

`gzip autotet.all.LD.sims.txt`

Next run the python2.7 script

`python conv.training.data.to.autotet.py`

This will output a large compressed file called "autotet.ld.data.npz". Use this file as training data for the network.

To train the network, first find a computer with a GPU and get keras and tensorflow installed.  This was done in python3.6
 
 `python3 train.autotet.mergenet.PYTHON3.py`
 
 The model trains fast.  It will print the mean of the training data.  Save that value (should be about 4.8).  
 After training the model it prints a json string that specifies the model structure.  You'll want to save that too.  Also the script saves the weights from the CNN in a file called '3rd.autotet.mergnet.weights'
 
 Next generate some new independent test data with ms.  Run
 
 `./run.ms.test.sh`
 
 and gzip the output as above (or skip, i've provided a file `autotet.test.data.LD.sims.txt.gz`)
 
 Finally test the CNN on this independent test data, (back in python2.7)
 
 `python final.autotet.test.results.py`
 
 This will print the R^2 and RMSE and a plot
