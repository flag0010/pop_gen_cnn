This directory builds a CNN to predict theta from coalescent sims, and determine the impact of various data reoganization techniques.

To run, first download ms and install executable in this directory. https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13

Then you can create a shell script by running:
`python make.sim.sh.py > sim.theta.sh`
 
 Then run the shell script:
 `./sim.theta.sh > theta.sims.txt`

Compress output with gzip:
`gzip theta.sims.txt`

Next extract a npz file by running: 
`python extract.data.set.for.training.py`

This creates the file theta_sim.npz, which is used to train the neural network.  Training is done in python 3:

`python3 train.nn.PYTHON3.py`

Then you can plot the figure from the paper using:
`python final.plot.py`
