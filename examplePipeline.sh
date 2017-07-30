#!/bin/bash

# This script describes each step of the pipeline for training FILET and running it on an actual
# real dataset. Example commands are given, and the input files required to get the pipeline
# started. Thus, if you clone the repository and cannot run through the whole pipeline then
# something is wrong and you may be missing one of the dependencies (see README file).

# Some steps in this pipeline may take a fair amount of time (e.g. several hours), so you may wish
# to run it overnight. On much larger datasets you may wish to use a compute cluster. and split
# up the computation. Anyway, here we go:

# Step 0: Simulate your training data (not shown here). This can be done with msmove, written by
# Daniel Garrigan and Anthony Geneva (https://github.com/geneva/msmove).

n1=20
n2=14
windowSize=10000
# Step 1: Make directories for intermediate files and final classification output
mkdir -p trainingSimsStats trainingSets classifier results

# Step 2: Calculate summary statistics for each input file in our training simulations.
# This step may take several hours for very large training sets. The input for twoPopnStats_forML
# (which is piped in via stdin) is ms-style simulation output. n1 and n2 give the sample sizes of
# the two population samples in these simulated data. The output file contains one tab-separated
# line for each simulation, this line gives the values of various population genetic summary
# statistics calculated within the simulated region.
# Note: the masking step described in our paper and used for  simulans/sechellia analysis is
# omitted here for simplicity.
#
# normalizeTwoPopnStats.py simply takes the output from twoPopnStats_forML and for statistics that
# scale with window size (e.g. pi, the number of segregating sites, the number of private alleles)
# the value of the statistic is divided by the number of sites in the window. The script takes
# two arguments: a file listing the fraction of missing data in each simulation replicate (just
# use the word None if no masking was performed), and the window size (10 kb in this example).
# The output is again a tab-separated list of values for each simulation.
gunzip trainingSims/*.msOut.gz # first we have to unzip the simulation output
for inFile in `ls trainingSims/ | grep .msOut` ; do cat trainingSims/$inFile | ./twoPopnStats_forML $n1 $n2 | python normalizeTwoPopnStats.py None $windowSize > trainingSimsStats/$inFile; done

# Step 3: Aggregate our training data into one labeled data matrix to be used for training.
# The first argument is the directory which contains all input files, and the second argument
# is the name of the output file. The input file names must be noMig.msOut (no migration, assigned
# label "0"), mig12.msOut (training examples for migration from population 2 into population 1,
# assigned label "1"), and mig21.msOut (training examples for migration from population 1 into
# population 2, assigned label "2"). The output file contains all of these feature vectors but
# with an additional column: the appropriate class lable.
python buildThreeClassTrainingSet.py trainingSimsStats/ trainingSets/threeClass.fvec

# Step 4: Train the classifier. The first two arguments are the paths to the training data and
# classifier pickle (see below) to be written, respectively. The remaining arguments are the names 
# of all statistics to be included in the classifier.
#
# Note: this step may take several hours for ~10000 training examples.
#
# Additional note: for this example I have omitted several statistics that are sensitive to the
# length of the sequence but for whcih this dependency is not linear (haplotype counts and IBS
# tract length statistics). It would be better to include these only if the training data have
# been masked in the same manner as the data the user wishes to classify.
python trainFiletClassifier.py trainingSets/threeClass.fvec classifier/threeClass.p pi1 hetVar1 ss1 private1 tajd1 ZnS1 pi2 hetVar2 ss2 private2 tajd2 ZnS2 Fst snn dxy_mean dxy_min gmin zx dd1 dd2 ddRank1 ddRank2

# Stpe 5: Meanwhile, we need some data to classify. The program pgStatsBedSubpop_forML does this
# for us. The arguments, in order, are the phased sequence data for population/species 1 in our 
# region or chromosome of interest (in fasta format), the same for population 2, a fasta file
# with one entry representing the inferred ancestral state for each position in our region
# of interest (can just use data from the reference genome if not using output information -- in
# this case be sure not to incorporate into your classifier either  Fay and Wu's thetaH or H
# statistics for either population), a BED file specifying the coordinates of the windows within
# our region of interest (which can be an entire chromosome/scaffold or a piece of one; all 
# coordinates here are zero-based and relative to the beginning of the region included in the
# fasta files--i.e. if our fasta files contain sequence for chr1:1001-2000 then the BED entry
# "chr1\t0\t100" corresponds to chr1:1001-1100), and the maximum fraction of missing data per
# site (e.g. in the example above, if >=50% of phased haplotypes are 'N' at a given site, then
# then entire site is ignored during the calculation of summary statistics.
# Note that for the fasta-formatted input files, each sequence must be the same length (i.e. the
# same chromosome/scaffold/region) though different files may have different numbers of sequences
# (i.e. pop1.fa and pop2.fa will have one entry for each of their haploid sample sizes n1 and n2,
# and anc.fa should have only one entry). The BED format is described here:
# https://genome.ucsc.edu/FAQ/FAQformat.html#format1
#
# The normalizePgStats.py script divides all measures of variation that scale with the number of
# sites (e.g. pi, the number of segregating sites, the number of private alleles) by the number
# of sites not filtered from the window (the numSites field from pgStatsBedSubpop_forML's output).
# Windows missing at least 75% of all sites are filtered out by this script.
#
# The output of this step is a tab-delimited file giving the coordinates for each window, the
# number of unmasked sites in the window and the values of various summary statistics.
#
# Note: the example files included in this script are from chr4 in the simulans-sechellia dataset.
# These data may not be of adequate quality for analysis, but are useful for testing the pipeline.
./pgStatsBedSubpop_forML dataToClassify/pop1.fa dataToClassify/pop2.fa dataToClassify/anc.fa dataToClassify/coords.bed 0.5 | python normalizePgStats.py $windowSize > featureVectorsToClassify/test.ss

# Step 6: Use classifier to produce classifications on the real data set.
# The first argument is the path to the classifier pickle, the second is the data we wish to
# classify (described above), and the third is the path where we wish to write our results.
python classifyChromosome.py classifier/threeClass.p featureVectorsToClassify/ results/

# Note about pickles: this pipeline saves our classifier as a "pickle" using scikit-learn's joblib
# library. Beware that pickles can contain executable python code, so if your pickles are tampered
# with there is the potential for execution of malicious code. Make sure your pickles cannot be
# edited by anyone else! For more information on pickles and joblib, see
# http://scikit-learn.org/stable/modules/model_persistence.html
#
# The joblib pickles can consist of many files, so it is best to save each pick;e in their own
# directory. Note that the formatting of pickles is system-dependent, so you may have to re-train
# the classifier on each machine where you wish to use it.
