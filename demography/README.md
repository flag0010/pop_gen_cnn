This directory contains the files necessary to perform the demographic inference for a 3-epoch model of population size changes described in our manuscript. We assume that you have already generated simulations using ms or a simulator that produces similarly-formatted output, and have recorded for each simulation the true demographic parameters as specified below. We now describe each file in this directory in turn.

1) The file extractDataSet.py is used to format data for training as follows:

`python extractDataSet.py $msFilePath $demogParamPath $maxSegsites $maxDistance $outFilePath`

The arguments are as follows:

msFilePath: the path to an ms-formatted simulation output file

demogParamPath: the path to a tab-delimited file that lists on each line the demographic parameters corresponding to the  simulation output. Thus, the ith line gives the parameters corresponding to the ith rep from demogParamPath. The order of these parameters are: present-day population size, time of recent size change, middle-epoch population size, time of ancient size changes, ancient population size. The population sizes are in numbers of individuals, while the times are in units of 4N_0 * generations, where N_0 is the present-day population size.

maxSegsites: an integer specifying the maximum number of segregating sites that will be included in any data set that the CNN will be applied to (including training, testing, validation, or real data sets).

maxDistance: a number between zero and 1 specifying the maximum distances between any two SNPs that will be encountered across all data sets. A value of one corresponds to the full length of the simulated chromosome.

outFilePath: the path to an npz-file (https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.savez.html) where the data will be stored until use by demogConvRegMerge.py (see below).


2) The file demogConvRegMerge.py reads in files produced by extractDataSet.py and performs training and testing as described in the manuscript.

`python demogConvRegMerge.py $convDim $convSize $poolSize $logTransformY $intAllele $sortRows $useDropout $lossThreshold $inDir $weightFileName $modelFileName $resultFileName`

The arguments are:

convDim: Specifies whether 1- or 2-dimensional convolutions (and pooling steps) are used. Must be "1d" or "2d".

convSize: An integer specifying the convolutional kernel size. Note that this script uses a square matrix for all 2D convolutions so both 1D and 2D kernel sizes are specified by a single value.

poolSize: An integer specifying the kernel size for pooling steps. Similar to convSize. (Was always set to 2 for the analyses presented in the manuscript.)

logTransformY: Specifies whether response variables are transformed (when set to "True") or not (when set to "False")

intAllele: Specifies whether alleles are encoded as 0/255 ("True") or -1/1 ("False").

sortRows: Specifies whether the haplotypes in the input matrix will be sorted by similarity.

useDropout: Specifies whether dropout steps are included ("True" or "False").

lossThreshold: Specifies a maximum threshold of the loss value (computed using the keras mean_squared_error function) that a model must achieve in order for its results to be written.

inDir: Input directory of .npz files produced by extractDataSet.py. Test and validation sets will both contain 10000 examples, so the file(s) included in this directory must contain a total of >20000 examples.

weightFileName: A path to an output file where the weights of the best-scoring neural network are to be written. This will be removed prior to termination of lossThreshold is not met. Together with the model file (below), this file can be used to load a network back into keras for future prediction.

modelFileName: A path to an output file that stores the architecture of the neural network. This will only be written if the lossThreshold is met. Together with the weights file (above), this file can be used to load a network back into keras for future prediction.

resultFileName: Path to a file where the true and predicted response values will be written for each test example.

3) msTools.py contains a function used by extractDataSet.py
