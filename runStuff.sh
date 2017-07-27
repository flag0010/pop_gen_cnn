#!/bin/bash

n1=20
n2=14
windowSize=10000
#mkdir -p trainingSimsStats trainingSets classifier results
#gunzip trainingSims/*.msOut.gz
#for inFile in `ls trainingSims/ | grep .msOut` ; do cat trainingSims/$inFile | ./twoPopnStats_forML $n1 $n2 | python normalizeTwoPopnStats.py None $windowSize > trainingSimsStats/$inFile; done
#gzip trainingSims/*.msOut
#python buildThreeClassTrainingSet.py trainingSimsStats/ trainingSets/threeClass.fvec
python trainFiletClassifier.py trainingSets/threeClass.fvec classifier/threeClass.p pi1 hetVar1 ss1 private1 tajd1 ZnS1 pi2 hetVar2 ss2 private2 tajd2 ZnS2 Fst snn dxy_mean dxy_min gmin zx dd1 dd2 ddRank1 ddRank2
./pgStatsBedSubpop_forML dataToClassify/pop1.fa dataToClassify/pop2.fa dataToClassify/anc.fa dataToClassify/coords.bed 0.5 | python normalizePgStats.py $windowSize > featureVectorsToClassify/test.ss
python classifyChromosome.py classifier/threeClass.p featureVectorsToClassify/ results/
