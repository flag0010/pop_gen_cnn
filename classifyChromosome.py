import sys, os
import numpy as np
import random
from sklearn.externals import joblib
from sklearn.ensemble import ExtraTreesClassifier

classifierPickleFileName, armDir, outDir = sys.argv[1:]
desiredFpRate = 0.05

statsToUse, header, grid_search = joblib.load(classifierPickleFileName)

statIndices = []
for i in range(1, len(header)):
    if header[i] in statsToUse or "all" in statsToUse:
        statIndices.append(i)

def probsToPreds(probs, desiredFpRate):
    preds = []
    for prob in probs:
        if prob[0] >= desiredFpRate:
            preds.append("0")
        elif prob[1] > prob[2]:
            preds.append("1")
        elif prob[2] > prob[1]:
            preds.append("2")
    return np.array(preds)

def writePreds(predictions, probs, coords, outFileName):
    with open(outFileName, "w") as outFile:
        for i in range(len(coords)):
            outLine = coords[i] + [predictions[i]] + list(probs[i])
            outFile.write("\t".join([str(x) for x in outLine]) + "\n")

for testSetFileName in os.listdir(armDir):
    sys.stderr.write("classifying " + testSetFileName + "\n")
    with open(armDir + "/" + testSetFileName) as testSetFile:
        currTestData = testSetFile.readlines()
        testHeader = currTestData[0].strip().split()
        assert testHeader[4:] == header[1:]
        currTestData.pop(0)#remove the header from the test data file

    testX = []
    coords = []
    for i in range(len(currTestData)):
        testDatum = currTestData[i].strip().split()
        currVector = []
        isBad = False
        for j in range(1, len(header)):
            if j in statIndices:
                if "nan" in testDatum[j+3] or "inf" in testDatum[j+3]:
                    isBad = True
                currVector.append(float(testDatum[j+3]))
        assert len(currVector) == len(statIndices)
        if not isBad:
            testX.append(currVector)
            coords.append(testDatum[:4])
    if len(testX) == 0:
        sys.stderr.write("\tskipping: not enough data!\n")
    else:
        testX = np.array(testX)
        probs = grid_search.predict_proba(testX)
        predictions = probsToPreds(probs, desiredFpRate)
        outFileName = outDir + "/%s.preds" %(testSetFileName)
        writePreds(predictions, probs, coords, outFileName)
        sys.stderr.write("\tdone.\n")
