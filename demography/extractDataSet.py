import numpy as np
import sys, os, gzip, time
from msTools import *

msInFileName, demogParamFileName, maxSnps, maxDistBtwnSnps, outFileName = sys.argv[1:]
maxSnps = int(maxSnps)
maxDistBtwnSnps = float(maxDistBtwnSnps)


def combineHapsAndSnpDists(haplotypeMatrices, distsBetweenSnpVectors, maxSnps, maxDistBtwnSnps):
    sampleSize = len(haplotypeMatrices[0])
    inputMatrices = np.empty((len(haplotypeMatrices), sampleSize+1, maxSnps+1), dtype='uint8')
    assert len(haplotypeMatrices) == len(distsBetweenSnpVectors)
    for i in range(len(haplotypeMatrices)):
        unpaddedLen = len(distsBetweenSnpVectors[i])
        #print unpaddedLen, len(haplotypeMatrices[i][0]), len(haplotypeMatrices[i])
        assert unpaddedLen == len(haplotypeMatrices[i][0]) + 1
        padLen = maxSnps - unpaddedLen
        # distsBetweenSnpVectors includes distances to edges of chromosome
        k = 0
        while k < len(distsBetweenSnpVectors[i]):
            inputMatrices[i][0][k] = int(round(255*distsBetweenSnpVectors[i][k]/maxDistBtwnSnps))
            k += 1
        while k < maxSnps+1:
            inputMatrices[i][0][k] = 0
            k += 1
        for j in range(sampleSize):
            k = 0
            while k < len(haplotypeMatrices[i][j]):
                inputMatrices[i][j+1][k] = 255*int(haplotypeMatrices[i][j][k])
                k += 1
            while k < maxSnps+1:
                inputMatrices[i][j+1][k] = 0
                k += 1
    return inputMatrices

def getDistancesBetweenSnps(positionVectors):
    distVectors = []
    for i in range(len(positionVectors)):
        currVec = []
        prevPos = 0.0
        for j in range(len(positionVectors[i])):
            currVec.append(positionVectors[i][j]-prevPos)
            prevPos = positionVectors[i][j]
        currVec.append(1.0-prevPos)
        distVectors.append(currVec)
    return distVectors

def readDemogParams(demogParamPath):
    params = []
    first = True
    with open(demogParamPath) as demogParamFile:
        for line in demogParamFile:
            if first:
                first = False
            else:
                params.append([float(x) for x in line.strip().split()])
    return params

haplotypeMatrices, positionVectors, demogParams = [], [], []
haplotypeMatrices, positionVectors = msOutToHaplotypeMatrix(msInFileName)
demogParams = readDemogParams(demogParamFileName)
assert len(positionVectors) == len(haplotypeMatrices)
distsBetweenSnpVectors = getDistancesBetweenSnps(positionVectors)

X = np.array(combineHapsAndSnpDists(haplotypeMatrices, distsBetweenSnpVectors, maxSnps, maxDistBtwnSnps))
y = np.array(demogParams)
np.savez_compressed(outFileName, X=X, y=y)