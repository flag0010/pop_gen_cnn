import math, random

def selectPopSize(meanFoldChange, maxFoldChange):
    foldChange = 1+randBoundedExp(meanFoldChange-1, maxFoldChange-1)
    if random.random() <= 0.5:
        return foldChange
    else:
        return 1/foldChange

def selectGridOfTimesLogScale(gridSize, maxTime):
    gridTimes = []
    interval = 1.0/(gridSize-1)
    for i in range(0, gridSize, 1):
        gridTimes.append(math.exp(interval * i * math.log(maxTime+1))-1)
    return gridTimes

def randBoundedExp(mean, upper):
    val = upper+1
    while val > upper:
        val = random.expovariate(1.0/mean)
    return val

def writeParamFile(paramOutLines, paramFilePath):
    with open(paramFilePath, "w") as paramFile:
        for line in paramOutLines:
            paramFile.write(line)
