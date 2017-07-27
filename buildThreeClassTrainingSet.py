import sys, os

def readAndLabelStats(statFileName, label):
    examples = []
    with open(statFileName) as statFile:
        first = True
        for line in statFile:
            if first:
                header = "classLabel\t"+line
                first = False
            else:
                examples.append(str(label)+"\t"+line)
    return header, examples

statDir, trainingSetFileName = sys.argv[1:]
headerH = {}
header, migrantTrainingSet = readAndLabelStats("%s/mig12.msOut" %(statDir), 1)
headerH[header] = 1
header, revMigrantTrainingSet = readAndLabelStats("%s/mig21.msOut" %(statDir), 2)
headerH[header] = 1
header, nonMigrantTrainingSet = readAndLabelStats("%s/noMig.msOut" %(statDir), 0)
headerH[header] = 1
if len(headerH) != 1:
    sys.exit("Not all headers are identical. AAAARRRGGGGHHHHHH!!!\n")
with open(trainingSetFileName, "w") as outFile:
    outFile.write(header)
    for line in migrantTrainingSet:
        outFile.write(line)
    for line in revMigrantTrainingSet:
        outFile.write(line)
    for line in nonMigrantTrainingSet:
        outFile.write(line)
