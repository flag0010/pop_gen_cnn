import sys

winSize=int(sys.argv[1])
missingDataThreshold=0.25

line = sys.stdin.readline()
first = True
colsToDiv, colsToDivTwice = [], []
while line:
    line = line.strip().split()
    if first:
        for i in range(len(line)):
            if line[i] in "pi1,private1,ss1,thetaH1,H1,pi2,private2,ss2,thetaH2,H2,dxy_mean,dxy_min,gmin".split(","):
                colsToDiv.append(i)
            elif line[i] in "hetVar1,hetVar2".split(","):
                colsToDivTwice.append(i)
            elif line[i] == "numSites":
                numSitesCol = i
        first = False
        print "\t".join(line)
    else:
        outLine = []
        numSites = float(line[numSitesCol])
        if numSites/winSize > missingDataThreshold:
            for i in range(len(line)):
                if "inf" in line[i] or "nan" in line[i]:
                    outLine.append(line[i])
                elif i in colsToDiv:
                    #outLine.append("%s,%f" %(line[i], float(line[i])/numSites))
                    outLine.append("%f" %(float(line[i])/numSites))
                elif i in colsToDivTwice:
                    outLine.append("%f" %(float(line[i])/numSites**2))
                else:
                    outLine.append(line[i])
            print "\t".join(outLine)
    line = sys.stdin.readline()
