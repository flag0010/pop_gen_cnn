import sys, os
from sklearn.externals import joblib
import numpy as np
import random
from sklearn.ensemble import ExtraTreesClassifier
from time import time
from operator import itemgetter
from sklearn.grid_search import GridSearchCV

trainingSetFileName, classifierPickleFileName = sys.argv[1:3]
statsToUse = sys.argv[3:]

XH = {}
with open(trainingSetFileName) as trainingSetFile:
    first = True
    for line in trainingSetFile:
        if first:
            header = line.strip().split()
            if header[0] != "classLabel":
                sys.exit("first column must be classLabel. AAAARRRRRRRRRGGGGGHHHH!!!\n")
            if "all" in statsToUse:
                indicesToKeep = range(1, len(header))
            else:
                indicesToKeep = []
                for i in range(1, len(header)):
                    if header[i] in statsToUse:
                        indicesToKeep.append(i)
            first = False
        else:
            line = line.strip().split()
            instance=[]
            isBad = False
            for i in indicesToKeep:
                if "nan" in line[i] or "inf" in line[i]:
                    isBad = True
                instance.append(float(line[i]))
            if not isBad:
                if not XH.has_key(line[0]):
                    XH[line[0]] = []
                XH[line[0]].append(instance)

minClassSize = min([len(XH[classLabel]) for classLabel in  XH.keys()])
X = []
y = []
for classLabel in sorted(XH.keys()):
    random.shuffle(XH[classLabel])
    for i in range(minClassSize):
        currVector = XH[classLabel][i]
        X.append(currVector)
        y.append(classLabel)
X = np.array(X)
sys.stderr.write("training set size after balancing: %s\n" %(len(y)))

# Utility function to report best scores
def report(grid_scores, n_top=3):
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
              score.mean_validation_score,
              np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")


sys.stderr.write("Checking accuracy when distinguishing among all %s classes\n" %(len(XH.keys())))

maxMaxFeatures=len(X[0])
param_grid_forest = {"max_depth": [3, 10, None],
              "max_features": [1, 3, int(maxMaxFeatures**0.5), maxMaxFeatures],
              "min_samples_split": [1, 3, 10],
              "min_samples_leaf": [1, 3, 10],
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

clf, mlType, paramGrid = ExtraTreesClassifier(n_estimators=100, random_state=0), "extraTreesClassifier", param_grid_forest

sys.stderr.write("Using %s\n" %(mlType))
grid_search = GridSearchCV(clf,param_grid=param_grid_forest,cv=10,n_jobs=10)
start = time()
grid_search.fit(X, y)
sys.stderr.write("GridSearchCV took %.2f seconds for %d candidate parameter settings.\n"
      % (time() - start, len(grid_search.grid_scores_)))
print "Results for %s" %(mlType)
report(grid_search.grid_scores_)
joblib.dump((statsToUse, header, grid_search), classifierPickleFileName)
