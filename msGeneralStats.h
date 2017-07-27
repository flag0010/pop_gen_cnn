/*msGeneralStats.h -- contains general analysis of Hudson format output */

#ifndef MSGS_INC
#define MSGS_INC
#include "stdio.h"


void imatrix_free(int **m, int nsam);
void cmatrix_free(char **m, int nsam);
int **imatrix(int nsam,int len);

//code from achaz
double * compute_HarmonicSums( int n );
double beta_FU1995( int i, double *HarmonicSums, int n  );
double sigma_ii_FU1995( int i, double *HarmonicSums, int n );
double sigma_ij_FU1995( int i, int j, double *HarmonicSums, int n );

int frequency( char allele,int site,int nsam,  char **list);
double nucdiv(int, int, char **);
double maxFDA( int nsam, int segsites, char **list);
double thetah(int, int, char **);
int segSites(int, int, char **);
void getSiteFreqSpec(int segsites, int nsam, char**list, int nSites, int *derived_counts);
int sortcmp (int *n1, int *n2);
int getHaplotypeFreqSpec(int segsites, int nsam, char **list, int *haplotype_counts);
double meanEHH(int segsites, double *posit, double delta, int nsam, char **list);
double meanREHH(int segsites, double *posit, double delta, int nsam, char **list);
void petrovHStatsWindow(int segsites, int nwins, double *posit, double *winsH12, double *winsH1, double *winsH2, int nsam, char **list);
double petrovH1(int *haplotype_counts, int nsam);
double petrovH2(int *haplotype_counts, int nsam);
double petrovH12(int *haplotype_counts, int nsam);
int nHaplotypes(int segsites, int nsam, char **list);
int biggerlist(int nsam,unsigned nmax,char **list );
int sampleSizeSite(int site, int nsam, char **list);
double ZnS(int segsites, int nsam, char** list);
double dij(int i, int j, int nsam, char** list);
double omega(int left, int segsites,int nsam, char** list);
double omegaWithTable(int left, int segsites, int nsam, char** list, double** dijTable);
double omegaMax(int segsites,int nsam, char** list);
double omegaCenter(int siteIdx , int segsites,int nsam, char** list);
double sigmaIJ(int i, int j, int n);
double sigmaII(int i, int n);
double sigmaBeta(int i,int n);
double sigmaAlpha(int n);
double achazNeutTestExponentWeights(int nsam, int segsites, char **list, int exponent1, int exponent2, double *harmonicSums);
double achazThetaExponentWeights(int nsam, int segsites, char **list, int exponent);
double achazThetaParabolicWeights(int nsam, int segsites, char **list, int exponent,double center);

double achazThetaHPi(int nsam, int segsites, char **list);
double achazTajimasDExtreme(int nsam, int segsites, char **list);
void nucdivWindow( int nwins, double *posit, double *output, int nsam, int segsites, char **list);
void tajdWindow(int nwins, double *posit, double *output, int nsam, int segsites, char **list);
void getSiteFreqSpecWindow(int segsites, int nsam, char**list, int nSites, int *derived_counts, double *pos, double low, double high);

double countAlleleDiffsForSnpPair(int nsam, int i, int j, char **list);
double sStarSnpDist(int nsam, int i, int j, int physLen, double *posit, char **list);
double sStar(int nsam, int segsites, int physLen, double *posit, char **list);

//subpopn declares
int frequencySub(char allele, int site, int startAllele, int stopAllele, char **list);
double nucdivSub( int nsam, int segsites, int startAllele, int stopAllele, char **list);
void privateSegSitesInTwoPopns(int segsites, int nsam, int stopAllele1, int *private1, int *private2, char **list);
double thetahSub( int nsam, int segsites, int startAllele, int stopAllele, char **list);
int nHaplotypesSub(int segsites, int nsam, int startAllele, int stopAllele, char **list);
double ZnSSub(int segsites, int nsam, int startAllele, int stopAllele, char** list);
double fst2Subs(int segsites, int nsam, int start1, int stop1, int start2, int stop2, char **list);
int segSitesSub(int segsites, int nsam, int startAllele, int stopAllele, char **list);
double Snn(int segsites,int nsam, int n1, int n2, char **list);
double xij_Snn(int segsites,int nsam, int seqIndex1, int n1, int n2, char **list);
double seqDist_Snn(int segsites, int index1,  int index2, char **list);
double meanRefDist(int segsites, int nsam, char **list);
double Dxy(int segsites,int nsam, int n1, int n2, char **list);
double Dxy_min(int segsites,int nsam, int n1, int n2, char **list);
double Dxy_mean(int segsites,int nsam, int n1, int n2, char **list);
void pairwiseDistances(int segsites,int nsam, double *dists, char **list);
double pairwiseDistRankAmongSampleRange(int segsites, int diffs, int firstSample, int numSamples, double *hetVar, char **list);
void Dxy_vector(int segsites,int nsam, int n1, int n2, double *v,char **list);

void nucdivSubWindow( int nwins, double *posit, double *output, int nsam, int segsites,int startAllele, int stopAllele, char **list);
void fst2SubsWindow(int nwins, double *posit, double *output,int segsites, int nsam, int start1, int stop1, int start2, int stop2, char **list);

double *pairwiseIBSVec1Popn(int segsites, int n1, int *vecLen, double *posit, char **list);
double *hetVec1Popn(int segsites, int n1, int physLen, int *vecLen, double *posit, char **list);
void reorderListIntoClusters(int *membership1, int g1Size, int *membership2, int g2Size, char **list);
int breakClusterAssignmentTie(double **hetMatrix, int targIndex, int *membership1, int m1Index, int *membership2, int m2Index);
int assignClusters(double *hetVec, int n1, int *g1Size, int *membership1, int *g2Size, int *membership2);
void clusterSeqsFromUnsortedHetVec(double *hetVec, int n1, int *g1Size, int *g2Size, char **list);
void statVecMoments(double *vec, int vecLen, double *mean, double *var, double *skew, double *kurt);
void statVecMinMedMax(double *vec, int vecLen, double *minVal, double *medVal, double *maxVal);
void pairwiseIBSHist2Popn(int segsites,int nsam, int n1, double *hist, int bins, double *posit, char **list);
double pairwiseIBSMax2Popn(int segsites,int nsam, int n1, double *posit, char **list);
double pairwiseIBSMeanWithin(int segsites,int first, int last, double *posit, char **list);


//twoSite
void sampleConfig(int i, int j, int nsam, char** list, int *config);
void printPairwiseSampleConfigs(int segSites, int nsam, char **list,double *posit, int nsites);
void printPairwiseSampleConfigs2Popn(int segsites, int nsam, int popnSize1, char **list, double *posit, int nsites);
void sampleConfig2Popn(int i, int j, int nsam, int popnSize1, char** list, int *config);

//From Hudson
double tajd(int nsam, int segsites, double sumk);
double a1f(int);
double a2f(int);
double b1f(int);
double b2f(int);
double c1f(double, double);
double c2f(int, double, double, double);
double e1f(double, double);
double e2f(double, double, double);

//utility
int cmp_doubles(const void *x, const void *y);

#endif
