//pgSummaryStats.h

#ifndef SUMSTAT_INC
#define SUMSTAT_INC

#include "miscCode/stringWrap.h"
#include "miscCode/sequenceMatrix.h"
#include "miscCode/vector.h"
#include <stdlib.h>
#include <stdio.h>


int numColumnsNotNedOutFromTo(struct sequenceMatrix *aSeqMat, int start, int end);
int numColumnsNotNedOutFromToBothMats(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start, int end);
double rEHHatFocalSnp(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *ancMat, int start, int stop, int focalSnpPos);
int getHaplotypeFreqSpec(struct sequenceMatrix* aSeqMat, int beg, int end, int *haplotype_counts);
double petrovH12(int *haplotype_counts, int nsam);
double petrovH1(int *haplotype_counts, int nsam);
double petrovH2(int *haplotype_counts, int nsam);
int nHaplotypes(struct sequenceMatrix* aSeqMat, int beg, int end);
int fixedDiffsPre(stringWrap** set, int size, int samplesize, stringWrap* outgdata);
double missingDataPre(stringWrap** array, int size, int samplesize);
double missingDataFromTo(struct sequenceMatrix* aSeqMat, int start, int end);
int frequencyFD( char allele,int site,int nsam,  char **list);
int ingroupSegSites(int segsites, int nsam, char **list);
int fixedDiffsAllele1(int segsites, int nsam, char **list);
int frequencyFDAllele1( char allele,int site,int nsam,  char **list);
double nucdivPre(stringWrap** set, stringWrap** array, int size);
double nucdiv(struct  sequenceMatrix *aSeqMat);
double nucdivFromTo( struct sequenceMatrix *aSeqMat, int start, int end);
double maxFDA(struct sequenceMatrix *aSeqMat);
double maxFDAFromTo(struct sequenceMatrix *aSeqMat, int start, int end);
int segSiteCountPre(stringWrap** set, int size);
int segSiteCountFromTo(struct sequenceMatrix *aSeqMat, int start, int end);
int segSiteCountSubpopFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start, int end);
void segSiteLocationsBiallelicFromTo(struct sequenceMatrix *aSeqMat, vector *segs, int start, int end);
void segSiteBiallelicLocations(struct sequenceMatrix *aSeqMat, vector *segs);
double thetaHPre(stringWrap** set, stringWrap** array,stringWrap** ancSet , int size);
double thetaHFromTo(struct sequenceMatrix *aSeqMat,struct sequenceMatrix *ancMat, int start, int end);
double thetaH(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *ancMat);
int segSiteCount(struct sequenceMatrix *aSeqMat);
void privateSegSitesInTwoPopnsFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int *private1, int *private2, int start, int end);
double tajD(struct sequenceMatrix *aSeqMat);
double tajDPre(stringWrap** set, stringWrap** array, int size, int samplesize);
double tajDFromTo(struct sequenceMatrix *aSeqMat,int start, int end);
double thetaWPre(stringWrap** set, int size, int samplesize);
double thetaWFromTo(struct sequenceMatrix *aSeqMat, int start, int end);
double thetaW(struct sequenceMatrix *aSeqMat);
double frequency( stringWrap *array, char aChar);
int count( stringWrap *array, char aChar);

void segSiteLocations(struct sequenceMatrix *aSeqMat, vector *segs);
void segSiteLocationsFromTo(struct sequenceMatrix *aSeqMat, vector *segs, int start, int end);
void segSiteLocationsSubpopFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, vector *segs, int start, int end);
char majorAlleleSite(struct sequenceMatrix *aSeqMat, int index);
int min_rec (struct sequenceMatrix *aSeqMat, int x);
double freqAlleleSite(struct sequenceMatrix *aSeqMat, char aChar, int index);
double rSquared(struct sequenceMatrix *aSeqMat, int indexA, int indexB);
double rSquaredOmega(struct sequenceMatrix *aSeqMat, int indexA, int indexB);
void rSquaredCounts(struct sequenceMatrix *aSeqMat, int indexA, int indexB);
double jointHeterozygosity(struct sequenceMatrix *aSeqMat, int indexA, int indexB);
double rSquared_2chroms(struct sequenceMatrix *aSeqMat,struct sequenceMatrix *bSeqMat, int indexA, int indexB);
void rSquaredCounts_2chroms(struct sequenceMatrix *aSeqMat,struct sequenceMatrix *bSeqMat, int indexA, int indexB);

int sampleSizeCleanIndex(struct sequenceMatrix *aSeqMat, int index);
char minorAlleleSite(struct sequenceMatrix *aSeqMat, int index);
double majorAlleleFreq(struct sequenceMatrix *aSeqMat, int index);
double siteAssoc(struct sequenceMatrix *aSeqMat, int indexA, int indexB);
void outputPolyMask(struct sequenceMatrix *aSeqMat, int start, int stop, FILE *outfile);
double dij(struct sequenceMatrix *aSeqMat, int indexA, int indexB);
double ZnSFromTo(struct sequenceMatrix *aSeqMat, int start, int stop);
double omegaFromTo(struct sequenceMatrix *aSeqMat, int start, int stop, int site, double **dijTable);
double omegaMaxFromTo(struct sequenceMatrix *aSeqMat, int start, int stop);
double omegaAtCenter(struct sequenceMatrix *aSeqMat, int start, int stop, double site);

double fstFromTo(struct sequenceMatrix *aSeqMat,struct sequenceMatrix *bSeqMat,struct sequenceMatrix *merge,int start,int stop);
void outputPolyMaskBed(struct sequenceMatrix *aSeqMat, int start, int stop, FILE *outfile);

double SnnFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop);
double xij_SnnFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, \
	int seqIndex1, int belongFlag, int start,int stop);

double seqDist_SnnFromTo( stringWrap *seq1,  stringWrap *seq2, int start, int stop);	
double DxyFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop);
double Dxy_minFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop);
double pairwiseDistRankAmongSampleRange(struct sequenceMatrix *aSeqMat, int diffs, double *hetVar, int start, int stop);
double pairwiseIBSMax2PopnFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop);
double pairwiseIBSMeanWithinFromTo(struct sequenceMatrix *aSeqMat, int start,int stop);	
double *hetVec1PopnFromTo(struct sequenceMatrix *aSeqMat, int *vecLen, int start, int end);
double *pairwiseIBSVec1PopnFromTo(struct sequenceMatrix *aSeqMat, int *vecLen, int start, int end);
double sStarSnpDist(struct sequenceMatrix *aSeqMat, int i, int j, char majAlleleI, char majAlleleJ);
double countAlleleDiffsForSnpPair(struct sequenceMatrix *aSeqMat, int i, int j, char majAlleleI, char majAlleleJ);
double sStarFromTo(struct sequenceMatrix *aSeqMat, int start, int end);
int breakClusterAssignmentTie(double **hetMatrix, int targIndex, int *membership1, int m1Index, int *membership2, int m2Index);
int assignClusters(double *hetVec, int n1, int *g1Size, int *membership1, int *g2Size, int *membership2);
void clusterSeqsFromUnsortedHetVec(double *hetVec, struct sequenceMatrix *fullSeqMat, struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat);

void statVecMoments(double *vec, int vecLen, double *mean, double *var, double *skew, double *kurt);
void statVecMinMedMax(double *vec, int vecLen, double *minVal, double *medVal, double *maxVal);

/// WARNING!
//functions below might be defunct!!!
int fixedDiffs(int segsites, int nsam, char **list);
double nucdivIn( int nsam, int segsites, char **list);
double tajd(int, int, double) ;
double hfay(int, int, char **);
double thetah(int, int, char **);
double a1f(int);
double a2f(int);
double b1f(int);
double b2f(int);
double c1f(double, double);
double c2f(int, double, double, double);
double e1f(double, double);
double e2f(double, double, double);
void sort_seqs( char **list,int left, int right);
void swap( void **p1,  void **p2);
int haploCount(char **list, int nsam);
int cmpr(const void *a, const void *b);
int cmp_doubles(const void *x, const void *y);

double piWindow(int nsam, int segsites, char **list, double *locs, double start, double stop);
int compare_doubles(const void *a,const void *b);

	
#endif

