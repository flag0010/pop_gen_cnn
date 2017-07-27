/*****************************
/ pgSummaryStats.c 
/ 
/ A. D. Kern
/
/ popGen statistics for
/ sequenceMatrix objects
******************************/


#include "pgSummaryStats.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "miscCode/vector.h"
#include "miscCode/numerical.h"

#define MAXMUTS 1000

int nHaplotypes(struct sequenceMatrix* aSeqMat, int beg, int end)
{
	int i;
	int j;
	int k;
	int haplotype_found;
	int allsame;
	
	int n_haplotypes = 0;
	char **haplotypes; //haplotypes[aSeqMat->sampleSize][end-beg];
	
	 if( ! ( haplotypes = (char **) malloc( (unsigned)( aSeqMat->sampleSize*sizeof( char* )) ) ) )
		perror("alloc error in haplotypes") ;
	for( i=0; i<aSeqMat->sampleSize; i++) {
		if( ! ( haplotypes[i] = (char *) malloc( (unsigned) ((end-beg)*sizeof( char )) )))
			perror("alloc error in haplotypes 2");
	}

	for(i=0; i<aSeqMat->sampleSize; i++)
	{
		haplotype_found = 0;
		for(j=0; j<n_haplotypes; j++)
		{
			allsame = 1;
			for(k=0; k<end-beg; k++)
			{
				if(haplotypes[j][k] != aSeqMat->matrix[i]->cString[beg+k])
				{
					if(haplotypes[j][k] == 'N' || haplotypes[j][k] == '-')
						haplotypes[j][k] = aSeqMat->matrix[i]->cString[beg+k];
												
					else if(aSeqMat->matrix[i]->cString[beg+k] != 'N' && aSeqMat->matrix[i]->cString[beg+k] != '-')
					{
						allsame = 0;
						break;	
					}
				}
			}
			
			if(allsame)
			{
				haplotype_found = 1;
				break;	
			}
		}
		
		if(!haplotype_found)
		{
			n_haplotypes++;
			for(j=0; j<end-beg; j++)
					haplotypes[n_haplotypes-1][j] = aSeqMat->matrix[i]->cString[beg+j];
		}
	}
	for(i=0;i<aSeqMat->sampleSize;i++){
		free(haplotypes[i]);
	}	
	free(haplotypes);
	return n_haplotypes;
}

//haplotype_counts[i] is the number of haplotypes found in exactly i+1 individuals
int getHaplotypeFreqSpec(struct sequenceMatrix* aSeqMat, int beg, int end, int *haplotype_counts)
{
        int i;
        int j;
        int k;
        int haplotype_found;
        int allsame;
        int freq;

        int n_haplotypes = 0;
        int haplotype_occurrences[aSeqMat->sampleSize];
        char **haplotypes; //haplotypes[aSeqMat->sampleSize][end-beg];
        
        if( ! ( haplotypes = (char **) malloc( (unsigned)( aSeqMat->sampleSize*sizeof( char* )) ) ) )
                perror("alloc error in haplotypes") ;
        for( i=0; i<aSeqMat->sampleSize; i++) {
                if( ! ( haplotypes[i] = (char *) malloc( (unsigned) (((end-beg)+1)*sizeof( char )) )))
                        perror("alloc error in haplotypes 2");
        }

        for(i=0; i<aSeqMat->sampleSize; i++)
        {
            haplotype_counts[i] = 0;
            haplotype_occurrences[i] = 0;
        }

        for(i=0; i<aSeqMat->sampleSize; i++)
        {
                haplotype_found = 0;
                for(j=0; j<n_haplotypes; j++)
                {
                        allsame = 1;
                        for(k=0; k<end-beg; k++)
                        {
                                if(haplotypes[j][k] != aSeqMat->matrix[i]->cString[beg+k])
                                {
                                        if(haplotypes[j][k] == 'N' || haplotypes[j][k] == '-')
                                                haplotypes[j][k] = aSeqMat->matrix[i]->cString[beg+k];
                                                                                                
                                        else if(aSeqMat->matrix[i]->cString[beg+k] != 'N' && aSeqMat->matrix[i]->cString[beg+k] != '-')
                                        {
                                                allsame = 0;
                                                break;  
                                        }
                                }
                        }

                        if(allsame)
                        {
                                haplotype_found = 1;
                                haplotype_occurrences[j]+=1;
                                break;
                        }
                }
                if(!haplotype_found)
                {
                        n_haplotypes++;
                        for(j=0; j<end-beg; j++)
                                        haplotypes[n_haplotypes-1][j] = aSeqMat->matrix[i]->cString[beg+j];
                        haplotypes[n_haplotypes-1][end-beg]='\0';
                        haplotype_occurrences[n_haplotypes-1]=1;
                }
        }

        for (i=0; i<n_haplotypes; i++)
        {
                freq = haplotype_occurrences[i];
                if (freq > 0 && freq <= aSeqMat->sampleSize)
                {
                        haplotype_counts[freq-1] += 1;
                }
        }

        for(i=0;i<aSeqMat->sampleSize;i++){
                free(haplotypes[i]);
        }       
        free(haplotypes);
        return n_haplotypes;
}

double petrovH1(int *haplotype_counts, int nsam)
{
    int hapFreq;
    double pi;
    double h1 = 0.0;

    for (hapFreq=nsam; hapFreq>0; hapFreq--)
    {
        pi = hapFreq/ (double)nsam;
        h1 += haplotype_counts[hapFreq-1]*pi*pi;
    }
    return h1;
}

double petrovH2(int *haplotype_counts, int nsam)
{
    int hapFreq;
    double pi;
    double h2 = 0.0;
    int first = 1;

    for (hapFreq=nsam; hapFreq>0; hapFreq--)
    {
        pi = hapFreq/ (double)nsam;
        if (haplotype_counts[hapFreq-1] > 0)
        {
            if (first)
            {
                first = 0;
                h2 += (haplotype_counts[hapFreq-1]-1)*pi*pi;
            }
            else
            {
                h2 += haplotype_counts[hapFreq-1]*pi*pi;
            }
        }
    }
    return h2;
}

double petrovH12(int *haplotype_counts, int nsam)
{
    int hapFreq, i;
    double pi;
    double part1 = 0.0;
    double part2 = 0.0;
    int totalAdded = 0;

    for (hapFreq=nsam; hapFreq>0; hapFreq--)
    {
        pi = hapFreq/ (double)nsam;
        for (i = 0;i < haplotype_counts[hapFreq-1];i++)
        {
            if (totalAdded < 2)
            {
                part1 += pi;
            }
            else
            {
                part2 += pi*pi;

            }
            totalAdded++;
        }
    }

    part1 = part1*part1;
    return part1+part2;
}

double maxFDA(struct sequenceMatrix *aSeqMat){
    double mfda = 0.0;
    mfda = maxFDAFromTo(aSeqMat,0,aSeqMat->length);
    return(mfda);
}

double nucdiv(struct sequenceMatrix *aSeqMat){
	double pi = 0.0;
	pi = nucdivFromTo(aSeqMat,0,aSeqMat->length);
	return(pi);	
}

// fixed differences (set is ingroup data)
int fixedDiffsPre(stringWrap** set, int size, int samplesize, stringWrap* outgdata)
{
	int i;
	int fixeddiffs = 0;

	for(i=0; i<size; i++)
	{
		if(set[i]->length == 1)
		{
			/* if outgdata != set[0] it's a fixed difference */
			if(outgdata->cString[0] != '-' && outgdata->cString[0] != 'N' && outgdata->cString[0] != set[i]->cString[0])
				fixeddiffs++;
		}
	}

	return fixeddiffs;
}

// proportion of missing data
double missingDataPre(stringWrap** array, int size, int samplesize)
{
	int i;
	int total;
	int present;

	total = samplesize * size;
	present = 0;

	for(i=0; i<size; i++)
		present += array[i]->length;

	return (double)(total - present) / (double)(total);
}

double missingDataFromTo(struct sequenceMatrix* aSeqMat, int start, int end)
{
	int i;
	stringWrap* array = stringWrapNew(aSeqMat->sampleSize);

	int total = (end-start) * aSeqMat->sampleSize;
	int present = 0;

	for(i=start; i<end; i++)
	{
		sequenceMatrix_siteArrayClean(aSeqMat, array, i);
		present += array->length;
	}

	stringWrapFree(array);
	return (double)(total - present) / (double)total;
}

/* calculates pi (nucdiv) using pre-calculated set and array */
double nucdivPre(stringWrap** set, stringWrap** array, int size)
{
	int i, j;
	double pi, p1, n, nnm1, sum;

	pi  = 0.0;
	sum = 0.0;

	for(i=0; i<size; i++)
	{
		if(set[i]->length > 1)
		{
			n = (double)(array[i]->length);
			nnm1 = n / (n-1.0);
			p1 = 0.0;
			pi = 0.0;

			for(j=0; j<set[i]->length; j++)
			{
				p1 = (double)(stringWrapCountChar(array[i], set[i]->cString[j])) / n;
				pi += p1*p1;
			}
			
			sum += (1.0-pi) * nnm1;
		}
	}

	return sum;
}

int breakClusterAssignmentTie(double **hetMatrix, int targIndex, int *membership1, int m1Index, int *membership2, int m2Index){
        double d1 = 0., d2 = 0.;
        int i;
        for (i=0; i<m1Index; i++)
        {
                d1 += hetMatrix[targIndex][membership1[i]];
        }
        d1 = d1/m1Index;

        for (i=0; i<m2Index; i++)
        {
                d2 += hetMatrix[targIndex][membership2[i]];
        }
        d2 = d2/m2Index;

        /*if (abs(d1 - d2) < 1e-9)
        {
                if (rand()/(double) RAND_MAX + 1 > 0.5)
                {
                        return 1;
                }
                else
                {
                        return 2;
                }
        }*/
        //else if (d1 > d2)
        if (d1 > d2)
        {
                return 2;
        }
        else
        {
                return 1;
        }
}

//returns 1 if we have to swap the resulting clusters (in order to make the first cluster the larger one
//which we assume to be the case later one).
int assignClusters(double *hetVec, int n1, int *g1Size, int *membership1, int *g2Size, int *membership2){
        int i, j, minI = -1, maxPairI = -1, maxPairJ = -1, pairIndex = 0, m1Index = 0, m2Index = 0;
        double minHet, maxHet = -1;
        double **hetMatrix;
        //printf("n1: %d\n", n1);
        hetMatrix = (double **) malloc( (unsigned)( n1*sizeof( double* )));
        for(i=0; i<n1; i++) hetMatrix[i] = (double *) malloc( (unsigned) (n1*sizeof( double )));
        for(i=0; i<n1-1;i++)
        {
                //printf("i: %d\n", i);
                for(j=i+1;j<n1;j++)
                {
                        hetMatrix[i][j] = hetVec[pairIndex];
                        hetMatrix[j][i] = hetVec[pairIndex];
                        //printf("hetVec[%d][%d]: %f; maxHet: %f\n", i, j, hetVec[pairIndex], maxHet);
                        if (hetVec[pairIndex] > maxHet)
                        {
                                //printf("maxHet: %f\n", maxHet);
                                maxHet = hetVec[pairIndex];
                                maxPairI = i;
                                maxPairJ = j;
                        }
                        pairIndex++;
                }
        }
        membership1[m1Index++] = maxPairI;
        membership2[m2Index++] = maxPairJ;
        //printf("%d: 1\n", maxPairI);
        //printf("%d: 2\n", maxPairJ);

        for(i=0; i<n1;i++)
        {
                if (i != maxPairI && i != maxPairJ)
                {
                        if (hetMatrix[i][maxPairI] > hetMatrix[i][maxPairJ] - 1e-9)
                        {
                                //printf("%d: 2, %f, %f\n", i, hetMatrix[i][maxPairI], hetMatrix[i][maxPairJ]);
                                membership2[m2Index++] = i;
                        }
                        else if (hetMatrix[i][maxPairI] < hetMatrix[i][maxPairJ] - 1e-9 || breakClusterAssignmentTie(hetMatrix, i, membership1, m1Index, membership2, m2Index) == 1)
                        {
                                //printf("%d: 1, %f, %f\n", i, hetMatrix[i][maxPairI], hetMatrix[i][maxPairJ]);
                                membership1[m1Index++] = i;
                        }
                        else
                        {
                                //printf("%d: 2, %f, %f\n", i, hetMatrix[i][maxPairI], hetMatrix[i][maxPairJ]);
                                membership2[m2Index++] = i;
                        }
                }
        }
        if (m1Index == 1 && m2Index > 1)
        {
                minHet = 10e10;
                //printf("not enough samples in group 1, going to take one from group 2\n");
                for (i=0; i<m2Index; i++)
                {
                        if (hetMatrix[membership1[0]][membership2[i]] < minHet)
                        {
                                minHet = hetMatrix[membership1[0]][membership2[i]];
                                minI = i;
                        }
                }
                membership1[m1Index++] = membership2[minI];
                //printf("put %d into group 1.\n", membership2[minI]);
                if (minI < m2Index-1)
                {
                        //printf("moved group 2's %dth element into the %dth slot.\n", m2Index-1, minI);
                        membership2[minI] = membership2[m2Index-1];
                }
                m2Index--;
        }
        else if (m2Index == 1 && m1Index > 1)
        {
                minHet = 10e10;
                //printf("not enough samples in group 2, going to take one from group 1\n");
                for (i=0; i<m1Index; i++)
                {
                        if (hetMatrix[membership2[0]][membership1[i]] < minHet)
                        {
                                minHet = hetMatrix[membership2[0]][membership1[i]];
                                minI = i;
                        }
                }
                membership2[m2Index++] = membership1[minI];
                //printf("put %d into group 2.\n", membership1[minI]);
                if (minI < m1Index-1)
                {
                        //printf("moved group 1's %dth element into the %dth slot.\n", m1Index-1, minI);
                        membership1[minI] = membership1[m1Index-1];
                }
                m1Index--;
        }
        for(i=0; i<n1; i++) free(hetMatrix[i]);
        free(hetMatrix);
        *g1Size = m1Index;
        *g2Size = m2Index;
        //for(i=0;i<*(g1Size);i++) printf("pop 1: %d\n", membership1[i]);
        //for(i=0;i<*(g2Size);i++) printf("pop 2: %d\n", membership2[i]);
        //printf("m1Index: %d; m2Index: %d\n", m1Index, m2Index);
        if (m1Index < m2Index)
        {
                //printf("gotta swap!!\n");
                return 1;
        }
        else
        {
                return 0;
        }
}

//clustering of sequences following Hammer et al. (2011)
//note: in current implementation, the lements in aSeqMat and bSeqMat (both the matrix and the names) are merely references to the elements of fullSeqMat (not clones!)
//so changes to fullSeqMat will effect aSeqMat and bSeqMat, and vice-versa
void clusterSeqsFromUnsortedHetVec(double *hetVec, struct sequenceMatrix *fullSeqMat, struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat)
{
	int g1Size, g2Size;
        int *membership1 = (int *) malloc (fullSeqMat->sampleSize*sizeof(int));
        int *membership2 = (int *) malloc (fullSeqMat->sampleSize*sizeof(int));
        int i, swap, tmp, *tmpArray;
        swap = assignClusters(hetVec, fullSeqMat->sampleSize, &g1Size, membership1, &g2Size, membership2);
        if (swap)
        {
                tmp = g1Size;
                g1Size = g2Size;
                g2Size = tmp;
                tmpArray = membership1;
                membership1 = membership2;
                membership2 = tmpArray;
        }
        sequenceMatrix_newByRef(g1Size, fullSeqMat->length, aSeqMat);
        for (i=0; i<g1Size; i++)
	{
		//printf("g1: %d\n", membership1[i]);
		aSeqMat->matrix[i] = fullSeqMat->matrix[membership1[i]];
		aSeqMat->names[i] = fullSeqMat->names[membership1[i]];
	}
        free(membership1);
        sequenceMatrix_newByRef(g2Size, fullSeqMat->length, bSeqMat);
        for (i=0; i<g2Size; i++)
	{
		//printf("g2: %d\n", membership2[i]);
		bSeqMat->matrix[i] = fullSeqMat->matrix[membership2[i]];
		bSeqMat->names[i] = fullSeqMat->names[membership2[i]];
	}
	assert(fullSeqMat->sampleSize == aSeqMat->sampleSize + bSeqMat->sampleSize);
        free(membership2);
}

//redundant with piWindow below but easier to use
double nucdivFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
	int i, j;
  	double pi, p1, n, nnm1,sum ;
	stringWrap *set, *array;

  	pi = 0.0 ;
	sum =0.0;
	set = stringWrapNew(aSeqMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		if(set->length > 1){
			sequenceMatrix_siteArrayClean(aSeqMat,array,i);
			n = (double) array->length;
			nnm1 = n / (n-1.0);
			p1 = 0.;
			pi = 0.;
			for(j=0; j < set->length; j++){
				p1 = (double) stringWrapCountChar(array,set->cString[j]) / n;
				pi += (p1*p1);
			}
			sum += (1.0-pi) * nnm1;
  		}
	}
	stringWrapFree(set);
	stringWrapFree(array);
	return(sum);
}

double maxFDAFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
        int i, j;
        double mfda, p1, n;
        stringWrap *set, *array;

        mfda = 0.0;
        set = stringWrapNew(aSeqMat->sampleSize);
        array = stringWrapNew(aSeqMat->sampleSize);
        for( i = start; i < end; i++){
                sequenceMatrix_siteSetClean(aSeqMat,set,i);
                if(set->length > 1){
                        sequenceMatrix_siteArrayClean(aSeqMat,array,i);
                        n = (double) array->length;
                        p1 = 0.;
                        for(j=0; j < set->length; j++){
                                p1 = (double) stringWrapCountChar(array,set->cString[j]) / n;
                                if (p1 < 1.0 && p1 > mfda){
                                    mfda = p1;
                                }
                        }
                }
        }
        stringWrapFree(set);
        stringWrapFree(array);
        return(mfda);
}

int numColumnsNotNedOutFromToBothMats(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start, int end){
        int i, numSites;
        stringWrap *array1, *array2;
        array1 = stringWrapNew(aSeqMat->sampleSize);
        array2 = stringWrapNew(bSeqMat->sampleSize);
	numSites = 0;

        for( i = start; i < end; i++){
                sequenceMatrix_siteArrayClean(aSeqMat,array1,i);
                sequenceMatrix_siteArrayClean(bSeqMat,array2,i);
                if (array1->length > 1 || array2->length > 1)
                {
                    numSites += 1;
                }
        }

        stringWrapFree(array1);
        stringWrapFree(array2);

        return numSites;
}

int numColumnsNotNedOutFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
        int i, numSites;
        stringWrap *array;

        array = stringWrapNew(aSeqMat->sampleSize);
	numSites = 0;

        for( i = start; i < end; i++){
                sequenceMatrix_siteArrayClean(aSeqMat,array,i);
                if (array->length > 1)
                {
                    numSites += 1;
                }
        }

        stringWrapFree(array);

        return numSites;
}

//Fay's Theta_H 
double thetaHPre(stringWrap** set, stringWrap** array,stringWrap** ancSet , int size)
{
	int i;
	double  p1, n, nnm1, sum;
	char derAllele;

	sum = 0.0;

	for(i=0; i<size; i++){
		if(set[i]->length == 2 && ancSet[i]->length == 1 && (set[i]->cString[0] == ancSet[i]->cString[0]  || set[i]->cString[1] == ancSet[i]->cString[0] )){
			//find derivedAllele
			derAllele = (ancSet[i]->cString[0] == set[i]->cString[0]) ? set[i]->cString[1] : set[i]->cString[0];
			n = (double) array[i]->length;
			nnm1 = n * (n - 1.0);
			p1 = (double) stringWrapCountChar(array[i],derAllele);
			sum += (p1*p1) / nnm1;
		}
	}

	return (2.0*sum);
}
//this needs an outgorup to unfold SFS. Assuming that the ancestral sequence is avail
double thetaHFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *ancMat, int start, int end){
	int i;
	double p1, n, nnm1,sum ;
	stringWrap *set, *array, *ancSet;
	char derAllele;

	sum =0.0;
	set = stringWrapNew(aSeqMat->sampleSize);
	ancSet= stringWrapNew(ancMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		sequenceMatrix_siteSetClean(ancMat,ancSet,i);
		//two alleles in ingroup and ancRecon at site?
		if(set->length == 2 && ancSet->length == 1 && (set->cString[0] == ancSet->cString[0]  || set->cString[1] == ancSet->cString[0] )){
			sequenceMatrix_siteArrayClean(aSeqMat,array,i);
			//find derivedAllele
			derAllele = (ancSet->cString[0] == set->cString[0]) ? set->cString[1] : set->cString[0];
			n = (double) array->length;
			nnm1 = n * (n - 1.0);
			p1 = (double) stringWrapCountChar(array,derAllele);
			sum += (p1*p1) / nnm1;
		}
	}
	stringWrapFree(set);
	stringWrapFree(ancSet);
	stringWrapFree(array);
	return(2.0*sum);
}

double thetaH(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *ancMat){
	return(thetaHFromTo(aSeqMat,ancMat,0,aSeqMat->length));
}

//Waterson's theta stuff
//
int segSiteCountPre(stringWrap** set, int size)
{
	int i, count;

	count = 0;
	for(i=0; i<size; i++)
	{
		if(set[i]->length > 1)
			count++;
	}

	return count;
}

int segSiteCountFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
	int i,count;
	stringWrap *set;

  	count = 0;
	set = stringWrapNew(aSeqMat->sampleSize);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		if(set->length > 1){
			count+=1;
  		}
	}
	stringWrapFree(set);
	return(count);
}

int segSiteCountSubpopFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start, int end){
	int i,count;
	stringWrap *set;

  	count = 0;
	set = stringWrapNew(aSeqMat->sampleSize+bSeqMat->sampleSize);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetSubpopClean(aSeqMat,bSeqMat,set,i);
		if(set->length > 1){
			count+=1;
  		}
	}
	stringWrapFree(set);
	return(count);
}

void privateSegSitesInTwoPopnsFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int *private1, int *private2, int start, int end){
	int i;
	stringWrap *set1, *set2;
	(*private1) = 0, (*private2) = 0;

	set1 = stringWrapNew(aSeqMat->sampleSize);
	set2 = stringWrapNew(bSeqMat->sampleSize);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set1,i);
		sequenceMatrix_siteSetClean(bSeqMat,set2,i);
		if (set1->length > 1 && set2->length <= 1){
			(*private1)++;
		}
		else if (set2->length > 1 && set1->length <= 1){
			(*private2)++;
		}
	}
	stringWrapFree(set1);
	stringWrapFree(set2);
}

int segSiteCount(struct sequenceMatrix *aSeqMat){
	return(segSiteCountFromTo(aSeqMat,0,aSeqMat->length));
}

void segSiteLocationsFromTo(struct sequenceMatrix *aSeqMat, vector *segs, int start, int end){
	int i;
	stringWrap *set;
	set = stringWrapNew(aSeqMat->sampleSize + 1);
	vectorInit(segs);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		if(set->length > 1){
			vectorAppend(segs,(double) i);
  		}
	}
	stringWrapFree(set);
}

void segSiteLocationsSubpopFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, vector *segs, int start, int end){
	int i;
	stringWrap *set;
	set = stringWrapNew(aSeqMat->sampleSize + bSeqMat->sampleSize + 1);
	vectorInit(segs);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetSubpopClean(aSeqMat,bSeqMat,set,i);
		if(set->length > 1){
			vectorAppend(segs,(double) i);
  		}
	}
	stringWrapFree(set);
}

//segSiteLocationsBiallelicFromTo 
void segSiteLocationsBiallelicFromTo(struct sequenceMatrix *aSeqMat, vector *segs, int start, int end){
	int i;
	stringWrap *set;
	set = stringWrapNew(aSeqMat->sampleSize + 1);
	vectorInit(segs);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		if(set->length == 2){
			vectorAppend(segs,(double) i);
  		}
	}
	stringWrapFree(set);
}
void segSiteLocations(struct sequenceMatrix *aSeqMat, vector *segs){
	segSiteLocationsFromTo(aSeqMat,segs,0,aSeqMat->length);
}

void segSiteBiallelicLocations(struct sequenceMatrix *aSeqMat, vector *segs){
	segSiteLocationsBiallelicFromTo(aSeqMat,segs,0,aSeqMat->length);
}

double thetaWPre(stringWrap** set, int size, int samplesize)
{
	int s;
	s = segSiteCountPre(set, size);
	return (double)(s) / a1f(samplesize);
}

double thetaWFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
	int s;	
	s = segSiteCountFromTo(aSeqMat,start,end);
	return((double)s/a1f(aSeqMat->sampleSize));
}

double thetaW(struct sequenceMatrix *aSeqMat){
	return(thetaWFromTo(aSeqMat,0,aSeqMat->length));
}

//wrapper for tajd function below
double tajD(struct sequenceMatrix *aSeqMat){
	return(tajd(aSeqMat->sampleSize, segSiteCount(aSeqMat), nucdiv(aSeqMat)));
}

double tajDPre(stringWrap** set, stringWrap** array, int size, int samplesize)
{
	double pi = nucdivPre(set, array, size);
	int     s = segSiteCountPre(set, size);
	return tajd(samplesize, s, pi);
}

double tajDFromTo(struct sequenceMatrix *aSeqMat,int start, int end){
	double pi=nucdivFromTo(aSeqMat,start,end);
	int s=segSiteCountFromTo(aSeqMat,start,end);
	return(tajd(aSeqMat->sampleSize,s,pi));
}

double frequency( stringWrap *array, char aChar){
	int i;
	double count=0.0;
  for( i=0; i<array->length; i++) count += ( array->cString[i] == aChar ? 1.0: 0.0 ) ;
  return( count/array->length);
}        

int count( stringWrap *array, char aChar){
	int i;
	int count=0;
  	for( i=0; i<array->length; i++) count += ( array->cString[i] == aChar ? 1.0: 0.0 ) ;
  	return( count);
}

int sampleSizeCleanIndex(struct sequenceMatrix *aSeqMat, int index){
	stringWrap *array;
	int count = 0;
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteArrayClean(aSeqMat,array,index);
	count = array->length;
	stringWrapFree(array);
	return(count);
}

double freqAlleleSite(struct sequenceMatrix *aSeqMat, char aChar, int index){
	stringWrap *array;
	double count;
	
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteArrayClean(aSeqMat,array,index);
	count = stringWrapCountChar(array,aChar);
	stringWrapFree(array);
	return(count/array->length);
}

//returns most frequency allele
char majorAlleleSite(struct sequenceMatrix *aSeqMat, int index){
	stringWrap *set, *array;
	int i, current , max, tmp;
	char c;
	
	current = max = 0;
	set= stringWrapNew(aSeqMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteSetClean(aSeqMat, set,index);
	sequenceMatrix_siteArrayClean(aSeqMat, array,index);
	for(i=0;i<set->length;i++){
		tmp = stringWrapCountChar(array,set->cString[i]);
		if(max < tmp){
			max = tmp;
			current = i;
		}
	}
	c = set->cString[current];
	stringWrapFree(set);
	stringWrapFree(array);
	return(c);
}

//returns most frequency allele
double majorAlleleFreq(struct sequenceMatrix *aSeqMat, int index){
	stringWrap *set, *array;
	int i, max, tmp;
	double p;
	
	max = 0;
	set= stringWrapNew(aSeqMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteSetClean(aSeqMat, set,index);
	sequenceMatrix_siteArrayClean(aSeqMat, array,index);
	for(i=0;i<set->length;i++){
		tmp = stringWrapCountChar(array,set->cString[i]);
		if(max < tmp){
			max = tmp;
		}
	}
	p = (double) max / array->length;
	stringWrapFree(set);
	stringWrapFree(array);
	return(p);
}

//returns least frequency allele
char minorAlleleSite(struct sequenceMatrix *aSeqMat, int index){
	stringWrap *set, *array, *tmpSW;
	int i, current , min, tmp;
	char c;
	
	current = 0;
	min = aSeqMat->sampleSize;
	set= stringWrapNew(aSeqMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteSetClean(aSeqMat, set,index);
	sequenceMatrix_siteArrayClean(aSeqMat, array,index);
	for(i=0;i<set->length;i++){
		tmp = stringWrapCountChar(array,set->cString[i]);
		if(min > tmp){
			min = tmp;
			current = i;
		}
	}
	c = set->cString[current];
	//make sure this isn't the same allele if freq = 0.5
	if(((float)min/array->length) == 0.5){
		tmpSW = stringWrapFindOther(set, majorAlleleSite(aSeqMat,index));
		c = tmpSW->cString[0];
		stringWrapFree(tmpSW);
	}
	stringWrapFree(set);
	stringWrapFree(array);
	return(c);
}

int min_rec (struct sequenceMatrix *aSeqMat, int x){
 // Calculate min # rec. events 
  	int a, b, c, e, gtest, flag = 0;
	int size, nsam;
	vector *locs;
	char majAlleleA, majAlleleB;
	
	size = segSiteCount(aSeqMat);
	nsam = aSeqMat->sampleSize;
	locs = vectorNew(size);
	segSiteLocations(aSeqMat,locs);

	if (size<2 || x >= (size-1))
		return (0);
	for (a=x+1; a<size; ++a) {
		for (b=x; b<a; ++b) {
			gtest = 0;
			majAlleleA = majorAlleleSite(aSeqMat,(int)vectorGetIndex(locs,a));
			majAlleleB = majorAlleleSite(aSeqMat,(int)vectorGetIndex(locs,b));
			for (e=0; e<nsam; ++e)
			if (aSeqMat->matrix[e]->cString[b] != majAlleleB && aSeqMat->matrix[e]->cString[a] != majAlleleA) {
				++gtest;
				break;
			}
			for (e=0; e<nsam; ++e)
			if (aSeqMat->matrix[e]->cString[b] != majAlleleB && aSeqMat->matrix[e]->cString[a] == majAlleleA) {
				++gtest;
				break;
			}
			for (e=0; e<nsam; ++e)
			if (aSeqMat->matrix[e]->cString[b] == majAlleleB && aSeqMat->matrix[e]->cString[a] != majAlleleA) {
				++gtest;
				break;
			}
			for (e=0; e<nsam; ++e)
			if (aSeqMat->matrix[e]->cString[b] == majAlleleB && aSeqMat->matrix[e]->cString[a] == majAlleleA) {
				++gtest;
				break;
			}       
			if (gtest == 4) {
				flag = 1;
				break;
			}
		}
		if (flag == 1)
			break;
	}
	if (a==size)
		return (0);
	else {
		c = min_rec(aSeqMat, a);
		return (1+c);
	}
}

double rSquared(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB, d, denom;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
//	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
//	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}

	pAB /= comp;
	pA /= comp;
	pB /= comp;
	d = pAB - (pA * pB);
	denom = sqrt(pA * (1.0 - pA) * pB * (1.0-pB));
	d /= denom;
	return(d*d);
}


double rSquaredOmega(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB, d, denom,sign;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
//	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
//	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}

	pAB /= comp;
	pA /= comp;
	pB /= comp;
	d = pAB - (pA * pB);
	if(d > 0){
		sign = 1.0;
	}
	else{
		sign = -1.0;
	}
	denom = sqrt(pA * (1.0 - pA) * pB * (1.0-pB));
	d /= denom;
	return (sign * (d*d));
}

void rSquaredCounts(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
	comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}

	pAB /= comp;
	pA /= comp;
	pB /= comp;
	printf("%f\t%f\t%f\t%d",pAB,pA,pB,comp);
}

//jointHeterozygosity returns p(1-p)q(1-q) the squared denom of r^2
double jointHeterozygosity(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, denom;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
//	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
//	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
		}
	}

	pA /= comp;
	pB /= comp;
	denom = pA * (1.0 - pA) * pB * (1.0-pB);
	return(denom);
}

double rSquared_2chroms(struct sequenceMatrix *aSeqMat,struct sequenceMatrix *bSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB, d, denom;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(bSeqMat,indexB);
//	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
//	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(bSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(bSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(bSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}

	pAB /= comp;
	pA /= comp;
	pB /= comp;
	d = pAB - (pA * pB);
	denom = sqrt(pA * (1.0 - pA) * pB * (1.0-pB));
	d /= denom;
	return(d*d);
}

void rSquaredCounts_2chroms(struct sequenceMatrix *aSeqMat,struct sequenceMatrix *bSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(bSeqMat,indexB);
	comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(bSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(bSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(bSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}

	pAB /= comp;
	pA /= comp;
	pB /= comp;
	printf("%f\t%f\t%f\t%d",pAB,pA,pB,comp);
}

double countAlleleDiffsForSnpPair(struct sequenceMatrix *aSeqMat, int i, int j, char majAlleleI, char majAlleleJ)
{
        int k, numDiffs = 0, numSames = 0;
        for (k=0; k<aSeqMat->sampleSize; k++)
        {
                if (aSeqMat->matrix[k]->cString[i] != 'N' && aSeqMat->matrix[k]->cString[j] != 'N')
                {
                    if ((aSeqMat->matrix[k]->cString[i] != majAlleleI && aSeqMat->matrix[k]->cString[j] == majAlleleJ) || (aSeqMat->matrix[k]->cString[i] == majAlleleI && aSeqMat->matrix[k]->cString[j] != majAlleleJ))
                    {
                        numDiffs++;
                    }
                    else
                    {
                        numSames++;
                    }
                }
        }
        if (numDiffs < numSames)
            return numDiffs;
        else
            return numSames;
}

double sStarSnpDist(struct sequenceMatrix *aSeqMat, int i, int j, char majAlleleI, char majAlleleJ)
{
        int numDiffs;
        numDiffs = countAlleleDiffsForSnpPair(aSeqMat, i, j, majAlleleI, majAlleleJ);
        if (numDiffs == 0)
        {
                return (i-j) + 5000;
        }
        else if (numDiffs < 5)
        {
                return -10000.0;
        }
        else
        {
                return -10e10;
        }
}

double sStarFromTo(struct sequenceMatrix *aSeqMat, int start, int end)
{
        int i,j,snpIndexI,snpIndexJ,segsites;
        char majAlleleI, majAlleleJ;
        double globalMaxVal, maxVal, currVal;
	double *sStarArray;
	vector *locs;

        globalMaxVal = 0.0;//gotta get me some snp locations or some shit
        segsites = segSiteCountFromTo(aSeqMat,start,end);
        locs = vectorNew(segsites);
        segSiteLocationsFromTo(aSeqMat,locs,start,end);
        sStarArray = (double *)malloc( sizeof(double)*segsites );
        sStarArray[0] = 0.0;
        for (snpIndexI=1; snpIndexI<segsites; snpIndexI++)
        {
                i = vectorGetIndex(locs,snpIndexI);
                majAlleleI = majorAlleleSite(aSeqMat,i);
                maxVal = 0.0;
                for (snpIndexJ=0; snpIndexJ<snpIndexI; snpIndexJ++)
                {
                        j = vectorGetIndex(locs,snpIndexJ);
                        majAlleleJ = majorAlleleSite(aSeqMat,j);
                        currVal = sStarArray[snpIndexJ] + sStarSnpDist(aSeqMat, i, j, majAlleleI, majAlleleJ);
                        //printf("%d, %d; %f, currVal: %f\n", snpIndexI, snpIndexJ, sStarArray[snpIndexJ], currVal);
                        if (currVal > maxVal)
                        {
                                maxVal = currVal;
                        }
                }
                sStarArray[snpIndexI] = maxVal;
                if (maxVal > globalMaxVal)
                {
                        globalMaxVal = maxVal;
                }
        }
        vectorFree(locs);
        free(sStarArray);
        return globalMaxVal;
}

//dij for ZnS
double dij(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB, d, denom, comp;
	int i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}
	if(comp < 1){
		return(0);
	}
	else{
		pA /= comp;
		pB /= comp;
		pAB /= comp;
		d = pAB - (pA * pB);
		denom = pA * (1.0 - pA) * pB * (1.0-pB);
		if(denom == 0)
			return(0);
		return((d*d)/denom);
	}
}

double ZnSFromTo(struct sequenceMatrix *aSeqMat, int start, int stop){
	int i, j, size;
	vector *locs;
	double sum= 0.0;
	
	size = segSiteCountFromTo(aSeqMat,start,stop);
	if(size < 2)
		return(0.0);
	locs = vectorNew(size);
	segSiteLocationsFromTo(aSeqMat,locs,start,stop);
	for(i = 0 ; i < size - 1;i++){
		for(j = i+1; j < size; j++){
			sum += dij(aSeqMat,vectorGetIndex(locs,i),vectorGetIndex(locs,j));
			//printf("test: %f\n",sum);
		}
	}
	vectorFree(locs);
	return ((2.0 / (double)(size * (size-1))) * sum);
}

double omegaFromTo(struct sequenceMatrix *aSeqMat, int start, int stop, int site, double **dijTable){
        int i,j, segsites, s;
        vector *locs;
        double sum,sumL,sumR,comp,denom;
        double numer;

        segsites = segSiteCountFromTo(aSeqMat,start,stop);
        s = segsites;
        sum = sumL = sumR =comp=denom=0;
        if(segsites < 3)
                return 0.0;

        locs = vectorNew(segsites);
        segSiteLocationsFromTo(aSeqMat,locs,start,stop);
        //calculate: 
        // sum for denom-- all pairwise r2
        // sumL and sumR for numerator -- blockwise r2
        for(i=0; i<segsites-1; i++){
                for(j=i+1; j<segsites; j++){
                        comp = dijTable[i][j];
                        if(i < site && j >= site)sum += comp;
                        if(i < site && j < site) sumL += comp;
                        if(i >= site && j >= site) sumR += comp;
                }
        }
        denom = sum * (1.0/(site*(s-site)));
        numer = 1.0 / ((site*(site-1)/2) + ((s-site)*(s-site-1)/2));
        numer *= (sumL+sumR);
	vectorFree(locs);
        return(numer/denom);
}

double rEHHatFocalSnp(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *ancMat, int start, int stop, int focalSnpPos)
{
    int i, j, k, l;
    int allsame, currSnpPos, lPos, foundFocalSnp, segsites;
    vector *locs;
    char jBase, kBase, ancBase, derBase;
    int aCount, cCount, gCount, tCount;
    int sameCountDer, sameCountAnc;
    int totalCountDer, totalCountAnc;
    sameCountDer = sameCountAnc = totalCountDer = totalCountAnc = foundFocalSnp = 0;

    assert(ancMat->sampleSize == 1);
    ancBase = sequenceMatrix_rowColumn(ancMat, 0, focalSnpPos);
    if (ancBase == 'A' || ancBase == 'C' || ancBase == 'G' || ancBase == 'T')
    {
        segsites = segSiteCountFromTo(aSeqMat,start,stop);
        locs = vectorNew(1000);
        segSiteLocationsFromTo(aSeqMat, locs, start, stop);
        for(i=0;i<segsites;i++)
        {
            currSnpPos = vectorGetIndex(locs, i);
            if (currSnpPos == focalSnpPos)
                foundFocalSnp = 1;
        }
        //if (!foundFocalSnp){
        //    return NAN;
        //}
        //figure out what the derived allele is
        aCount = cCount = gCount = tCount = 0;
        for(j=0; j < aSeqMat->sampleSize; j++)
        {
            jBase = sequenceMatrix_rowColumn(aSeqMat, j, focalSnpPos);
            if (jBase == 'A' && jBase != ancBase)
                aCount++;
            if (jBase == 'C' && jBase != ancBase)
                cCount++;
            if (jBase == 'G' && jBase != ancBase)
                gCount++;
            if (jBase == 'T' && jBase != ancBase)
                tCount++;
        }
        if (aCount > cCount && aCount > gCount && aCount > tCount)
            derBase = 'A';
        else if (cCount > aCount && cCount > gCount && cCount > tCount)
            derBase = 'C';
        else if (gCount > aCount && gCount > cCount && gCount > tCount)
            derBase = 'G';
        else
            derBase = 'T';
        for(j=0; j < aSeqMat->sampleSize; j++)
        {
            if (sequenceMatrix_rowColumn(aSeqMat, j, focalSnpPos) == derBase)
            {
                for(k=0;k<j;k++)
                {
                    if (sequenceMatrix_rowColumn(aSeqMat, k, focalSnpPos) == derBase)
                    {
                        allsame = 1;
                        //see if j and k are the same at all segsites in the range
                        for(l=0; l<segsites; l++)
                        {
                            lPos = vectorGetIndex(locs, l);
                            if (lPos != focalSnpPos)
                            {
                                jBase = sequenceMatrix_rowColumn(aSeqMat, j, lPos);
                                kBase = sequenceMatrix_rowColumn(aSeqMat, k, lPos);
                                if (jBase != 'N' && kBase != 'N' && jBase != kBase)
                                {
                                    allsame = 0;
                                    break;
                                }
                            }
                        }
                        if (allsame)
                        {
                            sameCountDer++;
                        }
                        totalCountDer++;
                    }
                }
            }
            else if (sequenceMatrix_rowColumn(aSeqMat, j, focalSnpPos) == ancBase)
            {
                for(k=0;k<j;k++)
                {
                    if (sequenceMatrix_rowColumn(aSeqMat, k, focalSnpPos) == ancBase)
                    {
                        allsame = 1;
                        //see if j and k are the same at all segsites in the range
                        for(l=0; l<segsites; l++)
                        {
                            lPos = vectorGetIndex(locs, l);
                            if (lPos != focalSnpPos)
                            {
                                jBase = sequenceMatrix_rowColumn(aSeqMat, j, lPos);
                                kBase = sequenceMatrix_rowColumn(aSeqMat, k, lPos);
                                if (jBase != 'N' && kBase != 'N' && jBase != kBase)
                                {
                                    allsame = 0;
                                    break;
                                }
                            }
                        }
                        if (allsame)
                        {
                            sameCountAnc++;
                        }
                        totalCountAnc++;
                    }
                }
            }
        }
    }
    return (sameCountDer/(float)totalCountDer) / (sameCountAnc/(float)totalCountAnc);
}

/*Get omega at known fixation position.  I.E. center of simulated chromosome*/

double omegaAtCenter(struct sequenceMatrix *aSeqMat, int start, int stop, double site){
	int i,j, segsites, loc_i, loc_j;
	vector *locs;
	double sum,sumL,sumR,comp,denom,denomL,denomR;
	double numer;
	//int temp_site = 0;

	segsites = segSiteCountFromTo(aSeqMat,start,stop);
	sum = sumL = sumR = comp = denom = denomL = denomR = 0;
	if(segsites < 3)
		return 0.0;
		
	locs = vectorNew(1000);
	segSiteLocationsFromTo(aSeqMat,locs,start,stop);
	
	/*Get snp that is nearest to our fixation*/

	/*min = (fabs(site - vectorGetIndex(locs, 0)));
	for (i = 0; i < segsites; ++i){
		temp_site = vectorGetIndex(locs, i);
		if (fabs(site - temp_site) <= min){
			siteIdx = i;
			min = fabs(site-temp_site);
			}
	}
	*/
	for(i=0; i<segsites-1; i++){
		loc_i = vectorGetIndex(locs, i);
		for(j=i+1; j<segsites; j++){
			loc_j = vectorGetIndex(locs, j);
			comp = dij(aSeqMat, loc_i, loc_j);
			if(loc_i < site && loc_j >= site){
				denom += 1;
				sum += comp;
			}
			if(loc_i < site && loc_j < site){
				denomL += 1;
				sumL += comp;
			}
			if(loc_i >= site && loc_j >= site){
				denomR += 1;
				sumR += comp;		
			}
		}
	}
	denom = sum * (1.0/(double)denom);
	numer = 1.0 / ((double)denomL + (double)denomR);
	numer *= (sumL+sumR);

	if (isnan(denom)){
		return 0.0;
	}
	else{
		return(numer/denom);	
	}
	
}

/*omegaMax -- goes through all possible site divisions to maximize omega
// Kim and Nielsen (2003)
*/

double omegaMaxFromTo(struct sequenceMatrix *aSeqMat, int start, int stop){
        int l, i, j, loc_i, loc_j, segsites;
        double max= 0;
        double tmp=0;
        vector *locs;
        double **dijTable;

	segsites = segSiteCountFromTo(aSeqMat,start,stop);
        dijTable = (double **)malloc( sizeof(double)*segsites );
	locs = vectorNew(segsites);
	segSiteLocationsFromTo(aSeqMat,locs,start,stop);

        if(segsites < 3)
                return(0);
        for(i=0; i<segsites-1; i++){
                dijTable[i] = (double *) malloc( sizeof(double)*segsites );
		loc_i = vectorGetIndex(locs, i);
                for(j=i+1; j<segsites; j++){
                        loc_j = vectorGetIndex(locs, j);
                        dijTable[i][j] = dij(aSeqMat, loc_i, loc_j);
                }
        }
        for(l=3;l<segsites-2;l++){
		tmp = omegaFromTo(aSeqMat,start,stop,l,dijTable);
                if(tmp > max){
                        max = tmp;
                }
        }

	for(i=0; i < segsites-1; i++){
                free(dijTable[i]);
        }
        free(dijTable);
	vectorFree(locs);

        return(max);
}

//returns FET PValue between two sites
double siteAssoc(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	int pA, pB, pAB, pa, pab, pAb, paB;
	int comp,i;
	double freqA, freqB;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
	comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}
	
	freqA = (double) pA / comp;
	freqB = (double) pB / comp;
	if(comp * (1.0-freqA) < 1.0 || comp * (1.0-freqB) < 1.0){
		return(-99.0);
	}
	else{
		paB = pB - pAB;
		pAb = pA - pAB;
		pa = comp - pA;
		pab = pa - paB;
		return(FsXctTst(pAB,pAb,paB,pab));
	}
}

//outputPolyMask -- writes mask files named for start position
// for use with sweepCoal2 etc.
void outputPolyMask(struct sequenceMatrix *aSeqMat, int start, int stop, FILE *outfile){
	int i, j, seqNumber, nstart, nstop, flag, action;

	nstart = nstop = action = flag = 0;
	seqNumber = 0;
	for(i=0;i<aSeqMat->sampleSize;i++){
		nstart = nstop = action = flag = 0;
		for(j = start; j < stop; j++){
			if (sequenceMatrix_rowColumn(aSeqMat,i,j) == 'N'){
				if (flag){
					nstop++;
				}
				else{
					nstart = j;
					nstop = j;
					flag = 1;
					if (action == 0){
						action = 1;
					}
				}
			}
			else{
				if (flag && action){
					fprintf(outfile,"%d %lf %lf\n",seqNumber,(float)nstart/aSeqMat->length, (float)nstop/aSeqMat->length);
					flag = 0;
				}
			}
		}
		if (flag && action){
			fprintf(outfile,"%d %lf %lf\n",seqNumber,(float)nstart/aSeqMat->length, (float)nstop/aSeqMat->length);
			flag = 0;
		}
		seqNumber++;
	}
	//finish
	fprintf(outfile,"\n//\n");
}

//outputPolyMaskBed -- same as above but uses start and stop to determine normalization
void outputPolyMaskBed(struct sequenceMatrix *aSeqMat, int start, int stop, FILE *outfile){
	int i, j, seqNumber, nstart, nstop, flag, action;

	nstart = nstop = action = flag = 0;
	seqNumber = 0;
	for(i=0;i<aSeqMat->sampleSize;i++){
		nstart = nstop = action = flag = 0;
		for(j = start; j < stop; j++){
			if (sequenceMatrix_rowColumn(aSeqMat,i,j) == 'N'){
				if (flag){
					nstop++;
				}
				else{
					nstart = j;
					nstop = j+1;
					flag = 1;
					if (action == 0){
						action = 1;
					}
				}
			}
			else{
				if (flag && action){
					fprintf(outfile,"%d %lf %lf\n",seqNumber,((float)(nstart-start))/(stop-start), (float)(nstop-start)/(stop-start));
					flag = 0;
				}
			}
		}
		if (flag && action){
			fprintf(outfile,"%d %lf %lf\n",seqNumber,((float)(nstart-start))/(stop-start), (float)(nstop-start)/(stop-start));
			flag = 0;
		}
		seqNumber++;
	}
	//finish
	fprintf(outfile,"\n//\n");
}

/******************** two seqMat things ************/
double fstFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, struct sequenceMatrix *merge, int start,int stop){
	double pi1, pi2, piTot, piW, f;
	pi1 = nucdivFromTo(aSeqMat,start,stop);
	pi2 = nucdivFromTo(bSeqMat,start,stop);
	if(pi1 == 0 && pi2 == 0){
		return(0.0);
	}
	piW = (pi1 * aSeqMat->sampleSize) + (pi2 * bSeqMat->sampleSize);
	piW /= (double) (aSeqMat->sampleSize + bSeqMat->sampleSize);
	piTot = nucdivFromTo(merge,start,stop);
	//printf("test: %f %f %f %f %d\n",pi1,pi2,piW,piTot,merge->sampleSize);
	f = (piTot - piW) / piTot;
	/*if (f<0)  //Allows the return of negative fsts
		return(0.0);
	else */	
		return(f);
}


//Snn -- Hudson's Snn statistic from 
double SnnFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop){
	double count = 0;
	int i;
	
	for(i=0;i<aSeqMat->sampleSize; i++){
		count += xij_SnnFromTo(aSeqMat,bSeqMat,i,1,start,stop);
	}
	for(i=0;i<bSeqMat->sampleSize; i++){
		count += xij_SnnFromTo(aSeqMat,bSeqMat,i,2,start,stop);
	}
	count /= (double) (aSeqMat->sampleSize + bSeqMat->sampleSize);
	return(count);
}

//xij counts the number of nearest neighbors in same population
double xij_SnnFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, \
	int seqIndex1, int belongFlag, int start,int stop){
	
	int i, n1, n2;
	double minWith = 667.0;
	double minBet  = 667.0;
	double tmp;
	int minCountW = 0;
	int minCountB = 0;
	
	n1 = aSeqMat->sampleSize;
	n2 = bSeqMat->sampleSize;
	if(belongFlag == 1){
		//get within min and count
		for(i = 0; i < n1; i++){
			if(i != seqIndex1){
				tmp = seqDist_SnnFromTo(aSeqMat->matrix[seqIndex1],aSeqMat->matrix[i],start,stop);
				if(tmp >= 0.0){
					if( tmp< minWith){
						minCountW = 1;
						minWith = tmp;
					}
					if( tmp == minWith)
						minCountW += 1;
				}
			}
		}
		//get between min and count
		for(i = 0; i < n2; i++){
			tmp = seqDist_SnnFromTo(aSeqMat->matrix[seqIndex1],bSeqMat->matrix[i],start,stop);
			if(tmp >= 0.0){
				if( tmp< minBet){
					minCountB = 1;
					minBet = tmp;
				}
				if( tmp == minBet)
					minCountB += 1;
			}
		}		
	}
	else{
		//get within min and count
		for(i = 0; i < n2; i++){
			if(i != seqIndex1){
				tmp = seqDist_SnnFromTo(bSeqMat->matrix[seqIndex1],bSeqMat->matrix[i],start,stop);
				if(tmp >= 0.0){
					if( tmp< minWith){
						minCountW = 1;
						minWith = tmp;
					}
					if( tmp == minWith)
						minCountW += 1;
				}
			}
		}
		//get between min and count
		for(i = 0; i < n1; i++){
			tmp = seqDist_SnnFromTo(bSeqMat->matrix[seqIndex1],aSeqMat->matrix[i],start,stop);
			if(tmp >= 0.0){
				if( tmp< minBet){
					minCountB = 1;
					minBet = tmp;
				}
				if( tmp == minBet)
					minCountB += 1;
			}
		}	
	}
	//printf("minWith: %f; minBet: %f; minCountW: %d; minCountB: %d\n", minWith, minBet, minCountW, minCountB);
	if(minWith < minBet){
		return(1.0);
	}
	if(minWith == minBet){
		return(minCountW / (double) (minCountW+minCountB));
	}
	return(0);

}

double seqDist_SnnFromTo( stringWrap *seq1,  stringWrap *seq2, int start, int stop){
	
	int i;
	double compareCount, nCount;
	double count = 0.0;
	char c1, c2;
	
	compareCount = 0.0;
	nCount = 0.0;
	
	for(i = start; i < stop; i++){
		c1 = seq1->cString[i];
		c2 = seq2->cString[i];
		if(c1 == 'N' || c2 == 'N'){
			nCount += 1; //uncertainty about states?
		}
		else{
			if(c1 != c2)
				count += 1;
		}
		compareCount +=1;
	}
	//arbitrary coverage requirement?
	if (nCount / compareCount > 0.75){
		//return large number
		return(-666.0);
	}
	return(count);
}

//Dxy statistic from Takahata and Nei 1985
double DxyFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop){
	int i,j;
	double sum = 0.0;
	double nComps = 0.0;
	double tmp;
	

	for(i=0;i<aSeqMat->sampleSize; i++){
		for(j=0;j<bSeqMat->sampleSize;j++){
			tmp = seqDist_SnnFromTo(aSeqMat->matrix[i],bSeqMat->matrix[j],start,stop);
			if (tmp >= 0.0)
			{
				sum += tmp;
				nComps++;
			}
		}
	}
	return(sum / nComps);
}

double Dxy_minFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop){
	int i,j;
	double min = 666666666.0;
	double tmpVal;
	

	for(i=0;i<aSeqMat->sampleSize; i++){
		for(j=0;j<bSeqMat->sampleSize;j++){
			tmpVal= seqDist_SnnFromTo(aSeqMat->matrix[i],bSeqMat->matrix[j],start,stop);
			if(tmpVal < min && tmpVal >= 0.0) min = tmpVal;
		}
	}
	return(min);
}

void statVecMoments(double *vec, int vecLen, double *mean, double *var, double *skew, double *kurt)
{
        double sum, sum2;
        int i;

        sum = 0;
        for (i=0; i<vecLen; i++)
        {
                sum += vec[i];
        }
        (*mean) = sum/(float) vecLen;
        sum = 0;
        for (i=0; i<vecLen; i++)
        {
                sum += pow(vec[i]-(*mean), 2);
        }
        (*var) = sum/(float) (vecLen-1);
        sum = 0;
        sum2 = 0;
        for (i=0; i<vecLen; i++)
        {
                sum += pow(vec[i]-(*mean), 3);
                sum2 += pow(vec[i]-(*mean), 4);
        }
        (*skew) = (sum/(float) vecLen) / pow((*var), 1.5);
        (*kurt) = (sum2/(float) vecLen) / pow((*var), 2);
}

//cmp_doubles -- for use with qsort
int cmp_doubles(const void *x, const void *y){
        double xx = *(double*)x;
        double yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

void statVecMinMedMax(double *vec, int vecLen, double *minVal, double *medVal, double *maxVal)
{
	double *vec2;
	int i;
	vec2 = malloc(sizeof(double) * vecLen);
	for(i=0;i<vecLen;i++) vec2[i]=vec[i];
        qsort(vec2, vecLen, sizeof(vec2[0]), cmp_doubles);
        (*minVal) = vec2[0];
        (*maxVal) = vec2[vecLen-1];
        if (vecLen % 2 == 1)
        {
                (*medVal) = vec2[vecLen/2];
        }
        else
        {
                (*medVal) = vec2[vecLen/2] + vec2[(vecLen/2)-1];
        }
	free(vec2);
}

//pairwiseDistRankAmongSampleRange takes a number of pairwise differences 
//and returns the number of pairs of samples this alignment that have a lesser or equal number of differences
double pairwiseDistRankAmongSampleRange(struct sequenceMatrix *aSeqMat, int diffs, double *hetVar, int start, int stop){
	int i, j, currDiffs, leCount, numComps, sumDiffs;
	double meanDiffs, numer;
	double *hetLs;
	sumDiffs = leCount = numComps = 0;
	(*hetVar)=0.0;
	
	hetLs = (double *)malloc(sizeof(double)*(aSeqMat->sampleSize*(aSeqMat->sampleSize-1)/2));
	for(i=0;i<aSeqMat->sampleSize-1;i++){
		for(j=i+1;j<aSeqMat->sampleSize;j++){
			currDiffs = seqDist_SnnFromTo(aSeqMat->matrix[i],aSeqMat->matrix[j],start,stop);
			if (currDiffs >= 0.0){
				if (currDiffs <= diffs){
					leCount += 1;
				}
				hetLs[numComps]=currDiffs;
				numComps += 1;
				sumDiffs += currDiffs;
                        }
		}
	}

	numer = 0.0;
	meanDiffs = sumDiffs/(float)numComps;
	for(i=0; i<numComps; i++){
		numer += pow((hetLs[i]-meanDiffs),2);
	}
	(*hetVar) = numer/((float)numComps-1);

	free(hetLs);
	return leCount/ (float) numComps;
}

double *pairwiseIBSVec1PopnFromTo(struct sequenceMatrix *aSeqMat, int *vecLen, int start, int end){
        double tmpLen, *vec;
        (*vecLen)=0;
        int i, j, k, tractStart, snpIndex, segsites;
	vector *locs;

        segsites = segSiteCountFromTo(aSeqMat,start,end);
        locs = vectorNew(segsites);
        segSiteLocationsFromTo(aSeqMat,locs,start,end);
        vec = (double *) malloc(sizeof(double) *(segsites+1)*aSeqMat->sampleSize*(aSeqMat->sampleSize-1)/2);
        for(i=0; i<aSeqMat->sampleSize-1;i++){
                for(j=i+1;j<aSeqMat->sampleSize;j++){
                        tractStart=start;
                        //now iterate across sites
                        for(snpIndex=0;snpIndex<segsites;snpIndex++){
				k = vectorGetIndex(locs, snpIndex);
				if(aSeqMat->matrix[i]->cString[k] != aSeqMat->matrix[j]->cString[k] && aSeqMat->matrix[i]->cString[k] != 'N' && aSeqMat->matrix[j]->cString[k] != 'N'){
					tmpLen = (k-tractStart)/(float)(end-start);
					vec[(*vecLen)++] = tmpLen;
					tractStart = k;
				}
			}
			tmpLen = (end-tractStart)/(float)(end-start);
			vec[(*vecLen)++] = tmpLen;
		}
	}
	vectorFree(locs);
	return vec;
}


double *hetVec1PopnFromTo(struct sequenceMatrix *aSeqMat, int *vecLen, int start, int end)
{
        (*vecLen) = 0;
        int i, j, k;
        double diffs;
        double *vec;

        vec = (double *) malloc(sizeof(double) *aSeqMat->sampleSize*(aSeqMat->sampleSize-1)/2);
        for(i=0; i<aSeqMat->sampleSize-1;i++){
                for(j=i+1;j<aSeqMat->sampleSize;j++){
                        diffs = 0.;
                        for(k=start;k<end;k++){
                                if(aSeqMat->matrix[i]->cString[k] != aSeqMat->matrix[j]->cString[k] && aSeqMat->matrix[i]->cString[k] != 'N' && aSeqMat->matrix[j]->cString[+k] != 'N'){
                                        diffs += 1;
                                }
                        }
                        vec[(*vecLen)++] = diffs/(double)(end-start);
                        //printf("het: %f\n", vec[*(vecLen)-1]);
                }
        }
        return vec;
}

///////////// IBS stuff

double pairwiseIBSMax2PopnFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop){
	int i, j,k, length, snpIndex, segsites;
	double tmpLen, max;
	double startL;
	char c1,c2;
	vector *locs;
	
	length = stop-start;
	max = 0.0;
        segsites = segSiteCountSubpopFromTo(aSeqMat,bSeqMat,start,stop);
        locs = vectorNew(segsites);
        segSiteLocationsSubpopFromTo(aSeqMat,bSeqMat,locs,start,stop);
	for( i=0; i<aSeqMat->sampleSize;i++){
		for(j=0;j<bSeqMat->sampleSize;j++){
			startL=0.0;
			//now iterate across sites
			for(snpIndex=0;snpIndex<segsites;snpIndex++){
                                k = vectorGetIndex(locs, snpIndex);
				c1 = aSeqMat->matrix[i]->cString[k];
				c2 = bSeqMat->matrix[j]->cString[k];
				if(c1 != 'N' && c2 != 'N' && c1 != c2){
					tmpLen = (float) (k-start) / length;
					tmpLen -= startL;
					//printf("%d %d %d %f\n", i, j, k, tmpLen);
					if(tmpLen > max){
						max = tmpLen;
					}
					startL = (float) (k-start) / length;
				}
			}	
			tmpLen = 1.0-startL;
			//printf("%d %d %d %f\n", i, j, k, tmpLen);
			if(tmpLen > max){
				max = tmpLen;
			}
		}
	}
	vectorFree(locs);
	return max;
}

double pairwiseIBSMeanWithinFromTo(struct sequenceMatrix *aSeqMat, int start,int stop){
	int i, j,k, length, snpIndex, segsites;
	double tmpLenSum,tmpLen;
	float comp=0.0;
	double startL;
	char c1,c2;
	vector *locs;
	
	tmpLenSum =0.0;
	length = stop-start;
        segsites = segSiteCountFromTo(aSeqMat,start,stop);
        locs = vectorNew(segsites);
        segSiteLocationsFromTo(aSeqMat,locs,start,stop);
	for( i=0; i<aSeqMat->sampleSize-1;i++){
		for(j=i+1;j<aSeqMat->sampleSize;j++){
			startL=0.0;
			//now iterate across sites
			for(snpIndex=0;snpIndex<segsites;snpIndex++){
                                k = vectorGetIndex(locs, snpIndex);
				c1 = aSeqMat->matrix[i]->cString[k];
				c2 = aSeqMat->matrix[j]->cString[k];
				if(c1 != 'N' && c2 != 'N' && c1 != c2){
					tmpLen = (float) (k-start) / length;
					tmpLen -= startL;
					tmpLenSum+=tmpLen;
					comp+=1.0;
					startL = (float) (k-start) / length;
				}
			}
			tmpLenSum += 1.0-startL;
			comp += 1.0;	
		}
	}
	vectorFree(locs);
	return tmpLenSum/comp;
}


///////////
////
///
/// still need to convert all below to stringWrap paradigm
///

////////////////////
/*   thetah - pi   */

/*
double hfay( int nsam, int segsites, char **list){
  int s, frequency( char, int, int, char**);
  double pi, p1, nd, nnm1  ;
  
  pi = 0.0 ;
  nd = nsam;
  nnm1 = nd/(nd-1.0) ;
  for( s = 0; s <segsites; s++){
    p1 = frequency('1', s,nsam,list)/nd ;
    pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ;
  }
  return( -pi ) ;
}

*/


/*fixed diffs  */
int frequencyFD( char allele,int site,int nsam,  char **list){
  int i, count=0;
  for( i=1; i<nsam; i++){
    count += ( list[i][site] == allele ? 1: 0 );
  }
  return( count);
} 

int fixedDiffs(int segsites, int nsam, char **list){
  int i, fd=0;
  char allele;

  for(i=0; i < segsites; i++){
    allele = list[0][i];
    fd += ((frequencyFD(allele, i, nsam, list) == 0) ? 1:0);
  }
  return(fd);
}

int frequencyFDAllele1( char allele,int site,int nsam,  char **list){
  int count=0;
  count += ( list[1][site] == allele ? 1: 0 ) ;
  return( count);
} 

int fixedDiffsAllele1(int segsites, int nsam, char **list){
  int i, fd=0;
  char allele;

  for(i=0; i < segsites; i++){
    allele = list[0][i];
    fd += (frequencyFDAllele1(allele, i, nsam, list) == 0 ? 1:0);
  }
  return(fd);
}

int ingroupSegSites(int segsites, int nsam, char **list){
  int i, ss = 0;
  
  for(i=0; i < segsites; i++){
    ss += (((frequencyFD('1', i, nsam, list) > 0) && (frequencyFD('1', i, nsam, list) < nsam - 1 )) ? 1:0);
  }
  return(ss);
}

double nucdivIn( int nsam, int segsites, char **list){
  int s;
  double pi, p1, nd, nnm1  ;

  pi = 0.0 ;

  nd = nsam - 1;
  nnm1 = nd/(nd-1.0) ;
  for( s = 0; s <segsites; s++){
    p1 = frequencyFD('1', s,nsam,list)/nd ;
    pi += 2.0*p1*(1.0 -p1)*nnm1 ;
  }
  return( pi ) ;
}

/*haploCount-- returns number of unique haplotypes.
  WARNGING-- sorts seqs in place */
int haploCount(char **list, int nsam){
  int i, count = 1;
  char *tmpSeq;

  
  qsort(list, nsam, sizeof(char *), cmpr); 
  tmpSeq = list[0];
  for(i = 0; i < nsam; i++){
    if (strcmp(list[i],tmpSeq)){
      count += 1;
      tmpSeq = list[i];
    }
  }
  return(count);
}

int cmpr(const void *a, const void *b) { 
 return strcmp(*(char **)a, *(char **)b);
}





	
/************************* tajima.c *************************************************************
 This program calculates Tajima's D when number of sequences, number of segregating sites,
   and average pairwise differences (pi) are known.  It also reports all the coefficients for Tajima's
   D (a1, a2, b1, b2, c1, c2, e1, e2). 
**************************************************************************************************/
double tajd(int nsam, int segsites, double sumk){

  double  a1, a2, b1, b2, c1, c2, e1, e2; 
  
  if( segsites == 0 ) return( 0.0) ;
  
  a1 = a1f(nsam);
  a2 = a2f(nsam);
  b1 = b1f(nsam);
  b2 = b2f(nsam);
  c1 = c1f(a1, b1);
  c2 = c2f(nsam, a1, a2, b2);
  e1 = e1f(a1, c1);
  e2 = e2f(a1, a2, c2);

  return( (sumk - (segsites/a1))/sqrt((e1*segsites) + ((e2*segsites)*(segsites
								      -1))) ) ;
}

double a1f(int nsam){
  double a1;
  int i;
  a1 = 0.0;
  for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
  return (a1);
}


double a2f(int nsam) {
  double a2;
  int i;
  a2 = 0.0;
  for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
  return (a2);
}


double b1f(int nsam){
  double b1;
  b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
  return (b1);
}


double b2f(int nsam){
  double b2;
  b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
  return (b2);
}


double e1f(double a1, double c1){
  double e1;
  e1 = c1/a1;
  return (e1);
}

double e2f(double a1, double a2, double c2){ 
  double e2;
  e2 = c2/((a1*a1)+a2);
  return (e2);
}

double c1f(double a1, double b1){
  double c1;
  c1 = b1 - (1/a1);
  return (c1);
}

double c2f(int nsam, double a1, double a2, double b2){
  double c2;
  c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
  return (c2);
}


//very general ptr swap
void swap( void **p1,  void **p2)
{
  void *pt = *p1;
  *p1 = *p2;
  *p2 = pt;
}

//for sorting stuff
int compare_doubles(const void *a,const void *b){
	double *pa = (double *) a;
	double *pb = (double *) b;
	if ((*pa - *pb) > 0 ){
		return 1;
	}
	else{
		return - 1;
	}
}

