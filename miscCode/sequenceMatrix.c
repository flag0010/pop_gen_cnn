/* sequenceMatrix.c -- object for doing popGen etc. */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "sequenceMatrix.h"
#include "../pgSummaryStats.h"


struct sequenceMatrix *sequenceMatrix_new(int sampleSize, int length){
	struct sequenceMatrix *aSeqMat;
	stringWrap **aMat, **nameVect;
	int i;
	
	aSeqMat = malloc(sizeof(sequenceMatrix));
	assert(aSeqMat);
	aMat = malloc(sampleSize * sizeof(stringWrap *));
	assert(aMat);
	nameVect = malloc(sampleSize * sizeof(stringWrap *));
	assert(nameVect);
	for(i=0; i < sampleSize; i++){
		aMat[i] = NULL;//stringWrapNew(length);	
		nameVect[i] = NULL; //stringWrapNew(100);		
	}
	aSeqMat->matrix = aMat;
	aSeqMat->names = nameVect;
	aSeqMat->sampleSize = sampleSize;
	aSeqMat->length = length;
	
	return(aSeqMat);
}

void sequenceMatrix_newByRef(int sampleSize, int length, struct sequenceMatrix *aSeqMat){
	stringWrap **aMat, **nameVect;
	int i;
	
	assert(aSeqMat);
	aMat = malloc(sampleSize * sizeof(stringWrap *));
	assert(aMat);
	nameVect = malloc(sampleSize * sizeof(stringWrap *));
	assert(nameVect);
	for(i=0; i < sampleSize; i++){
		aMat[i] = NULL;//stringWrapNew(length);	
		nameVect[i] = NULL; //stringWrapNew(100);		
	}
	aSeqMat->matrix = aMat;
	aSeqMat->names = nameVect;
	aSeqMat->sampleSize = sampleSize;
	aSeqMat->length = length;
}

void sequenceMatrix_free(struct sequenceMatrix *aSeqMat){
	int i;
	for(i=0; i < aSeqMat->sampleSize; i++){
		stringWrapFree(aSeqMat->matrix[i]);
		stringWrapFree(aSeqMat->names[i]);
	}
	free(aSeqMat);
}

/*fasta import */
struct sequenceMatrix *sequenceMatrix_importFasta(char *fileName){
	int initSeqs, initLength, count;
	FILE *infile;
	struct sequenceMatrix *newSeqMat;
	stringWrap *line;
	
	//open file
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	//initialization
	initSeqs = 400;
	initLength = 30000000;
	count = -1;
	newSeqMat = sequenceMatrix_new(initSeqs, initLength);
	//go line by line through file
	line = stringWrapNew(SMALL);
	while(stringWrapReadLine(line,infile) != EOF){

		if(line->cString[0] == '>'){
			newSeqMat->names[++count] = stringWrapSlice(line,1,stringWrapLength(line) - 1);
			newSeqMat->matrix[count] = stringWrapNew(0);		
		}
		else{
			stringWrapChomp(line);
			stringWrapAppend(newSeqMat->matrix[count],line);
		}
	}
	newSeqMat->sampleSize = count + 1;
	newSeqMat->length = newSeqMat->matrix[0]->length;
	return(newSeqMat);
}

/* ///////
//
/ Utilities
/
///////// */

//Ns out any alignment column with too much missing data
void sequenceMatrix_NOutSitesWithTooMuchMissingData(struct sequenceMatrix *aSeqMat, double maxFracMissingData){
	int i, j, numberOfNs;
	for(j=0;j<aSeqMat->length;j++){
		numberOfNs = 0;
	        for(i=0;i<aSeqMat->sampleSize;i++){
			if (sequenceMatrix_rowColumn(aSeqMat,i,j) == 'N'){
				numberOfNs++;
			}
		}
		if ((double) numberOfNs / (double) aSeqMat->sampleSize > maxFracMissingData) {
	                for(i=0;i<aSeqMat->sampleSize;i++){
				aSeqMat->matrix[i]->cString[j] = 'N';
			}
		}
                //else{
		//	printf("keeping position %d; %d of %d Ns\n", j+1,numberOfNs,aSeqMat->sampleSize);
		//}
	}
}

void sequenceMatrix_NOutSitesWithTooMuchMissingDataBothMats(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, double maxFracMissingData){
	int i, j, numberOfNs;
	for(j=0;j<aSeqMat->length;j++){
		numberOfNs = 0;
	        for(i=0;i<aSeqMat->sampleSize;i++){
			if (sequenceMatrix_rowColumn(aSeqMat,i,j) == 'N'){
				numberOfNs++;
			}
		}
	        for(i=0;i<bSeqMat->sampleSize;i++){
			if (sequenceMatrix_rowColumn(bSeqMat,i,j) == 'N'){
				numberOfNs++;
			}
		}
		if ((double) numberOfNs / (double) (aSeqMat->sampleSize + bSeqMat->sampleSize) > maxFracMissingData) {
	                for(i=0;i<aSeqMat->sampleSize;i++){
				aSeqMat->matrix[i]->cString[j] = 'N';
			}
			for(i=0;i<bSeqMat->sampleSize;i++){
				bSeqMat->matrix[i]->cString[j] = 'N';
			}
		}
                //else{
		//	printf("keeping position %d; %d of %d Ns\n", j+1,numberOfNs,aSeqMat->sampleSize);
		//}
	}
}

//siteArray -- places alignment column into pre-alloced stringWrap
void sequenceMatrix_siteArray(struct sequenceMatrix *aSeqMat, stringWrap *dst, int index){
	int i;
	stringWrapInit(dst);
	for(i=0;i<aSeqMat->sampleSize;i++){
		stringWrapAppendSingle(dst, aSeqMat->matrix[i]->cString[index]);
	}
}

//siteArrayClean -- cleans out Ns
void sequenceMatrix_siteArrayClean(struct sequenceMatrix *aSeqMat, stringWrap *dst, int index){
	int i;
	char test = '\0';
	
	stringWrapInit(dst);
	for(i=0;i<aSeqMat->sampleSize;i++){
		test = sequenceMatrix_rowColumn(aSeqMat,i,index);
		if( test != 'N')
			stringWrapAppendSingle(dst, test);
	}
}

//siteSet is the uniq characters at an index
void sequenceMatrix_siteSet(struct sequenceMatrix *aSeqMat, stringWrap *dst, int index){
	int i;
	
	stringWrapInit(dst);
	for(i=0;i<aSeqMat->sampleSize;i++){
		if(stringWrapHasChar(dst,aSeqMat->matrix[i]->cString[index]) == 0){
			stringWrapAppendSingle(dst, aSeqMat->matrix[i]->cString[index]);
		}
	}
}
//clean version of above
void sequenceMatrix_siteSetClean(struct sequenceMatrix *aSeqMat, stringWrap *dst, int index){
	int i;
	stringWrapInit(dst);
	for(i=0;i<aSeqMat->sampleSize;i++){
		if(stringWrapHasChar(dst,aSeqMat->matrix[i]->cString[index]) == 0 && aSeqMat->matrix[i]->cString[index] != 'N'){
			stringWrapAppendSingle(dst, aSeqMat->matrix[i]->cString[index]);
		}
	}
}

void sequenceMatrix_siteSetSubpopClean(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, stringWrap *dst, int index){
	int i;
	stringWrapInit(dst);
	for(i=0;i<aSeqMat->sampleSize;i++){
		if(stringWrapHasChar(dst,aSeqMat->matrix[i]->cString[index]) == 0 && aSeqMat->matrix[i]->cString[index] != 'N'){
			stringWrapAppendSingle(dst, aSeqMat->matrix[i]->cString[index]);
		}
	}
	for(i=0;i<bSeqMat->sampleSize;i++){
		if(stringWrapHasChar(dst,bSeqMat->matrix[i]->cString[index]) == 0 && bSeqMat->matrix[i]->cString[index] != 'N'){
			stringWrapAppendSingle(dst, bSeqMat->matrix[i]->cString[index]);
		}
	}
}

//siteArray -- places subset of sequence from seq_i into pre-alloced stringWrap
void sequenceMatrix_sitesFromSequence(struct sequenceMatrix *aSeqMat, stringWrap *dst, int alleleNumberI, int index, int length){
	int i, start;
	start = index;
	stringWrapInit(dst);
	for(i=0;i<length;i++){
		stringWrapAppendSingle(dst, aSeqMat->matrix[alleleNumberI]->cString[start++]);
	}
}


//simplified access

char sequenceMatrix_rowColumn(struct sequenceMatrix *aSeqMat,int r, int c){
	return(aSeqMat->matrix[r]->cString[c]);
}

void sequenceMatrix_rowColumn_set(struct sequenceMatrix *aSeqMat,int r, int c, char x){
	aSeqMat->matrix[r]->cString[c] = x;
}

void sequenceMatrix_printFasta(struct sequenceMatrix *aSeqMat){
	int i,j,count;
	for(i=0;i<aSeqMat->sampleSize;i++){
		printf(">%s\n",aSeqMat->names[i]->cString);
		count = 60;
		for(j=0;j<aSeqMat->length;j++){
			printf("%c",sequenceMatrix_rowColumn(aSeqMat,i,j));
			count -= 1;
			if(count == 0){
				printf("\n");
				count = 60;
			}
		}
		printf("\n");
	}
}

//merge two into one
struct sequenceMatrix *sequenceMatrix_merge(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat){
	struct sequenceMatrix *new;
	int n,i;
	n = aSeqMat->sampleSize + bSeqMat->sampleSize;
	new = sequenceMatrix_new(n,aSeqMat->length);
	for(i=0;i<aSeqMat->sampleSize;i++){
		new->matrix[i] = aSeqMat->matrix[i];
		new->names[i] = aSeqMat->names[i];
	}
	n=aSeqMat->sampleSize;
	for(i=0;i<bSeqMat->sampleSize;i++){
		new->matrix[i+n] = bSeqMat->matrix[i];
		new->names[i+n] = bSeqMat->names[i];
	}
	return(new);
}

//merge three into one
struct sequenceMatrix *sequenceMatrix_merge3(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat,struct sequenceMatrix *cSeqMat){
	struct sequenceMatrix *new;
	int n,i,count;
	n = aSeqMat->sampleSize + bSeqMat->sampleSize + cSeqMat->sampleSize;
	count = 0;
	new = sequenceMatrix_new(n,aSeqMat->length);
	for(i=0;i<aSeqMat->sampleSize;i++){
		new->matrix[count] = aSeqMat->matrix[i];
		new->names[count++] = aSeqMat->names[i];
	}
	for(i=0;i<bSeqMat->sampleSize;i++){
		new->matrix[count] = bSeqMat->matrix[i];
		new->names[count++] = bSeqMat->names[i];
	}
	for(i=0;i<cSeqMat->sampleSize;i++){
		new->matrix[count] = cSeqMat->matrix[i];
		new->names[count++] = cSeqMat->names[i];
	}
	return(new);
}


