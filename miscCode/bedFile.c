/* bedFile utilities
/
/
/ Andrew Kern
*/

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "assert.h"
#include "bedFile.h"

struct bedFile *bedFile_new(){
	struct bedFile *newBed;
	
	newBed = malloc(sizeof(bedFile));
	if(newBed == NULL){
		fprintf(stderr,"didn't make bedFile malloc\n");
		exit(1);
	}
	newBed->next = NULL;

	return(newBed);
}

void bedFile_free(struct bedFile *aBedFile){
	free(aBedFile->data);
	aBedFile->next = NULL;
	free(aBedFile);
}

//not sure if I'll ever use the multiple format, but what the heck -- this adds bedFiles as a linked list
void addBedFile(struct bedFile *parent,struct bedFile *new){
	parent->next = new;
}

//counts number of separate bed samples in bedFile list
int bedFileCount(struct bedFile *aBedFile){
	int i = 0;
	struct bedFile *curPtr;
	
	curPtr = aBedFile;
	while(curPtr != NULL){
		i++;
		curPtr = curPtr->next;
	}
	return(i);
}

//inRangeBedEl asks if a site is >= start < end
int inRangeBedEl(struct bedEl *data, long int site){
	if( (site >= data->chromStart) && (site < data->chromStart)){
		return(1);
	}
	else{
		return(0);
	}
}


//bedFileImport-- reads a file and stores info into pre-alloc'd data for 5 column
//returns bedElNumber
int bedFileImport5(char *fileName, struct bedEl *data){
	FILE *infile;
	long int chromStart, chromEnd;
	int score, j;
	char chrom[1001], name[1001];

	/* open file, errors? */
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	/* go through infile and get bed info*/
	j = 0;
	while (fscanf(infile, "%s %ld %ld %s %d %*s", chrom, &chromStart, &chromEnd, name, &score) != EOF){
		data[j].chrom = chrom;
		data[j].chromStart = chromStart;
		data[j].chromEnd = chromEnd;
		data[j].name = name;
		data[j].score = score;
		j += 1;
	}
	fclose(infile);
	return(j);
}

//bedFileImport-- reads a file and stores info into pre-alloc'd data for 3 column
//returns bedElNumber
int bedFileImport3(char *fileName, struct bedEl *data){
	FILE *infile;
	long int chromStart, chromEnd;
	int  j;
	char chrom[10001];

	/* open file, errors? */
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	/* go through infile and get bed info*/
	j = 0;
	chromEnd = 0;
	while (fscanf(infile, "%s %ld %ld ", chrom, &chromStart, &chromEnd) != EOF){
		data[j].chrom = chrom;
		data[j].chromStart = chromStart;
		data[j].chromEnd = chromEnd;
		j += 1;
	}
	fclose(infile);
	return(j);
}
