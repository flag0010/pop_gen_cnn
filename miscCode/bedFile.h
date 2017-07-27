/* bedFile.h -- for dealing with the bed file format 
/
/
/
/ Also contains the definition of a bedElement */


#ifndef BED_INC
#define BED_INC


#define MAXBEDELS 10000000


 /*bedEl definition- can currently deal with BED9 format */
struct bedEl{
  	long int chromStart, chromEnd, thickStart, thickEnd;
	int score; 
	char *chrom, *name, strand;
}bedEl;

struct bedFile{
	int bedElNumber;
	struct bedEl data[MAXBEDELS];
	struct bedFile *next;
}bedFile;


struct bedFile *bedFile_new();
void bedFile_free(struct bedFile *aBedFile);
void addBedFile(struct bedFile *parent,struct bedFile *new);
int bedFileCount(struct bedFile *aBedFile);
int bedFileImport3(char *fileName, struct bedEl *data);
int bedFileImport5(char *fileName, struct bedEl *data);
int inRangeBedEl(struct bedEl *data, long int site);
#endif
