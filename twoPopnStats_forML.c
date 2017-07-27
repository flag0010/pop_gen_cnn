/******* maskedStats.c ********
for calculating sample stats from MS output 
after it has been filtered by msMask
********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"

#define LINEBUF 1000000
int maxsites = 100000 ;
void usage();

int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i,  howmany  ;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int   segsites, count  , n1, n2, iss,h, private1, private2 ;
	double pi1, pi2 , th,  z1,z2,ztot, f, snn, dxy_min, dxy_mean, H, tajD, gmin,dd1,dd2,ddRank1,ddRank2,ibsMaxBetween,ibsMaxWithin1,ibsMaxWithin2,hetVar1,hetVar2;
	char dum[20], astr[100] ;
	int starFlag=0;
	int migOption=0;
//	int bins = 10;
//	double hist[bins];

	if( argc > 1 ) { 
		n1 = atoi( argv[1] ) ; 
		n2 = atoi( argv[2] ) ; 
		if(argc==4){
			if(argv[3][1] == 'c') migOption=1;
		}
	}
	else{
		usage();
	}

/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, LINEBUF, pfin);


	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;

	//print header line
	printf("pi1\thetVar1\tss1\tprivate1\tthetaH1\ttajd1\tH1\tHapCount1\tZnS1\tpi2\thetVar2\tss2\tprivate2\tthetaH2\ttajd2\tH2\tHapCount2\tZnS2\tFst\tsnn\tdxy_mean\tdxy_min\tgmin\tzx\tdd1\tdd2\tddRank1\tddRank2\tibsMaxB\tibsMean1\tibsMean2" ) ;
//	for(i=0;i<bins;i++){
//		printf("\tibs[%d]",i);
//	}
	printf("\n");
	
	
	while( howmany-count++ ) {

/* read in a sample */
		do {
			fgets( line, LINEBUF, pfin);
		}while ( line[0] != '/' );
		if(line[2] == '*'){
			starFlag = 1;
		}
		else{ 
			starFlag = 0;
		}
		if(migOption==0)starFlag=1;
		
		fscanf(pfin,"  segsites: %d", &segsites );
		if( segsites >= maxsites){
			maxsites = segsites + 10 ;
			posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
			biggerlist(nsam,maxsites, list) ;
		}
		if( segsites > 0) {
			fscanf(pfin," %s", astr);

			for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
			for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
		}
		/* analyse sample ( do stuff with segsites and list) */

		//between popn stats
		snn = Snn(segsites,nsam,n1,n2,list);
		dxy_min = Dxy_min(segsites,nsam,n1,n2,list);
		dxy_mean = Dxy_mean(segsites,nsam,n1,n2,list);
		f = fst2Subs(segsites,nsam,0,n1,n1,n1+n2,list);
		gmin = dxy_min / dxy_mean; 
		ztot = ZnSSub( segsites,  nsam, 0,nsam, list);
		
		//gotta do single pop pi here before the dds
		pi1 = nucdivSub(nsam,segsites,0,n1,list);
		pi2 = nucdivSub(nsam,segsites,n1,nsam,list);
		dd1 = dxy_min / pi1;
		dd2 = dxy_min / pi2;
		ddRank1 = pairwiseDistRankAmongSampleRange( segsites, dxy_min, 0, n1, &hetVar1, list);
		ddRank2 = pairwiseDistRankAmongSampleRange( segsites, dxy_min, n1, n2, &hetVar2, list);

		//single popn stats
                privateSegSitesInTwoPopns(segsites, nsam, n1, &private1, &private2, list);
		iss = segSitesSub(segsites,nsam,0,n1,list);
		th = thetahSub(nsam, segsites,0,n1, list) ;
		h = nHaplotypesSub(segsites,nsam,0,n1,list);
		z1 = ZnSSub( segsites,  nsam, 0,n1, list);
		H = th-pi1;
		tajD = tajd(nsam,iss,pi1);
		if(starFlag==1)
			printf("%lf\t%lf\t%d\t%d\t%lf\t%f\t%f\t%d\t%f\t",pi1, hetVar1, iss, private1, th,tajD,H , h, z1) ;
		iss = segSitesSub(segsites,nsam,n1,nsam,list);
		th = thetahSub(nsam, segsites,n1,nsam, list) ;
		h = nHaplotypesSub(segsites,nsam,n1,nsam,list);
		z2 = ZnSSub( segsites,  nsam, n1,nsam, list);
		H = th-pi2;
		tajD = tajd(nsam,iss,pi2);
		if(starFlag==1)
			printf("%lf\t%lf\t%d\t%d\t%lf\t%f\t%f\t%d\t%f\t",pi2, hetVar2, iss, private2, th,tajD,H , h, z2) ;

		//now print those two popn stats
		if(starFlag==1)
			printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",f,snn,dxy_mean,dxy_min,gmin,(z1+z2)/2.0/ztot,dd1,dd2,ddRank1,ddRank2) ;

		if(starFlag==1){
			ibsMaxBetween = pairwiseIBSMax2Popn(segsites,nsam,n1,posit,list);
			ibsMaxWithin1 = pairwiseIBSMeanWithin(segsites,0, n1, posit, list);
			ibsMaxWithin2 = pairwiseIBSMeanWithin(segsites,n1, nsam, posit, list);
			printf("\t%f\t%f\t%f",ibsMaxBetween,ibsMaxWithin1,ibsMaxWithin2);
		}
//		//IBS tracts
//		if(starFlag==1){
//				pairwiseIBSHist2Popn(segsites,nsam, n1, hist, bins, posit, list);
//				for(i=0;i<bins;i++){
//					printf("\t%f",hist[i]);
//				}
//		}
		printf("\n");
	}
	return(0);
}

void usage(){
	printf("maskedStatsSubpop n1 n2\n");
	printf("returns analysis of Hudson style output assuming two subpops of size n1 and n2\n");
	printf("options:\n\t-c <condition on migration>\n");
	exit(1);
}

