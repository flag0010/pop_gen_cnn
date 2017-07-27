//pgStats for a bedFile rather than windows
// calculates fst and related

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "miscCode/stringWrap.h"
#include "miscCode/sequenceMatrix.h"
#include "pgSummaryStats.h"
#include "miscCode/bedFile.h"


void usage();

int main(int argc, char *argv[]){
	struct sequenceMatrix *aSeqMat, *ancSeqMat, *bSeqMat;
	int i, bedElNumber, h, ss, private1, private2;
	long int numSites;
	struct bedEl data[80000];
	double hetVar1, hetVar2;
	double pi1,pi2, theta_h, z1,z2,ztot,f,snn,tajD, H;
	double gmin,dxy_mean,dxy_min,dd1,dd2,ddRank1,ddRank2;
	double ibsMaxBetween, ibsMaxWithin1, ibsMaxWithin2;
	float maxFracMissingData;
	struct sequenceMatrix *merge;
	
	if(argc < 6){
		usage();
		exit(1);
	}

	//open fastaFile and bedFile
	aSeqMat = sequenceMatrix_importFasta(argv[1]);
	bSeqMat = sequenceMatrix_importFasta(argv[2]);
	ancSeqMat = sequenceMatrix_importFasta(argv[3]);
	bedElNumber = bedFileImport3(argv[4], data);
	maxFracMissingData = atof(argv[5]);
	
	sequenceMatrix_NOutSitesWithTooMuchMissingDataBothMats(aSeqMat,bSeqMat,maxFracMissingData);
	merge = sequenceMatrix_merge(aSeqMat, bSeqMat);


	//print header
	printf("chrom\tchromStart\tchromEnd\tnumSites\tpi1\thetVar1\tss1\tprivate1\tthetaH1\ttajd1\tH1\tHapCount1\tZnS1\t");
	printf("pi2\thetVar2\tss2\tprivate2\tthetaH2\ttajd2\tH2\tHapCount2\tZnS2\t");
	printf("Fst\tsnn\tdxy_mean\tdxy_min\tgmin\tzx\tdd1\tdd2\tddRank1\tddRank2\tibsMaxB\tibsMean1\tibsMean2\n");
	

	//loop through beds; the adjustments to end are to honor the zero indexed half open bed convention
	for(i=0;i<bedElNumber;i++){
		numSites = numColumnsNotNedOutFromToBothMats(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd);
		dxy_min = Dxy_minFromTo(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ddRank1 = pairwiseDistRankAmongSampleRange(aSeqMat, dxy_min, &hetVar1, data[i].chromStart, data[i].chromEnd+1);
		ddRank2 = pairwiseDistRankAmongSampleRange(bSeqMat, dxy_min, &hetVar2, data[i].chromStart, data[i].chromEnd+1);
		printf("%s\t%ld\t%ld\t%ld\t",data[i].chrom,data[i].chromStart,data[i].chromEnd, numSites);
		//ingroup 1 stats
		pi1 = nucdivFromTo(aSeqMat, data[i].chromStart, data[i].chromEnd+1);
		privateSegSitesInTwoPopnsFromTo(aSeqMat, bSeqMat, &private1, &private2, data[i].chromStart, data[i].chromEnd+1);
		ss = segSiteCountFromTo(aSeqMat, data[i].chromStart, data[i].chromEnd+1);
		theta_h = thetaHFromTo(aSeqMat,ancSeqMat,data[i].chromStart, data[i].chromEnd+1);
		h = nHaplotypes(aSeqMat,data[i].chromStart, data[i].chromEnd+1);
		z1 = ZnSFromTo(aSeqMat,data[i].chromStart, data[i].chromEnd+1);
		tajD = tajd(aSeqMat->sampleSize,ss,pi1);
		H = theta_h-pi1;
		printf("%f\t%f\t%d\t%d\t%f\t%f\t%f\t%d\t%f\t",pi1,hetVar1,ss,private1,theta_h,tajD,H,h,z1);
		
		
		//ingroup 2 stats	
		pi2 = nucdivFromTo(bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ss = segSiteCountFromTo(bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		theta_h = thetaHFromTo(bSeqMat,ancSeqMat,data[i].chromStart, data[i].chromEnd+1);
		h = nHaplotypes(bSeqMat,data[i].chromStart, data[i].chromEnd+1);
		z2 = ZnSFromTo(bSeqMat,data[i].chromStart, data[i].chromEnd+1);
		tajD = tajd(bSeqMat->sampleSize,ss,pi2);
		H = theta_h-pi2;
		printf("%f\t%f\t%d\t%d\t%f\t%f\t%f\t%d\t%f\t",pi2,hetVar2,ss,private2,theta_h,tajD,H,h,z2);
		
		f = fstFromTo(aSeqMat,bSeqMat,merge,data[i].chromStart, data[i].chromEnd+1);
		snn = SnnFromTo(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		dxy_mean = DxyFromTo(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		gmin = dxy_min / dxy_mean; 
		dd1 = dxy_min / pi1;
		dd2 = dxy_min / pi2;
		
		ztot = ZnSFromTo(bSeqMat,data[i].chromStart, data[i].chromEnd+1);
		printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t",f,snn,dxy_mean,dxy_min,gmin,(z1+z2)/2.0/ztot,dd1,dd2,ddRank1,ddRank2) ;

		ibsMaxBetween =pairwiseIBSMax2PopnFromTo(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ibsMaxWithin1= pairwiseIBSMeanWithinFromTo( aSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ibsMaxWithin2= pairwiseIBSMeanWithinFromTo( bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		printf("%f\t%f\t%f",ibsMaxBetween,ibsMaxWithin1,ibsMaxWithin2);
		
		printf("\n");
	}
	sequenceMatrix_free(aSeqMat);
	sequenceMatrix_free(bSeqMat);
	free(merge); //individual pieces freed above
	return(0);
}	

void usage(){
	printf("pgStatsBedSubpop_forML ingroupFastaFile secondPopulationFasta ancestorFastaFile bedFile maxFractionMissingData\n");
}
