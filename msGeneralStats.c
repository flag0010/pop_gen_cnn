#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"
/* allocates space for gametes (character strings) */


void cmatrix_free(char **m, int nsam)
{
	int i;
	for(i=0;i<nsam;i++) free(m[i]);
	free(m);
}

void imatrix_free(int **m, int nsam)
{
	int i;
	for(i=0;i<nsam;i++) free(m[i]);
	free(m);
}

char **
	cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) )
		perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) )))
			perror("alloc error in cmatric. 2");
	}
	return( m );
}

int **
	imatrix(nsam,len)
	int nsam, len;
{
	int i,ii;
	int **m;

	if( ! ( m = (int **) malloc( (unsigned)( nsam*sizeof( int* )) ) ) )
		perror("alloc error in imatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (int *) malloc( (unsigned) (len*sizeof( int )) )))
			perror("alloc error in imatric. 2");
		for(ii=0;ii<len;ii++)m[i][ii]=0;
	}
	return( m );
}
int
	biggerlist(nsam, nmax, list )
	int nsam ;
unsigned nmax ;
char ** list ;
{
	int i;

	for( i=0; i<nsam; i++){
		list[i] = (char *)realloc( list[i],nmax*sizeof(char) ) ;
		if( list[i] == NULL ) perror( "realloc error. bigger");
	}
	return(0);
}                        


//below we've corrected the sampleSizes for missing data
double nucdiv( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;


	for( s = 0; s <segsites; s++){
		nd = sampleSizeSite(s,nsam,list);
		if (nd > 1){
			nnm1 = nd/(nd-1.0) ;
			p1 = frequency('1', s,nsam,list)/nd ;
			pi += 2.0*p1*(1.0 -p1)*nnm1 ;
		}
	}
	return( pi ) ;
}

double maxFDA( int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        double mfda, p1, nd;

        mfda = 0.0 ;


        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        p1 = frequency('1', s,nsam,list)/nd ;
                        if (p1 < 1.0 && p1 > mfda){
                            mfda = p1;
                        }
                }
        }
        return( mfda ) ;
}

//fills a vector size l with values of nucdiv in "windows"
void nucdivWindow( int nwins, double *posit, double *output, int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	int wcount = 0;
	double pi, p1, nd, nnm1  ;
	double start, end, delta;
	start = 0;
	delta = 1.0 / nwins;
	end = delta;
	
	while(start < 1.0){
		pi = 0.0 ;
		for( s = 0; s <segsites; s++){
			if(posit[s] <=end && posit[s] > start){
				nd = sampleSizeSite(s,nsam,list);
				if (nd > 1){
					nnm1 = nd/(nd-1.0) ;
					p1 = frequency('1', s,nsam,list)/nd ;
					pi += 2.0*p1*(1.0 -p1)*nnm1 ;
				}
			}
		}
		output[wcount++]=pi;
		start += delta;
		end += delta;
	}
}


//fills a vector size l with values of Tajima's D in "windows"
void tajdWindow(int nwins, double *posit, double *output, int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        int wcount = 0;
        double pi, p1, nd, nnm1  ;
        double start, end, delta;
        int segsitesInWin;
        start = 0;
        delta = 1.0 / nwins;
        end = delta;

        while(start < 1.0){
                pi = 0.0 ;
                segsitesInWin = 0;
                for( s = 0; s <segsites; s++){
                        if(posit[s] <=end && posit[s] > start){
                                nd = sampleSizeSite(s,nsam,list);
                                if (nd > 1){
                                        segsitesInWin += 1;
                                        nnm1 = nd/(nd-1.0) ;
                                        p1 = frequency('1', s,nsam,list)/nd ;
                                        pi += 2.0*p1*(1.0 -p1)*nnm1 ;
                                }
                        }
                }
                output[wcount++]=tajd(nsam,segsitesInWin,pi);
                start += delta;
                end += delta;
        }
}


double achazThetaExponentWeights(int nsam, int segsites, char **list, int exponent)
{
        int s, frequency( char, int, int, char**);
        double thetaA, i, wi, nd, nnm1, wSum ;

        wSum = 0.0 ;
        for ( i = 1; i < nsam; i++){
            wi = pow(i,exponent);
            wSum += wi;
        }

        thetaA = 0.0 ;
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        nnm1 = nd/(nd-1.0);
                        i = frequency('1', s,nsam,list);
                        wi = pow(i,exponent);
                        thetaA += wi*i;
                }
        }

        thetaA = thetaA/wSum;
        return thetaA;
}

double achazThetaParabolicWeights(int nsam, int segsites, char **list, int exponent,double center)
{
        int s, frequency( char, int, int, char**);
        double thetaA, i, wi, nd, nnm1, wSum ;

        wSum = 0.0 ;
        for ( i = 1; i < nsam; i++){
            wi = pow((center - i),exponent);
            wSum += wi;
        }

        thetaA = 0.0 ;
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        nnm1 = nd/(nd-1.0);
                        i = frequency('1', s,nsam,list);
                        wi = pow((center - i), (double) exponent);
                        thetaA += wi*i;
                }
        }

        thetaA = thetaA/wSum;
        return thetaA;
}

//pi with a twist of theta H
double achazThetaHPi(int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        double thetaA, i, wi, nd, nnm1, wSum ;

        wSum = 0.0 ;
        for ( i = 1; i < nsam; i++){
            wi = i*(nsam-i);
            wSum += wi;
        }

        thetaA = 0.0 ;
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        nnm1 = nd/(nd-1.0);
                        i = frequency('1', s,nsam,list);
                        wi = 1/(i*i); //i^2 is tw4
                        thetaA += wi*i;
                }
        }

        thetaA = thetaA/wSum;
        return thetaA;
}

//pi minus a sumary stat that is highest when your SFS is U-shaped
double achazTajimasDExtreme(int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        double pi, upsideDownPi, i, wi1, wi2, nd, nnm1, wSum1, wSum2 ;

        wSum1 = 0.0 ;
        wSum2 = 0.0;
        for ( i = 1; i < nsam; i++){
            wi1 = (nsam-i);
            wi2 = pow(((nsam/2.0)-i),2)/i;
            wSum1 += wi1;
            wSum2 += wi2;
        }

        pi = 0.0 ;
        upsideDownPi = 0.0;
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        nnm1 = nd/(nd-1.0);
                        i = frequency('1', s,nsam,list);
                        wi1 = (nsam-i);
                        wi2 = pow(((nsam/2.0)-i),2)/i;
                        pi += wi1*i;
                        upsideDownPi += wi2*i;
                }
        }

        pi = pi/wSum1;
        upsideDownPi = upsideDownPi/wSum2;
        return pi-upsideDownPi;
}

double sigmaAlpha(int n)
{
    int i;
    int sum = 0;
    for (i=1; i<n;i++)
    {
        sum += 1.0/(double)i;
    }
    return sum;
}

double sigmaBeta(int i,int n)
{
    double part1,part2,part3;

    part1 = (double)(2*n)/(double)((n-i+1)*(n-i));
    part2 = sigmaAlpha(n+1) - sigmaAlpha(i);
    part3 = 2.0/(double)(n-i);

    return part1*part2-part3;
}

double sigmaII(int i, int n)
{
    if (2*i < n)
    {
        return sigmaBeta(i+1,n);
    }
    else if (2*i == n)
    {
        return 2.0*(sigmaAlpha(n) - sigmaAlpha(i))/(double)(n-i) - (1.0/(double)(i*i));
    }
    else
    {
        return sigmaBeta(i,n) - (1.0/(double)(i*i));
    }
}

//copied from Achaz
double * compute_HarmonicSums( int n ){
	double *HarmonicSums;
	int i;
	HarmonicSums = (double *)malloc( sizeof(double)*n );
	if(!HarmonicSums)
		fprintf(stderr, "compute_HarmonicSums: malloc error for HarmonicSums\n"), exit(3);
	HarmonicSums[0]=0;
	i=1;
	while(i<n){
		HarmonicSums[i] = HarmonicSums[i-1] + 1.0/i;
		i++;
	}
	return HarmonicSums;
}

//these next 3 functions are copied from Achaz because my versions have a bug somewhere and I am lazy; these seem to work
double beta_FU1995( int i, double *HarmonicSums, int n  ){

	double ai = HarmonicSums[i-1], an = HarmonicSums[n-1];
	double beta=0;
	double nloci= (double)n;

	beta = 2.0 * nloci * ( an + (1.0/nloci) - ai )/(  (nloci-i+1.0 )*(nloci-i) ) - 2.0/(nloci - i);

	return beta;
}

double sigma_ii_FU1995( int i, double *HarmonicSums, int n ){
	double nloci= (double)n;
	double sigma_ii=0;
	double ai = HarmonicSums[i-1], an = HarmonicSums[n-1];
	if( 2*i < n )
	{
		sigma_ii = beta_FU1995( i+1, HarmonicSums, n  );
	}
	else
	{
		if( 2*i == n  )
		{
			sigma_ii = 2.0*(an-ai)/(nloci - i) - 1.0/(i*i);
		}
		else
		{
			sigma_ii = beta_FU1995( i , HarmonicSums, n  ) - 1.0/(i*i);
		}
	}
	return sigma_ii;
}

double sigma_ij_FU1995( int i, int j, double *HarmonicSums, int n ){
	double nloci= (double)n;
	double sigma_ij=0;

	if(i==j){
		return sigma_ii_FU1995( i, HarmonicSums, n );
	}
 	if(i<j){
		int tmp=i;
		i=j;
		j=tmp;
	}
	double  ai=HarmonicSums[i-1],
	       aj=HarmonicSums[j-1],
	       an=HarmonicSums[n-1];
	if( i+j < n )
	{
		sigma_ij = ( beta_FU1995( i+1, HarmonicSums, n  ) -  beta_FU1995( i, HarmonicSums, n  ) ) / 2.0;
	}
	else
	{
		if( i+j == n  )
		{
			sigma_ij  =  ((an - ai)/(nloci - i) + (an - aj)/(nloci - j)) 
			           - ( ( beta_FU1995( i, HarmonicSums, n   ) +  beta_FU1995( j+1 , HarmonicSums, n ) )  / 2.0 ) 
				   - (1.0/(i*j));

		}
		else
		{
			sigma_ij = (( beta_FU1995( j, HarmonicSums, n  ) -  beta_FU1995( j+1, HarmonicSums, n ) )/2.0) - (1.0/(i*j));
		}

	}
	return sigma_ij;

}

double sigmaIJ(int i, int j, int n)
{
    double part1,part2,part3,part4;
    int tmp;

    if (i<j)
    {
        tmp = i;
        i = j;
        j = tmp;
    }
    if (i+j < n)
    {
        return ((sigmaBeta(i+1,n) - sigmaBeta(i,n))/2.0);
    }
    else if (i+j == n)
    {
        part1 = (sigmaAlpha(n) - sigmaAlpha(i)) / (double)(n-i);
        part2 = (sigmaAlpha(n) - sigmaAlpha(j)) / (double)(n-j);
        part3 = (sigmaBeta(i,n) + sigmaBeta(j+1,n)) / 2.0;
        part4 = 1.0/(double)(i*j);
        return part1 + part2 - part3 - part4;
    }
    else
    {
        //printf("%d,%d,sigmaIJ=%f\n",i,j,((sigmaBeta(j,n) - sigmaBeta(j+1,n))/2.0 - 1.0/(double)(i*j)));
        return ((sigmaBeta(j,n) - sigmaBeta(j+1,n))/2.0 - 1.0/(double)(i*j));
    }
}

double achazNeutTestExponentWeights(int nsam, int segsites, char **list, int exponent1, int exponent2, double *harmonicSums)
{
        int i, j, s, frequency( char, int, int, char**);
        double nd, T, Tnum, alphaN, w1sum, w2sum, betaNpart1, betaNpart2, betaN, thetaEstimate, thetaSquaredEstimate, a2;
        int *afs;
        double *w1;
        double *w2;
        double *omega;

        //initializing afs, and weight arrays, and weight sums
        //in each array, element index zero remains uninitialized: not counting monomorphic sites
        afs = (int *) malloc(sizeof(int)*nsam);
        w1 = (double *) malloc(sizeof(double)*nsam);
        w2 = (double *) malloc(sizeof(double)*nsam);
        w1sum = 0.0;
        w2sum = 0.0;
        for (i = 1; i < nsam; i++){
            afs[i] = 0;
            w1[i] = pow(i,exponent1);
            w2[i] = pow(i,exponent2);
            w1sum += w1[i];
            w2sum += w2[i];
        }

        omega = (double *) malloc(sizeof(double)*nsam);
        //initializing omega, again not using the element at index 0
        for (i = 1; i < nsam; i++){
            omega[i] = (w1[i]/w1sum) - (w2[i]/w2sum);
        }

        //compute afs
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        i = frequency('1', s,nsam,list);
                        afs[i] += 1;
                }
        }

        //calculate neutrality test statistic
        Tnum = 0.0 ;
        betaNpart1 = 0.0;
        betaNpart2 = 0.0;
        alphaN = 0.0;
        for ( i = 1; i < nsam; i++){

            Tnum += omega[i]*(double)(i*afs[i]);

            alphaN += i*omega[i]*omega[i];
            betaNpart1 += (double)(i*i)*omega[i]*omega[i]*sigma_ii_FU1995(i,harmonicSums,nsam);
            for (j=i+1; j<nsam; j++)
            {
                betaNpart2 += (double)(2*i*j)*omega[i]*omega[j]*sigma_ij_FU1995(i,j,harmonicSums,nsam);
                //printf("omega[i]=%f,omega[j]=%f,betaNpart2: i=%d,j=%dn=%d,sigmaIJ=%f\n", omega[i],omega[j],i,j,nsam,sigma_ij_FU1995(i,j,harmonicSums,nsam));
            }
        }

        //borrowing some code here from Achaz to get estimates of theta (using thetaW) and theta^2
        a2 = 0;
        for(i=1;  i< nsam;  i++){
		a2 += 1.0/(double)(i*i);
	}
	thetaEstimate =  (double) segsites / harmonicSums[nsam-1];
	thetaSquaredEstimate =  (double) segsites*( (double)segsites-1.0 )/(harmonicSums[nsam-1]*harmonicSums[nsam-1] + a2 );

        betaN=betaNpart1+betaNpart2;
        //printf("betaNpart1=%f,betaNpart2=%f\tbetaN=%f\n",betaNpart1,betaNpart2,betaN);
        //printf("alphaN=%f\n",alphaN);
        //printf("root=%f\n",pow(alphaN*thetaEstimate+(betaN*thetaSquaredEstimate),0.5));
        T = Tnum / sqrt(alphaN*thetaEstimate+(betaN*thetaEstimate*thetaSquaredEstimate));

        free(afs);
        free(w1);
        free(w2);
        free(omega);
        return T;
}


/* Fay's theta_H  */
//again corrected for sampleSize variation
double	thetah( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd  ;

	pi = 0.0 ;

	nd = nsam;
	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list) ;
		nd = sampleSizeSite(s,nsam,list);
		if(nd > 1){
			pi += (p1*p1)/( nd*(nd-1.0) ) ; 
		}
	}
	return(pi*2.0) ;
}


int frequency( char allele,int site,int nsam,  char **list){
	int i, count=0;
	for( i=0; i<nsam; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
	return( count);
}        

//sampleSizeSite -- returns the sampleSize at a site corrected for missing data
int sampleSizeSite(int site, int nsam, char **list){
	return(nsam - frequency('N',site,nsam,list));	
}


//gets the right number of segSites when there are Ns
int segSites(int segsites, int nsam, char **list){
	int i, ss = 0;

	for(i=0; i < segsites; i++){
		ss += (((frequency('1', i, nsam, list) > 0) && (frequency('0', i, nsam, list) >0 )) ? 1:0);
	}
	return(ss);
}

//gets the derived site frequency spectrum; fixations treated as monomorphic
//derived_counts must be an int array of length nsam
//derived_counts[i] is the fraction of sites with derived allele present in i chromosomes
//number of monomorphic sites is sorted in derived_counts[0]
void getSiteFreqSpec(int segsites, int nsam, char**list, int nSites, int *derived_counts)
{
        int i;
        int freq;

        for (i=0; i<nsam; i++)
        {
            derived_counts[i] = 0;
        }

        int polycount = 0;
        for (i=0; i<segsites; i++)
        {
                freq = frequency('1', i, nsam, list);
                if (freq > 0 && freq < nsam)
                {
                        polycount++;
                        derived_counts[freq] += 1;
                }
        }
        derived_counts[0] = nSites-polycount;
}

//gets the derived site frequency spectrum; fixations treated as monomorphic
//derived_counts must be an int array of length nsam
//derived_counts[i] is the fraction of sites with derived allele present in i chromosomes
//number of monomorphic sites is sorted in derived_counts[0]
void getSiteFreqSpecWindow(int segsites, int nsam, char**list, int nSites, int *derived_counts, double *pos, double low, double high)
{
        int i;
        int freq;

        for (i=0; i<nsam; i++)
        {
            derived_counts[i] = 0;
        }

        int polycount = 0;
        for (i=0; i<segsites; i++)
        {
		if (pos[i] >low && pos[i] <= high)
		{
                	freq = frequency('1', i, nsam, list);
                	if (freq > 0 && freq < nsam)
                	{
                        	polycount++;
                        	derived_counts[freq] += 1;
                	}
		}
        }
        derived_counts[0] = nSites-polycount;
}

//counts the number of haplotypes, and gets their frequencies (stored in haplotype_counts)
//haplotype_counts must be an int array of length nsam
//haplotype_counts[i] is the number of haplotypes found in exactly i+1 individuals
int getHaplotypeFreqSpec(int segsites, int nsam, char **list, int *haplotype_counts)
{
        int i;
        int j;
        int k;
        int haplotype_found;
        int allsame;
        int freq;

        int n_haplotypes = 0;
        char haplotypes[nsam][segsites+1];
        int haplotype_occurrences[nsam];

        for(i=0; i<nsam; i++)
        {
            haplotype_counts[i] = 0;
            haplotype_occurrences[i] = 0;
        }

        for(i=0; i<nsam; i++)
        {
                haplotype_found = 0;
                for(j=0; j<n_haplotypes; j++)
                {
                        allsame = 1;
                        for(k=0; k<segsites; k++)
                        {
                                if(haplotypes[j][k] != list[i][k])
                                {
                                        if(haplotypes[j][k] == 'N' )
                                                haplotypes[j][k] = list[i][k];

                                        else if(list[i][k] != 'N')
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
                        for(j=0; j<segsites; j++)
                        {
                                haplotypes[n_haplotypes-1][j] = list[i][j];
                        }
                        haplotypes[n_haplotypes-1][segsites]='\0';
                        haplotype_occurrences[n_haplotypes-1]=1;
                }
        }

        for (i=0; i<n_haplotypes; i++)
        {
                freq = haplotype_occurrences[i];
                if (freq > 0 && freq <= nsam)
                {
                        haplotype_counts[freq-1] += 1;
                }
        }

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

//gets H12, H1, and H2 in windows
void petrovHStatsWindow(int segsites, int nwins, double *posit, double *winsH12, double *winsH1, double *winsH2, int nsam, char **list)
{
	int i;
	int j;
	int k;
	int haplotype_found, freq;
	int allsame;
	float start, delta, end;
	int n_haplotypes = 0;
	char haplotypes[nsam][segsites+1];
        int haplotype_occurrences[nsam];
        int haplotype_counts[nsam];
        int wcount = 0;
	start = 0;
        delta = 1.0 / nwins;
        end = delta;

        while(start < 1.0)
        {
                for(i=0; i<nsam; i++)
                {
                    haplotype_counts[i] = 0;
                    haplotype_occurrences[i] = 0;
                }
                for(i=0; i<nsam; i++)
                {
                    haplotype_found = 0;
                    for(j=0; j<n_haplotypes; j++)
                    {
                        allsame = 1;
                        for(k=0; k<segsites; k++)
                        {
                            if (posit[k] > start && posit[k] <= end)
                            {
                                if(haplotypes[j][k] != list[i][k])
                                {
                                    if(haplotypes[j][k] == 'N' )
                                        haplotypes[j][k] = list[i][k];
                                    else if(list[i][k] != 'N')
                                    {
                                        allsame = 0;
                                        break;
                                    }
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
                        for(j=0; j<segsites; j++)
                        {
                            haplotypes[n_haplotypes-1][j] = list[i][j];
                        }
                        haplotypes[n_haplotypes-1][segsites]='\0';
                        haplotype_occurrences[n_haplotypes-1]=1;
                    }
                }
                for (i=0; i<n_haplotypes; i++)
                {
                    freq = haplotype_occurrences[i];
                    if (freq > 0 && freq <= nsam)
                    {
                        haplotype_counts[freq-1] += 1;
                    }
                }
                winsH12[wcount]=petrovH12(haplotype_counts,nsam);
                winsH1[wcount]=petrovH1(haplotype_counts,nsam);
                winsH2[wcount]=petrovH2(haplotype_counts,nsam);
                wcount++;
                start += delta;
                end += delta;
                n_haplotypes = 0;
        }
}

double meanEHH(int segsites, double *posit, double delta, int nsam, char **list)
{
	int i, j, k, l;
	int allsame;
	float start, end;
        int sameCount;
        int totalCount;
        float ehhSum = 0;
        float ehhCount = 0;
        for(i=0;i<segsites;i++)
        {
            sameCount = 0;
            totalCount = 0;
            start = posit[i]-delta;
            end = posit[i]+delta;
            if (start >= 0 && end <= 1)
            {
                for(j=0; j<nsam; j++)
                {
                    if (list[j][i] == '1')
                    {
                        for(k=0;k<j;k++)
                        {
                            if (list[k][i] == '1')
                            {
                                allsame = 1;
                                //see if j and k are the same at all segsites in the range
                                for(l=0; l<segsites; l++)
                                {
                                    if (posit[l] >= start && posit[l] <= end && l != i)
                                    {
                                        if (list[j][l] != list[k][l])
                                        {
                                            allsame = 0;
                                            break;
                                        }
                                    }
                                }
                                if (allsame)
                                {
                                    sameCount++;
                                }
                                totalCount++;
                            }
                        }
                    }
                }
                if (totalCount > 0)
                {
                    ehhSum += sameCount/(float)totalCount;
                }
                ehhCount += 1;
            }
        }
        return ehhSum/ehhCount;
}

double meanREHH(int segsites, double *posit, double delta, int nsam, char **list)
{
	int i, j, k, l;
	int allsame;
	float start, end;
        int sameCountDer, sameCountAnc;
        int totalCountDer, totalCountAnc;
        float ehhSum = 0;
        float ehhCount = 0;
        for(i=0;i<segsites;i++)
        {
            sameCountDer = 0;
            totalCountDer = 0;
            sameCountAnc = 0;
            totalCountAnc = 0;
            start = posit[i]-delta;
            end = posit[i]+delta;
            if (start >= 0 && end <= 1)
            {
                for(j=0; j<nsam; j++)
                {
                    if (list[j][i] == '1')
                    {
                        for(k=0;k<j;k++)
                        {
                            if (list[k][i] == '1')
                            {
                                allsame = 1;
                                //see if j and k are the same at all segsites in the range
                                for(l=0; l<segsites; l++)
                                {
                                    if (posit[l] >= start && posit[l] <= end && l != i)
                                    {
                                        if (list[j][l] != list[k][l])
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
                    else if (list[j][i] == '0')
                    {
                        for(k=0;k<j;k++)
                        {
                            if (list[k][i] == '0')
                            {
                                allsame = 1;
                                //see if j and k are the same at all segsites in the range
                                for(l=0; l<segsites; l++)
                                {
                                    if (posit[l] >= start && posit[l] <= end && l != i)
                                    {
                                        if (list[j][l] != list[k][l])
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
                if (totalCountAnc > 0 && totalCountDer > 0 && sameCountAnc > 0)
                {
                    ehhSum += (sameCountDer/(float)totalCountDer) / (sameCountAnc/(float)totalCountAnc);
                }
                ehhCount += 1;
            }
        }
        return ehhSum/ehhCount;
}

//counts number of haplotypes
int nHaplotypes(int segsites, int nsam, char **list)
{
	int i;
	int j;
	int k;
	int haplotype_found;
	int allsame;
	
	int n_haplotypes = 0;
	char haplotypes[nsam][segsites+1];
	
	for(i=0; i<nsam; i++)
	{
		haplotype_found = 0;
		for(j=0; j<n_haplotypes; j++)
		{
			allsame = 1;
			for(k=0; k<segsites; k++)
			{
				if(haplotypes[j][k] != list[i][k])
				{
					if(haplotypes[j][k] == 'N' )
						haplotypes[j][k] = list[i][k];
												
					else if(list[i][k] != 'N')
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
			for(j=0; j<segsites; j++)
					haplotypes[n_haplotypes-1][j] = list[i][j];
                        haplotypes[n_haplotypes-1][segsites]='\0';
		}
	}
		
	return n_haplotypes;
}

double countAlleleDiffsForSnpPair(int nsam, int i, int j, char **list)
{
	int k;
	int numDiffs = 0, numSames=0;
	for (k=0; k<nsam; k++)
	{
		if (list[k][i] != 'N' && list[k][j] != 'N')
		{
			if (list[k][i] != list[k][j])
			{
				numDiffs++;
			}
			else
			{
				numSames++;
			}
		}
	}
	if (numDiffs<numSames)
		return numDiffs;
	else
		return numSames;
}

double sStarSnpDist(int nsam, int i, int j, int physLen, double *posit, char **list)
{
	int numDiffs;
	numDiffs = countAlleleDiffsForSnpPair(nsam, i, j, list);
	//printf("%d, %d, %f; numDiffs: %d\n", i, j, posit[i]-posit[j], numDiffs);
	if (numDiffs == 0)
	{
		return (posit[i]-posit[j]) * physLen + 5000;
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

double sStar(int nsam, int segsites, int physLen, double *posit, char **list)
{
	int i,j;
	double globalMaxVal, maxVal, currVal;
	double sStarArray[segsites];
	sStarArray[0] = 0.0;
	globalMaxVal = 0.0;
	for (i=1; i<segsites; i++)
	{
		maxVal = 0.0;
		for (j=0; j<i; j++)
		{
			currVal = sStarArray[j] + sStarSnpDist(nsam, i, j, physLen, posit, list);
			//printf("%d, %d; %f, currVal: %f\n", i, j, sStarArray[j], currVal);
			if (currVal > maxVal)
			{
				maxVal = currVal;
			}
		}
		sStarArray[i] = maxVal;
		if (maxVal > globalMaxVal)
		{
			globalMaxVal = maxVal;
		}
	}
	return globalMaxVal;
}

double dij(int i, int j, int nsam, char** list){
	
	double pi  = 0.0;
	double pj  = 0.0;
	double pij = 0.0;
	double count = 0.0;
	int k;
	for(k=0; k<nsam; k++){
		if(list[k][i] != 'N' && list[k][j] != 'N'){
			if(list[k][i] == '1')
				pi++;
			
			if(list[k][j] == '1')
				pj++;
			
			if(list[k][i] == '1' && list[k][j] == '1')
				pij++;	
			count += 1;
		}
	}
	if (count == 0){
		return(0);
	}
	else{
		pi  /= count;
		pj  /= count;
		pij /= count;
	
		double Dij = pij - (pi*pj);
	
		return (Dij*Dij) / ((pi*(1.0-pi)) * (pj*(1.0-pj)));	
	}
}


/* implements the ZnS statistic of J Kelly (1997), "A Test of Neutrality Based */
/* on Interlocus Associations." Genetics, 146: 1197-1206.                      */
double ZnS(int segsites, int nsam, char** list){

	if(segsites < 2)
		return 0.0;

	int i,j, s;
	double sum = 0.;
	
	s = segSites(segsites, nsam, list);
	for(i=0; i<segsites-1; i++){
                if (frequency('1', i, nsam, list) < nsam){
			for(j=i+1; j<segsites; j++){
        	        	if (frequency('1', j, nsam, list) < nsam){
					sum += dij(i, j, nsam, list);
				}
			}
		}
	}
	return (2.0 / (double)(s * (s-1))) * sum;
}

/*omega statistic from Kim and Nielsen (2003)
** not robust to missing data
*/
double omegaWithTable(int left, int segsites,int nsam, char** list, double** dijTable){
	int i,j, s;
	double sum,sumL,sumR,comp,denom;
	double numer;
	
	sum = sumL = sumR =comp=denom=0;
	if(segsites < 3)
		return 0.0;
	s = segSites(segsites, nsam, list);
	//calculate: 
	// sum for denom-- all pairwise r2
	// sumL and sumR for numerator -- blockwise r2
	for(i=0; i<segsites-1; i++){
		for(j=i+1; j<segsites; j++){
			comp = dijTable[i][j];
			if(i < left && j >= left)sum += comp;
			if(i < left && j < left) sumL += comp;
			if(i >= left && j >= left) sumR += comp;	
		}
	}
	denom = sum * (1.0/(left*(s-left)));
	numer = 1.0 / ((left*(left-1)/2) + ((s-left)*(s-left-1)/2));
	numer *= (sumL+sumR);
	//printf("n/d: %f d: %f n %f sumL: %f sumR: %f term: %f left: %d\n",numer/denom,denom,numer,sumL,sumR,1.0 / ((left*(left-1)/2) + ((s-left)*(s-left-1)/2)),left);
	return(numer/denom);
}

/*omega statistic from Kim and Nielsen (2003)
** not robust to missing data
*/
double omega(int left, int segsites,int nsam, char** list){
        int i,j, s;
        double sum,sumL,sumR,comp,denom;
        double numer;

        sum = sumL = sumR =comp=denom=0;
        if(segsites < 3)
                return 0.0;
        s = segSites(segsites, nsam, list);
        //calculate: 
        // sum for denom-- all pairwise r2
        // sumL and sumR for numerator -- blockwise r2
        for(i=0; i<segsites-1; i++){
                for(j=i+1; j<segsites; j++){
                        comp = dij(i, j, nsam, list);
                        if(i < left && j >= left)sum += comp;
                        if(i < left && j < left) sumL += comp;
                        if(i >= left && j >= left) sumR += comp;
                }
        }
        denom = sum * (1.0/(left*(s-left)));
        numer = 1.0 / ((left*(left-1)/2) + ((s-left)*(s-left-1)/2));
        numer *= (sumL+sumR);
        return(numer/denom);
}

/*Get omega at known fixation position. I.E. center of simulation*/
double omegaCenter(int siteIdx , int segsites,int nsam, char** list){
	int i,j, s;
	double sum,sumL,sumR,comp,denom;
	double numer;
	
	sum = sumL = sumR =comp=denom=0;
	if(segsites < 3)
		return 0.0;
	s = segSites(segsites, nsam, list);
	//calculate: 
	// sum for denom-- all pairwise r2
	// sumL and sumR for numerator -- blockwise r2
	for(i=0; i<segsites-1; i++){
		for(j=i+1; j<segsites; j++){
			comp = dij(i, j, nsam, list);
			if(i < siteIdx && j >= siteIdx)sum += comp;
			if(i < siteIdx && j < siteIdx) sumL += comp;
			if(i >= siteIdx && j >= siteIdx) sumR += comp;	
		}
	}
	denom = sum * (1.0/(siteIdx*(s-siteIdx)));
	numer = 1.0 / ((siteIdx*(siteIdx-1)/2) + ((s-siteIdx)*(s-siteIdx-1)/2));
	numer *= (sumL+sumR);
//	printf("d: %f n %f sumL: %f sumR: %f term: %f left: %d\n",denom,numer,sumL,sumR,1.0 / ((left*(left-1)/2) + ((s-left)*(s-left-1)/2)),left);
	if (isnan(denom)){
        return 0.0;
    }
    else{
        return(numer/denom);    
    }
}


/*omegaMax -- goes through all possible site divisions to maximize omega
// Kim and Nielsen (2003)
// This version builds a table of all pairwise r^2 vals for faster downstream computation
*/
double omegaMax(int segsites,int nsam, char** list){
	int l, i, j;
	double max= 0;
	double tmp=0;
	double **dijTable;

	dijTable = (double **)malloc( sizeof(double)*segsites );

	if(segsites < 3)
		return(0);
	for(i=0; i<segsites-1; i++){
		dijTable[i] = (double *) malloc( sizeof(double)*segsites );
		for(j=i+1; j<segsites; j++){
			dijTable[i][j] = dij(i, j, nsam, list);
                        //printf("dijTable[i][j]: %f\n", dijTable[i][j]);
		}
	}
	for(l=3;l<segsites-2;l++){
		tmp = omegaWithTable(l, segsites, nsam,list, dijTable);
		if(tmp > max){
			max = tmp;
		}
	}

	for(i=0; i < segsites-1; i++){
                free(dijTable[i]);
        }
        free(dijTable);

	return(max);
}

/////////////////////////
//Two Site Utils
//////
////
//
//sampleConfig-- fills vector with the sample configuration
// for sites i and j
void sampleConfig(int i, int j, int nsam, char** list, int *config){
	int p1, p2, x11, k;
	p1 = p2 = x11 = 0;

	for(k=0; k<nsam; k++){
		if(list[k][i] != 'N' && list[k][j] != 'N'){
			if(list[k][i] == '1')
				p1++;
			
			if(list[k][j] == '1')
				p2++;
			
			if(list[k][i] == '1' && list[k][j] == '1')
				x11++;	
		}
	}
	if(p1 > p2){
		config[0] = p1;
		config[1] = p2;
	}
	else{
		config[0] = p2;
		config[1] = p1;
	}
	config[2] = x11;
}

void printPairwiseSampleConfigs(int segsites, int nsam, char **list, double *posit, int nsites){
	int i,j, config[3];

	if(segsites < 2){
		return;
	}
	else{
		for(i=0; i<segsites-1; i++){
			for(j=i+1; j<segsites; j++){
				sampleConfig(i, j, nsam, list,config);
				printf("%d %d %d %d %d\n",nsam,config[0],config[1],config[2],(int)floor((posit[j]-posit[i]) * nsites));
			}	
		}
	}
}

//sampleConfig2Popn- fills vector with the sample configuration
// for sites i and j; 6 dimensions for 2 population 2 site sample config
void sampleConfig2Popn(int i, int j, int nsam, int popnSize1, char** list, int *config){
	int p1, p2, x11,p3,p4,y11, k;
	p1 = p2 = x11 = p3 = p4 = y11 = 0;

	for(k=0; k<nsam; k++){
		if(list[k][i] != 'N' && list[k][j] != 'N'){
			if(k < popnSize1){
				if(list[k][i] == '1')
					p1++;
			
				if(list[k][j] == '1')
					p2++;
			
				if(list[k][i] == '1' && list[k][j] == '1')
					x11++;	
			}
			else{
				if(list[k][i] == '1')
					p3++;
			
				if(list[k][j] == '1')
					p4++;
			
				if(list[k][i] == '1' && list[k][j] == '1')
					y11++;
			}
		}
	}
//	if(p1 > p2){
		config[0] = p1;
		config[1] = p2;
//	}
//	else{
//		config[0] = p2;
//		config[1] = p1;
//	}
	config[2] = x11;
//	if(p3 > p4){
		config[3] = p3;
		config[4] = p4;
//	}
//	else{
//		config[3] = p4;
//		config[4] = p3;
//	}
	config[5] = y11;
}

void printPairwiseSampleConfigs2Popn(int segsites, int nsam, int popnSize1, char **list, double *posit, int nsites){
	int i,j, config[6];

	if(segsites < 2){
		return;
	}
	else{
		for(i=0; i<segsites-1; i++){
			for(j=i+1; j<segsites; j++){
				sampleConfig2Popn(i, j, nsam, popnSize1,list,config);
				printf("%d %d %d %d %d %d %d %d %d\n",popnSize1,nsam-popnSize1,config[0],config[1],config[2],config[3],config[4],config[5],(int)floor((posit[j]-posit[i]) * nsites));
			}	
		}
	}
}
/****************
///   Sub popn versions
******************/

//frequencySub-- allows for arbitrary allele indexes as you would need for sub pops
int frequencySub(char allele, int site, int startAllele, int stopAllele, char **list){
	int i, count=0;
	for( i=startAllele; i<stopAllele; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
	return( count);
}


//gets the right number of segSites when there are Ns
int segSitesSub(int segsites, int nsam, int startAllele, int stopAllele, char **list){
	int i, ss = 0;

	for(i=0; i < segsites; i++){
		ss += (((frequencySub('1', i, startAllele, stopAllele, list) > 0) && (frequencySub('0', i, startAllele, stopAllele, list) >0 )) ? 1:0);
	}
	return(ss);
}

//gets the right number of segSites when there are Ns
void privateSegSitesInTwoPopns(int segsites, int nsam, int stopAllele, int *private1, int *private2, char **list){
	int i, isSeg1, isSeg2;
	(*private1) = 0;
	(*private2) = 0;

	for(i=0; i < segsites; i++){
		isSeg1 = 0;
		isSeg2 = 0;
		if ((frequencySub('1', i, 0, stopAllele, list) > 0) && (frequencySub('0', i, 0, stopAllele, list) > 0))
			isSeg1 = 1;
		if ((frequencySub('1', i, stopAllele, nsam, list) > 0) && (frequencySub('0', i, stopAllele, nsam, list) > 0))
			isSeg2 = 1;
		if (isSeg1 && !isSeg2)
			(*private1)++;
		else if (isSeg2 && !isSeg1)
			(*private2)++;
	}
}


//sampleSizeSiteSub -- returns the sampleSize at a site corrected for missing data
//in startAllele to stopAllele rows
int sampleSizeSiteSub(int site, int nsam, int startAllele, int stopAllele, char **list){
	return(stopAllele - startAllele  - frequencySub('N',site,startAllele, stopAllele,list));	
}

double *hetVec1Popn(int segsites, int n1, int physLen, int *vecLen, double *posit, char **list)
{
	(*vecLen) = 0;
	int i, j, k;
	double diffs;
	double *vec;

	vec = (double *) malloc(sizeof(double) *n1*(n1-1)/2);
	for(i=0; i<n1-1;i++){
		for(j=i+1;j<n1;j++){
			diffs = 0.;
			for(k=0;k<segsites;k++){
				if(list[i][k] != list[j][k] && list[i][k] != 'N' && list[j][k] != 'N'){
					diffs += 1;
				}
			}
			vec[(*vecLen)++] = diffs/(double)physLen;
		}
	}
	return vec;
}

void reorderListIntoClusters(int *membership1, int g1Size, int *membership2, int g2Size, char **list){
	int i, j = 0;
	char *tmp;
	for (i=0; i<g1Size; i++)
	{
		//printf("checking to see if we have to swap member %d which is %d\n", i, membership1[i]);
		if (membership1[i] >= g1Size)
		{
			//printf("indeed we must (%d >= %d)! Looking for a suitable parter starting with %d which is %d\n", membership1[i], g1Size, j, membership2[j]);
			while (membership2[j] >= g1Size)
			{
				j++;
				//printf("nope! moving on to %d which is %d\n", j, membership2[j]);
			}
			if (j >= g2Size) fprintf(stderr, "reorderListIntoClusters has a bug!\n"), exit(-1);
			//printf("okay we have a winner! (%d, which is %d)\n", j, membership2[j]);
			tmp = list[membership1[i]];
			list[membership1[i]] = list[membership2[j]];
			list[membership2[j]] = tmp;
			j++;
		}
	}
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
	int i, j, minI, maxPairI = -1, maxPairJ = -1, pairIndex = 0, m1Index = 0, m2Index = 0;
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
	*(g1Size) = m1Index;
	*(g2Size) = m2Index;
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

void clusterSeqsFromUnsortedHetVec(double *hetVec, int n1, int *g1Size, int *g2Size, char **list){
	int *membership1 = (int *) malloc (n1*sizeof(int));
	int *membership2 = (int *) malloc (n1*sizeof(int));
	int i, swap, tmp, *tmpArray;
	swap = assignClusters(hetVec, n1, g1Size, membership1, g2Size, membership2);
	if (swap)
	{
		//printf("swapping!!\n");
		tmp = (*g1Size);
		(*g1Size) = (*g2Size);
		(*g2Size) = tmp;
		tmpArray = membership1;
		membership1 = membership2;
		membership2 = tmpArray;
	}
	//printf("before:\n");
	//for (i=0; i<n1; i++) printf("%s\n", list[i]);
	//for(i=0;i<*(g1Size);i++) printf("g1: %d\n", membership1[i]);
	//for(i=0;i<*(g2Size);i++) printf("g2: %d\n", membership2[i]);
	reorderListIntoClusters(membership1, *g1Size, membership2, *g2Size, list);
	//printf("after:\n");
	//for (i=0; i<n1; i++) printf("%s\n", list[i]);
	free(membership1);
	free(membership2);
}

double nucdivSub( int nsam, int segsites, int startAllele, int stopAllele, char **list){
	int s;
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;


	for( s = 0; s <segsites; s++){
		nd = sampleSizeSiteSub(s,nsam,startAllele,stopAllele,list);
		if (nd > 1){
			nnm1 = nd/(nd-1.0) ;
			p1 = frequencySub('1', s,startAllele,stopAllele,list)/nd ;
			pi += 2.0*p1*(1.0 -p1)*nnm1 ;
		}
	}
	return( pi ) ;
}


//fills a vector size l with values of nucdiv in "windows". for subpops
void nucdivSubWindow( int nwins, double *posit, double *output, int nsam, int segsites,int startAllele, int stopAllele, char **list)
{
	int s, frequency( char, int, int, char**);
	int wcount = 0;
	double pi, p1, nd, nnm1  ;
	double start, end, delta;
	start = 0;
	delta = 1.0 / nwins;
	end = delta;
	
	while(start < 1.0){
		pi = 0.0 ;
		for( s = 0; s <segsites; s++){
			if(posit[s] <=end && posit[s] > start){
				nd = sampleSizeSiteSub(s,nsam,startAllele,stopAllele,list);
				if (nd > 1){
					nnm1 = nd/(nd-1.0) ;
					p1 = frequencySub('1', s,startAllele,stopAllele,list)/nd ;
					pi += 2.0*p1*(1.0 -p1)*nnm1 ;
				}
			}
		}
		output[wcount++]=pi;
		start += delta;
		end += delta;
	}
}


void fst2SubsWindow(int nwins, double *posit, double *output,int segsites, int nsam, int start1, int stop1, int start2, int stop2, char **list){
	double h1[nwins], h2[nwins], hTot[nwins], hW;

	int i;
	
	nucdivSubWindow(nwins,posit,h1,nsam, segsites, start1, stop1, list);
	nucdivSubWindow(nwins,posit,h2,nsam, segsites,start2,stop2,list);
	nucdivWindow(nwins,posit,hTot,nsam, segsites,list);
	
	for(i=0;i<nwins;i++){
		hW = (((stop1 - start1) * h1[i]) + ((stop2-start2) * h2[i])) / (double) ((stop1-start1) +(stop2-start2));

		if(h1[i] == 0.0 && h2[i] == 0.0)output[i]=0.0;
		else output[i] = (hTot[i]-hW) / hTot[i];
	}
	
}

/* Fay's theta_H  */
//again corrected for sampleSize variation
double thetahSub( int nsam, int segsites, int startAllele, int stopAllele, char **list){
	int s;
	double pi, p1, nd  ;

	pi = 0.0 ;

	nd = nsam;
	for( s = 0; s <segsites; s++){
		p1 = frequencySub('1', s,startAllele,stopAllele,list) ;
		nd = sampleSizeSiteSub(s,nsam,startAllele,stopAllele,list);
		if(nd > 1 && p1 != (stopAllele-startAllele)){
			pi += (p1*p1)/( nd*(nd-1.0) ) ; 
		}
	}
	return(pi*2.0) ;
}

//counts number of haplotypes
int nHaplotypesSub(int segsites, int nsam, int startAllele, int stopAllele, char **list)
{
	int i;
	int j;
	int k;
	int haplotype_found;
	int allsame;
	
	int n_haplotypes = 0;
	char haplotypes[ stopAllele - startAllele][segsites];

	for(i=startAllele; i<stopAllele; i++)
	{
		haplotype_found = 0;
		for(j=0; j<n_haplotypes; j++)
		{
			allsame = 1;
			for(k=0; k<segsites; k++)
			{
				if(haplotypes[j][k] != list[i][k])
				{
					if(haplotypes[j][k] == 'N' )
						haplotypes[j][k] = list[i][k];
												
					else if(list[i][k] != 'N')
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
			for(j=0; j<segsites; j++)
					haplotypes[n_haplotypes-1][j] = list[i][j];
		}
	}
		
	return n_haplotypes;
}

double dijSub(int i, int j, int nsam, int startAllele, int stopAllele, char** list){
	
	double pi  = 0.0;
	double pj  = 0.0;
	double pij = 0.0;
	double count = 0.0;
	double Dij = 0.0;
	double denom =0.0;
	int k;
	for(k=startAllele; k<stopAllele; k++){
		if(list[k][i] != 'N' && list[k][j] != 'N'){
			if(list[k][i] == '1')
				pi++;
			
			if(list[k][j] == '1')
				pj++;
			
			if(list[k][i] == '1' && list[k][j] == '1')
				pij++;	
			count += 1;
		}
	}
	if (count == 0){
		return(0);
	}
	else{
		pi  /= count;
		pj  /= count;
		pij /= count;
	
		Dij = pij - (pi*pj);
	    denom = (pi*(1.0-pi)) * (pj*(1.0-pj));
		if(denom == 0)
			return(0);
		return ((Dij*Dij) / denom);	
	}
}

/* implements the ZnS statistic of J Kelly (1997), "A Test of Neutrality Based */
/* on Interlocus Associations." Genetics, 146: 1197-1206.                      */
double ZnSSub(int segsites, int nsam, int startAllele, int stopAllele, char** list){

	if(segSitesSub(segsites,nsam,startAllele,stopAllele,list) < 2){
		return(0);
	}

	int i,j, s;
	double sum = 0.;
	
	s = segSitesSub(segsites, nsam, startAllele,stopAllele,list);
	for(i=0; i<segsites-1; i++){
		for(j=i+1; j<segsites; j++){
			sum += dijSub(i, j, nsam,startAllele,stopAllele, list);
		}
	}
	
	return (2.0 / (double)(s * (s-1))) * sum;
}

//fst--
double fst2Subs(int segsites, int nsam, int start1, int stop1, int start2, int stop2, char **list){
	double h1, h2, hTot, hW;
	double f;
	
	h1 = nucdivSub(nsam, segsites, start1, stop1, list);
	h2 = nucdivSub(nsam, segsites,start2,stop2,list);
	hTot = nucdiv(nsam, segsites,list);
	hW = (((stop1 - start1) * h1) + ((stop2-start2) * h2)) / (double) ((stop1-start1) +(stop2-start2));
	//printf("%f %f %f %f\n",h1,h2,hTot,hW);
	if(h1 == 0.0 && h2 == 0.0)
		return(0.0);
	f = (hTot-hW) / hTot;
//	if(f < 0)
//		return(0.0); going to allow this to return negative values... might be better for estimation?
	return(f);
	
}
//Dxy statistic from Takahata and Nei 1985
double Dxy(int segsites,int nsam, int n1, int n2, char **list){
	int i,j,tmp;
	double sum = 0.0;
	double tmpDist;
	int denom = 0;
	
	tmp=n1+n2;
	for(i=0;i<n1;i++){
		for(j=n1;j<tmp;j++){
			tmpDist = seqDist_Snn(segsites,i,j,list);
			if (tmpDist >= 0.0){
				sum += tmpDist;
				denom++;
			}
		}
	}
	return(sum / (float) denom);
}

double Dxy_mean(int segsites,int nsam, int n1, int n2, char **list){
	int i,j,tmp;
	double sum = 0.0;
	double ncomps = 0.0;
	double tmpDist;
	tmp=n1+n2;

	for(i=0;i<n1;i++){
		for(j=n1;j<tmp;j++){
			tmpDist = seqDist_Snn(segsites,i,j,list);
			if (tmpDist >= 0.0){
				sum += tmpDist;
				ncomps++;
			}
		}
	}
	return(sum /ncomps);
}

//Dxy_min statistic used in Garrigan's Gmin
double Dxy_min(int segsites,int nsam, int n1, int n2, char **list){
	int i,j,tmp;
	double min = 666666666.0;
	double tmpVal;
	
	tmp=n1+n2;
	for(i=0;i<n1;i++){
		for(j=n1;j<tmp;j++){
			tmpVal = seqDist_Snn(segsites,i,j,list);
			if(tmpVal < min && tmpVal >= 0.0) min = tmpVal;
		}
	}
	return(min);
}

//pairwiseDistRankAmongSampleRange takes a number of pairwise differences and a range of rows in the alignment
//and returns the number of pairs of samples within this range that have a lesser or equal number of differences
double pairwiseDistRankAmongSampleRange(int segsites, int diffs, int firstSample, int numSamples, double *hetVar, char **list){
	int i, j, currDiffs, leCount, numComps, sumDiffs;
	double numer, meanDiffs;
	double *hetLs;
	sumDiffs = leCount = numComps = 0;
	(*hetVar)=0.0;
	
	hetLs = (double *)malloc(sizeof(double)*(numSamples*(numSamples-1)/2));
	for(i=firstSample; i<firstSample+numSamples-1 ;i++){
		for(j=i+1; j<firstSample+numSamples ;j++){
			currDiffs = seqDist_Snn(segsites, i, j, list);
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
	return leCount/(float)numComps;
}

//Dxy_vector-- fills a vector of dxy values sorted
//TODO: filter out large values
void Dxy_vector(int segsites,int nsam, int n1, int n2, double *v,char **list){
	int i,j,tmp;
	int count = 0;
	
	tmp=n1+n2;
	for(i=0;i<n1;i++){
		for(j=n1;j<tmp;j++){
			v[count] = seqDist_Snn(segsites,i,j,list);
			if (v[count] < 0.0){
				v[count] = NAN;
			}
			count++;
		}
	}
	qsort(v,count, sizeof(v[0]),cmp_doubles);
}

//cmp_doubles -- for use with qsort
int cmp_doubles(const void *x, const void *y){
	double xx = *(double*)x;
	double yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}


//pairwiseDistances-- fills an array of doubles with all pairwise dists per site
void pairwiseDistances(int segsites,int nsam, double *dists, char **list){
	int i,j, count;
	
	count = 0;
	for(i=0;i<nsam-1;i++){
		for(j=i+1;j<nsam;j++){
			if(segsites == 0) dists[count++] = 0.0;
			else dists[count++] = (double) seqDist_Snn(segsites,i,j,list) / segsites;
			if (dists[count-1] < 0) dists[count-1] = NAN;
		}
	}
}


//******************
//Snn -- Hudson's Snn statistic from Hudson (2000)
double Snn(int segsites,int nsam, int n1, int n2, char **list){
	double count = 0;
	int i;
	
	for(i=0;i<nsam; i++){
		count += xij_Snn(segsites,nsam,i,n1,n2,list);
	}
	count /= (double) nsam;
	return(count);
}

//this counts the actually proportions of nearest neighbors in same pop
double xij_Snn(int segsites,int nsam, int seqIndex1, int n1, int n2, char **list){
	
	int i;
	double minWith = 666666666.0;
	double minBet  = 666666666.0;
	double tmp;
	int minCountW = 0;
	int minCountB = 0;
	
	if(seqIndex1 < n1){
		//get within min and count
		for(i = 0; i < n1; i++){
			if(i != seqIndex1){
				tmp = seqDist_Snn(segsites,seqIndex1,i,list);
				if (tmp >= 0.0){
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
		for(i = n1; i < nsam; i++){
			tmp = seqDist_Snn(segsites,seqIndex1,i,list);
			if (tmp >= 0.0){
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
		for(i = n1; i < nsam; i++){
			if(i != seqIndex1){
				tmp = seqDist_Snn(segsites,seqIndex1,i,list);
				if (tmp >= 0.0){
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
			tmp = seqDist_Snn(segsites,seqIndex1,i,list);
			if (tmp >= 0.0){
				if( tmp< minBet){
					minCountB = 1;
					minBet = tmp;
				}
				if( tmp == minBet)
					minCountB += 1;
			}
		}
	}

	if(minWith < minBet){
		return(1.0);
	}
	if(minWith == minBet){
		return(minCountW / (double) (minCountW+minCountB));
	}
	return(0);

}

//seqDist_Snn is the metric for Snn; note that a lot of functions use this thing
//and it might not behave the way they want when there are Ns in the alignment.
double seqDist_Snn(int segsites, int index1,  int index2, char **list){
	
	int i;
	double compareCount, nCount;
	double count = 0.0;
	char c1, c2;
	
	compareCount = 0.0;
	nCount = 0.0;
	
	for(i = 0; i < segsites; i++){
		c1 = list[index1][i];
		c2 = list[index2][i];
		if(c1 == 'N' || c1 == 'N'){
			nCount += 1; //uncertainty about states?
		}
		else{
			if(c1 != c2)
				count += 1;
		}
		compareCount +=1;
	}
	//arbitrary coverage requirement? disabled for now
	if (nCount / compareCount > 1.0){
		//return large negative number
		return(-666.0);
	}
	return(count);
}

//meanRefDist-- calculates the mean dist of all alleles to a ref--sequence 1
double meanRefDist(int segsites, int nsam, char **list){
	int i;
	double sum = 0.0;
	double denom = 0.0;
	double tmp;
	
	for(i = 1; i < nsam; i++){
		tmp = seqDist_Snn(segsites,0,i,list);
		if (tmp >= 0.0){
			sum += tmp;
			denom++;
		}
	}
	return(sum/denom);
}


//Pairwise IBS stuff
double *pairwiseIBSVec1Popn(int segsites, int n1, int *vecLen, double *posit, char **list){
        double tmpLen,start;
        (*vecLen)=0;
	int i, j, k;
	double *vec;

	vec = (double *) malloc(sizeof(double) * (segsites+1)*n1*(n1-1)/2);
        for(i=0; i<n1-1;i++){
                for(j=i+1;j<n1;j++){
                        start=0.0;
                        //now iterate across sites
                        for(k=0;k<segsites;k++){
                                if(list[i][k] != list[j][k] && list[i][k] != 'N' && list[j][k] != 'N'){
                                        tmpLen = posit[k]-start;
                                        vec[(*vecLen)++] = tmpLen;
                                        start = posit[k];
                                }
                        }
                        tmpLen = 1.0-start;
                        vec[(*vecLen)++] = tmpLen;
                }
        }
	return vec;
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

void pairwiseIBSHist2Popn(int segsites,int nsam, int n1, double *hist, int bins, double *posit, char **list){
	double binWidth = 1.0 / ((float) bins);
	int i, j,k,sumComp=0;
	double tmpLen,start;
	
	for(i=0;i<bins;i++) hist[i] = 0.0;
	for( i=0; i<n1;i++){
		for(j=n1;j<nsam;j++){
			start=0.0;
			//now iterate across sites
			for(k=0;k<segsites;k++){
                                if(list[i][k] != list[j][k] && list[i][k] != 'N' && list[j][k] != 'N'){
					tmpLen = posit[k]-start;
					hist[(int)round(tmpLen/binWidth)]+= 1;
					sumComp +=1;
					start = posit[k];
				}
			}
			tmpLen = 1.0-start;
			hist[(int)round(tmpLen/binWidth)]+= 1;
			sumComp +=1;
		}
	}

	for (i=0;i<bins;i++) hist[i] = hist[i]/sumComp;
}

double pairwiseIBSMax2Popn(int segsites,int nsam, int n1, double *posit, char **list){
	int i, j,k;
	double tmpLen,start, max;
	
	max = 0.0;
	for( i=0; i<n1;i++){
		for(j=n1;j<nsam;j++){
			start=0.0;
			//now iterate across sites
			for(k=0;k<segsites;k++){
                                if(list[i][k] != list[j][k] && list[i][k] != 'N' && list[j][k] != 'N'){
					tmpLen = posit[k]-start;
					//printf("%d %d %d %f\n", i, j, k, tmpLen);
					if(tmpLen > max){
						max = tmpLen;
					}
					start = posit[k];
				}
			}
			tmpLen = 1.0-start;
			//printf("%d %d %d %f\n", i, j, k, tmpLen);
			if(tmpLen > max){
				max = tmpLen;
			}
		}
	}
	return max;
}

double pairwiseIBSMeanWithin(int segsites,int first, int last, double *posit, char **list){
	int i, j,k;
	double tmpLenSum,start;
	float comp =0.0;
	tmpLenSum =0.0;
	for( i=first; i<last-1;i++){
		for(j=i+1;j<last;j++){
			start=0.0;
			//now iterate across sites
			for(k=0;k<segsites;k++){
                                if(list[i][k] != list[j][k] && list[i][k] != 'N' && list[j][k] != 'N'){
					tmpLenSum += posit[k]-start;
					comp += 1.0;
					start = posit[k];
				}
			}
			tmpLenSum += 1.0-start;
			comp += 1.0;
		}
	}
	return tmpLenSum/comp;
}

//From Hudson's ms package
/************************* tajima.c *************************************************************
 This program calculates Tajima's D when number of sequences, number of segregating sites,
   and average pairwise differences (pi) are known.  It also reports all the coefficients for Tajima's
   D (a1, a2, b1, b2, c1, c2, e1, e2). 
**************************************************************************************************/


	double
tajd(int nsam, int segsites, double sumk)
{

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

double a1f(int nsam)
{
double a1;
int i;
a1 = 0.0;
for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
return (a1);
}


double a2f(int nsam) 
{
double a2;
int i;
a2 = 0.0;
for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
return (a2);
}


double b1f(int nsam)
{
double b1;
b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
return (b1);
}


double b2f(int nsam)
{
double b2;
b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
return (b2);
}


double e1f(double a1, double c1) 
{
double e1;
e1 = c1/a1;
return (e1);
}

double e2f(double a1, double a2, double c2)
{ 
double e2;
e2 = c2/((a1*a1)+a2);
return (e2);
}


double c1f(double a1, double b1)
{
double c1;
c1 = b1 - (1/a1);
return (c1);
}


double c2f(int nsam, double a1, double a2, double b2)
{
double c2;
c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
return (c2);
}
