NOTE: LEAVING THIS FOR HISTORICAL PURPOSES. UNDER SOME SCENARIOS CONVERTING FROM 0/1 TO -1/1 IS BENEFICIAL. TYPICALLY WHEN THE MEAN OF EACH MATRIX IS NEAR 0.5. HOWEVER AS ALLELE FREQS ARE COMMONLY SKEWED, THE MEAN WILL BE CLOSER TO ZERO, AND THEN 0/1 ENCODING MAY BE BETTER. SO CHECK YOUR DATA BEFORE RECODING. THE SIMS HERE DO NOT REFLECT 'NORMAL' COALESCENTS, AND SHOULD NOT BE OVERINTERPRETED REGARDING SNP ENCODING

This is a simple net that predicts Tajima's D and demonstrates the impact of different techniques used to transform the training data. 

Figure_1.png contains the results.  The lower the RMSE the better.  The lines are the average of 10 runs.  Note y-axis is log scaled to pick up small diffs at lower end.

The red line is the SNP matrix as it is normally portrayed (indv on rows, SNPs on cols).  This forces convolutions to go indv-wise rather than snp-wise. It works poorly, so we can just forget about that tactic. Then the other 4 colors all utilze a transposed SNP matrix (SNPs on rows), which is clearly a big improvement.

The blue line is just this matrix transposed with input SNPs as 0/1.  Green is transposed and SNPs changed from 0/1 to -1/1.  Black is transposed, SNPs as 0/1, and the indv in the matrix sorted by genetic similarity.   And finally magenta is transposed, resorted, and SNPs -1/1, i.e. the full enchilada.  Magenta is the best, but sorting indv. by genetic sim also does a lot of good.  

I've isolated the functions that do these data transformations and put them in a file called `genet.data.matrix.prep.tricks.py`.
