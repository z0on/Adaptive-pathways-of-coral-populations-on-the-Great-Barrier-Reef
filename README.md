Scripts and data to accompany this preprint:
M. V. Matz, M. J. H. van Oppen, L. K. Bay and E. A. Treml (2017) Adaptive pathways of coral populations on the Great Barrier Reef. bioRxiv:

Summary

We genotyped five populations of coral Acropora millepora along the Great Barrier Reef using 2bRAD (http://ecogeno.weebly.com/uploads/7/6/2/2/76229469/wang12_2b-rad.pdf) and inferred population sizes and migration rates using dadi (https://bitbucket.org/gutenkunstlab/dadi). These results were compared to biophysical model of larval transport and used to build a multilocus metapopulation adaptation model in SLiM (https://messerlab.org/slim/). 

Data generation

2bRAD reads were processed using scirpts and pipeline described at https://github.com/z0on/2bRAD_GATK. The cleaned reads were mapped to the genome of Acropora digitifera and quality-filtered (using GATK's Variant Quality Score Recalibration procedure, http://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr) based on true SNPs identified as a match between genotyping replicates for several individuals. 

The final accuracy (match between replicates) was >98% at >95% genotyping rate (fraction of non-missing data). Heterozygote discovery rate (proportion of matching heterozygote calls among replicates) was >95%.

Main genotyping results in VCF format (25433 biallelic SNPs):
Amil_mapped2adig_cleaned.vcf.zip
In the sample names the initial letter identifies the source population:
	W: Wilkie Island
	S: Sudbury reef
	O: Orpheus Island
	M: Magnetic Island (Nelli Bay)
	K: North Keppel Island

Data analysis overview (see detailed_walkthrough_amilGBR.txt for details)

1. Excluding potential sites under selection. 
The VCF file was reformatted into Bayescan format using PDGspider and Bayescan was fun on the whole dataset. Thirteen loci with q-value less than 0.05 were excluded from the subsequent analysis.

2. Pairwise Fst and ADMIXTURE.
Pairwise Weir and Cockerham Fst between populations were calculated using vcftools. The dataset was thinned to include only SNPs on average 5 kb 9and at least 2.5 kb) apart, choosing SNPs with highest alternative allele frequency. This dataset contained and leaving 11426 SNPs. It was reformatted into bed format (plink) and fed to ADMIXTURE. 

3. dadi
The full dataset was thinned again as described in the previous paragraph but this time the SNPs were chosen without regard to allele frequency, not to distort the allele frequency spectrum. The thinned VCF was reformatted into dadi format and 120 bootstrap replicates of it were created by resampling A.digitifera genome contigs with replacement. The data were analyzed by a series of 1d and 2d dadi models to determine effective historical population sizes and pairwise migration rates. The models were compared by summarizing AIC differences across bootstrap replicates.

4. multi-QTL metapopulation model in SLiM.
The model was developed based on the SLiM recipe “Quantitative genetics and phenotypically-based fitness”. The model simulates Fisher-Wright populations with discreet generations. At the start of the simulation, populations are established at specified population sizes and pairwise migration rates (genetic replacement rates), and all QTLs in all individuals are given a mutation with the effect size drawn from a normal distribution with mean zero and specified standard deviation, to create standing genetic variation. The phenotype of each individual is calculated as the sum of QTL effects plus random noise to simulate desired heritability. Then, fitness of each individual is calculated based on the difference between the individual’s phenotype (thermal optimum), temperature of the environment, and the setting for phenotypic plasticity, modeled as the standard deviation of the Gaussian slope of fitness decline with increasing distance between phenotype and environment. Then, parents are chosen to produce the next generation according to their fitness; parents for immigrant individuals are chosen from among individuals in the source population. New mutations at QTLs happened at the specified rate when transitioning to the next generation and the effect of a new mutation replaced the previous QTL effect.

Environment is modeled as identical trends across populations with constant population-specific offsets. The trends can be any combination of linear, cyclical and random components. 

To better model population dynamics, we implemented linear scaling of the population size and immigration rates with the population’s mean fitness post pre-adaptation period (the initial 2000 generations with no linear change in environment). In this way, a population declining in fitness shrinks in size and stops contributing migrants to other populations

--------
SCRIPTS:
--------

Accessory scripts (run without arguments to see full usage description):
thinner.pl : for thinning a VCF file
vcf2dadi.pl : to convert VCF to dadi format
dadiBoot.pl	: generates bootstrapped dadi datasets by resampling entries in the "Gene" column
removeBayescanOutliers.pl : removes SNPs with Bayescan qvalue less than specified.


Command-line dadi scripts:
projections_calc_auto.py : calculates number of segregating site for different projection values
1d_spectrum_plot.py : plots 1d frequency spectrum to pdf
2d_spectrum_plot.py	: plots 1d frequency spectrum to pdf
onegrowth_auto.py :	one-population model with a single growth period
twogrowth_auto.py :	one-population model with a two growth periods
s2m_auto.py	: two-population s2m model with split into two different sizes with asymmetrical migration
twogrowth2d_null.py : two-population model with size change but no split - a null model to compare to s2m to prove that populations are demographically distinct. 
s2mRM_auto.py : extension of s2m model with a change in migration rates in the last 0.01 time units.

SLiM:
SLIM_WSOMKmodel_withNeutral.txt : model main code
slimPlots.R	: script to plot fitness and phenotype trends from SLiM output
