# _ecast_
Maximum likelihood Estimate of human Contamination And Split Time in low-coverage archaic samples.

## Table of contents
* [General info](#general-info)
* [Necessary scripts](#necessary-scripts)
* [Necessary files](#necessary-files)
* [Generate input](#generate-input)
* [Run estimate](#run-estimate)
* [List of files](#list-of-files)


## General info
The likelihood is calculated using coalescent simulations that generates the expectations for the number of ancestral and derived sites for all possibile class of sites, using _M-Y-A-V-C-D_ high-coverage genomes, where _M_ is modern human Mbuti, _Y_ is modern human Yoruba, _A_ is Altai Neandertal, _V_ is Vindija Neandertal, _C_ is Chagyrskaya Neandertal and _D_ is Denisova. The simulations also include modern human, European which experienced 2.5% gene flow from Neandertals around 70 kya.  


## Necessary scripts

* `bam2cecast.py` prepare the input file. It can also generate a pileup with quite a few options at informative sites.
   (some standard python libraries are required, but you can check the script.)

*  `cecast.R` runs the estimate with R and requires the following libraries: argparse, parallel, RColorBrewer and fields   


## Other necessary files 

*  `infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz` bed file with the informative sites and the infos necessary to create input files for cecast.

This is an example of the bed file with the informative sites (~9.9 millions) which is necessary to generate the input file, like the lineage assignment. The file should be indexed with ```tabix -p bed infosites*.gz``` before running the `bam2cecast.py`.

```
zcat infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz | head -n 10
1	754162	754163	C	T	CBGO_avcdmy	T	C	T	T	T	T	T	T	T	T
1	755211	755212	T	G	CBGO_avcdmy	G	G	G	G	G	G	G	?	G	G
1	755219	755220	G	A	CBGO_avcdmy	A	A	?	A	A	?	A	?	A	A
1	755238	755239	G	T	CBG_avcdmy	T	?	?	T	T	?	T	?	T	T
1	769885	769886	T	C	CBGO_avcdmy	?	C	C	?	C	C	?	C	C	?
1	776806	776807	G	A	CBGO_avcdmy	A	A	?	A	A	A	A	?	A	A
1	778112	778113	C	A	BGO_avcdmy	A	A	A	A	A	A	A	A	A	A
1	779452	779453	G	A	CGO_avcdmy	A	A	A	A	A	A	A	A	A	?
1	780033	780034	C	A	CGO_avcdmy	A	A	A	A	A	A	A	A	A	A
1	785502	785503	A	G	CBGO_avcdmy	G	G	G	A	G	G	G	G	G	G
```

## Generate input
I describe the pipeline by using a few examples: 1) Denisova4 and 2) Oase1.  

```
infosites=infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz
BAM=/mnt/sequencedb/PopGen/cesare/ecast/0_bams/Denisova4.chrs.rmdup.L30MQ1.mapL.noIndels.bam
bam2cecast.py -l 32 -L 100 -RISE -d 3,3 -s ${infosites} -m 25 -b 10 ${BAM} > testsample_L32MQ25_deam3.tabs
bam2cecast.py -l 34 -L 100 -RISE -s ${sites} -m ${m} -b 10 ${bam} > ${i}_${r}_L34MQ${m}_w_hum.tabs

```

The `-IRES` option means that it filters for indels (`-I`), randomly samples one read (`-R`), excludes reads that are neither ancestral or derived in the INFO file (`-E`) and does the strand-specific orientation sampling (`-S`). Option `-m` is the mappping quality, which can be skept if the mappability by length filter (desfribed in de Filippo et al. 2018) is applied before. The input will be a lineage assignment like table but in blocks of 10Mb specified by the `-b` option.


## Run estimate
```
cecast -p Ancient lin_den4.tabs
#t 	t_low 	t_high 	pop 	c 	c_low 	c_high 	logLike 	c_source	nseq	boot	comments 	
92	72	132	den	0.66	0.63	0.6852	-195.3644	AncientCtrl	16456	data		
Warning message:
(-p) population: Sample AncientCtrl is taken but Ancient was specified. 
```
where
* `-p Ancient` uses our _Ancient Control_ as source of human contamination. There are also other humans samples from the SGDP, of which Yoruba is the default. Notice that this produces the *Warning message* since the name specified with the `-p` option does not perfectly match. There is also the option to use simulations by specifing `-p sim` (or `sims`, `simulations`, `Simulat`, etc... all case insensitive), but this has lower precision than using real samples since simulations are always (a bit or a lot) different from the truth. 

and with the columns of the output:
* _**t**_, _**t_low**_ and _**t_up**_ are the estimate, the lower and upper intervals of split time in kya using 4x29 as scaling coalescent unit (29 are the years per generation).
* _**pop**_ is the lineage or population from which splits (see genealogy in fig X)
* _**c**_, _**c_low**_ and _**c_up**_ are the estimate, the lower and upper intervals of human contamination.
* _**logLike**_ is the log of the maximum likelihood values which is the point estimate.
* _**c_source**_ is the source of contamination specified by the otion `-p`. 
* _**nseq**_ is the number of sequences/reads. 
* _**chr**_ is the chromosome(s) used. 
* _**boot**_ is the type of bootstrap (`-b` option) generated from uncertainity in data, history and both. Eventually, it reports the fraction of bootstraps supporting different split populations. 
* _**comments**_ reports whether there are other likelihhod values that do not pass the likelihood ratio test, i.e. when the difference between the max and these values is lower than 3.4. As for _boot_ column, it reports the fraction of different populations split.  


In order to have a more precise estimate of contamination, we can change the values of the grid to look for with the option `-c`. The default `-c 0,1,0.01` is to search from 0% to 100% in steps of 1%. This is the argument passed to the R-function _seq()_; therefore, to search in steps of 0.1% use `-c 0,1,0.001`. However, this will slow down the procedure because there are 10 times more values for the contamination to be considered. 
Nevertheless, since the estimates are between 63% and 69%, we could ask to look between 60% and 70% with `-c 0.6,0.7,0.001`. This will be as fast as before: both with 101 total values to consider. 
In order to save some time, it is possible to constrain on the split time range with the option `-t`. This should be done only after a first estimate so that the constrain is within the plausible ranges, which in this example would be between 60 and 150 kya. 

```
cecast -p AncientCtrl -c 0.6,0.7,0.001 -t 60,150 lin_den4.tabs
#t 	t_low 	t_high 	pop 	c 	c_low 	c_high 	logLike 	c_source	nseq	boot	comments 	
92	72	132	den	0.665	0.6334	0.685	-195.2035	AncientCtrl	16456	data
```

Notice that both estimates are slightly different, but the second one is better since the likelihood is slightly higher, although not signifincatly so. This can be done when there is the need of knowing the precise estimate of contamination, for instance to use contamination to correct some other estimates.   

It is possible to change the source of human contamination with `-p`. There are samples from the SGDP (see help message for the full list) and the ancient control. The former are using genotypes and the latter is using reads of short lengths to mimic ancient DNA.
It is also possible to use simulations as source of contamination (i.e. a simulated European). However, in my experience this always leads to underestimates since the simulated source of contamination is more likely to differ from the real one.  
If nothing is specified the default is to use Yoruba as contaminant. 

Here it is an excample of using different source of contamination. The following bash command loop over the different sources and store the different estimates (more friendly) in results_c_source.tsv.    
```
OUT=results_c_source.tsv
SAMPLES='Gambian French Sardinian Japanese Papuan  AncientCtrl Simulations'
cecast lin_den4.tabs > ${OUT}
for s in ${SAMPLES} ; do cecast -p ${s} lin_den4.tabs | sed 1d; done >> ${OUT}

column -t ${OUT}
#t   t_low  t_high  pop  c     c_low   c_high  logLike    c_source       nseq   chr   boot  comments
132  92     180     den  0.62  0.5748  0.65    -221.0021  S_Yoruba-1     16456  1:22  data
132  74     177     den  0.62  0.58    0.66    -219.8928  S_Gambian-2    16456  1:22  data
112  72     147     den  0.64  0.6048  0.67    -216.4896  S_French-2     16456  1:22  data
112  72     140     den  0.64  0.6148  0.67    -215.311   S_Sardinian-2  16456  1:22  data
112  72     147     den  0.64  0.61    0.67    -216.2334  S_Japanese-2   16456  1:22  data
92   72     117     den  0.66  0.64    0.6852  -214.9521  S_Papuan-14    16456  1:22  data
92   72     132     den  0.66  0.6348  0.6852  -195.3644  AncientCtrl    16456  1:22  data
217  182    247     den  0.51  0.47    0.54    -288.4961  simulation     16456  1:22  data
```

By changing the source of contaminant there is some change for both split time (varing of 40%)  and contamination. To determine which is the better estimate, it is sufficient to compare the likelihoods since the estimates uses the same number of sequences. If you filter the data, and therefore generate a different input file, you should not compare the likelihood. 

Using an interesting sample, Oase1. Run the estimantes with different human sources and for different chromosomes. 

```
## SAMPLES and INFO are defined above. 
BAM=/mnt/sequencedb/PopGen/cesare/ecast/0_bams/Oase1_strict_filter_nuclear_varsites.bam 
OUT=results_oase.tsv
IN=lin_oase.tabs
bam2lineage.py -b 10 -s ${INFO} -m 25 -l 30 -SIRE ${BAM} > ${IN}

echo -e "#c_source\tt\tt_low\tt_high\tpop\tc\tc_low\tc_high\tlogLike\tComments" > ${OUT}
cecast ${IN} > ${OUT}

for s in ${SAMPLES} ; do cecast -p ${s} ${IN} | sed 1d; done >> ${OUT}

s=S_Gambian-2
for k in {1..22} X ; do cecast -k ${k} -p ${s} ${IN} | sed 1d ; done >> ${OUT}

```


## List of files

- **cecast.R** and **ecast.R** 

run the estimates. The latter is the older version.

- **cecast_exp_boot.RData** 

cointains the expectations for calculating the likelihood. Old file was **ecast_expected_counts.RData**)

-- simulations

- **run_simulations.R** 

runs the simulations beased on the demographic history estimated by Fabrizio Mafessoni.

- **power_analyses.R** 

perform the power analyses using simulations or real data as test dataset.

- **ecast_functions.R** 

contains all the functions necessary for the project. It is very useful for debugging.

