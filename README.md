# _cecast_
CoEstimating Contamination And Split Time


## Table of contents
* [General info](#general-info)
* [Necessary scripts](#necessary-scripts)
* [Necessary files](#necessary-files)
* [Generate input](#generate-input)
* [Run estimate](#run-estimate)
* [List of files](#list-of-files)


## General info
cecast is a maximum likelihood estimate of human contamination and split time in low-coverage archaic samples.
The likelihood is calculated using coalescent simulations that generates the expectations for the number of ancestral and derived sites for all possibile class of sites, using _M-Y-A-V-C-D_ high-coverage genomes, where _M_ is modern human Mbuti, _Y_ is modern human Yoruba, _A_ is Altai Neandertal, _V_ is Vindija Neandertal, _C_ is Chagyrskaya Neandertal and _D_ is Denisova. The simulations also include modern human, European which experienced 2.5% gene flow from Neandertals around 70 kya.  


## Necessary scripts

* `bam2cecast.py` prepare the input file from bam file. It can also generate a pileup with quite a few options at informative sites.
   (some standard python libraries are required, but you can check the script.)

* `vcf2cecast.py` prepare the input file from vcf file. It also does something else. 

*  `cecast.R` runs the estimate with R and requires the following libraries: argparse, parallel, RColorBrewer, fields, and jsonlite.   


## Other necessary files 

*  `cecast_funtions.R` functions to be loaded.

*  `infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz` bed file with the informative sites and the infos necessary to create input files for cecast.
   The file should be indexed with ```tabix -p bed infosites*.gz``` before running the `bam2cecast.py`.

This is an example of the bed file with the informative sites (~9.9 millions) which is necessary to generate the input file, like the lineage assignment. Columns 1,2 and 3 are the chromosome, start and end coordinates, columns 4 and 5 the ancestral and derived alleles, column 6 specify the ancestry type and the lineage, and the rest of the columns are the alleles in 10 modern humans.  

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
I describe the pipeline by using an example of 100,000 sequences from Mezmaskaya1 bam file that can be downloaded from [here](https://bioinf.eva.mpg.de/SpAl/downloads/example.bam): 

```
infosites=infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz
bam2cecast.py -l 30 -L 100 -IRES -s ${infosites} -m 25 -b 10 example.bam > example_L30MQ25.tabs
bam2cecast.py -l 30 -L 100 -IREF -s ${infosites} -m 25 -b 10 example.bam > example_L30MQ25_sf_loose.tabs
bam2cecast.py -l 30 -L 100 -IRES -d 3,3 -s ${infosites} -m 25 -b 10 example.bam > example_L30MQ25_deam3.tabs
```

The `-IRES` option means that it filters for indels (`-I`), randomly samples one read (`-R`), excludes reads that are neither ancestral or derived in the INFO file (`-E`) and does the strand-specific orientation sampling (`-S`). Option `-m` is the mappping quality, which can be skept if the mappability by length filter, desfribed in [de Filippo et al. 2018](https://bioinf.eva.mpg.de/MapL/), is applied before. The input will be a lineage assignment like table but in blocks of 10Mb specified by the `-b` option.


## Run estimate
```
cecast data/example_L30MQ25.tabs 
#t t_low t_high pop c      c_low  c_high logLike 	c_source nseq chr  boot                 comments 	
90 75    90     vin 0.0664	0.0486 0.0933 -131.5178 Yoruba   1108 1:22 data:vin=52%,cha=48% 57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
```

with the columns of the output:
* _**t**_, _**t_low**_ and _**t_up**_ are the estimate, the lower and upper intervals of split time in kya using 4x29 as scaling coalescent unit (29 are the years per generation).
* _**pop**_ is the lineage or population from which the test sample splits.
* _**c**_, _**c_low**_ and _**c_up**_ are the estimate, the lower and upper intervals of human contamination.
* _**logLike**_ is the log of the maximum likelihood values which is the point estimate.
* _**c_source**_ is the source of contamination specified by the otion `-p`. 
* _**nseq**_ is the number of sequences/reads. 
* _**chr**_ is the chromosome(s) used. 
* _**boot**_ is the type of bootstrap (`-b` option) generated from uncertainity in data, history and both. Eventually, it reports the fraction of bootstraps supporting different split populations. 
* _**comments**_ reports whether there are other likelihhod values that do not pass the likelihood ratio test, i.e. when the difference between the max and these values is lower than 3.4. As for _boot_ column, it reports the fraction of different populations split.  

Here, follows a description of some options but look at the help of cecast to see all options.

* `-p` is the source of human contamination (dfault _Yoruba_). There are also other nine human samples from the SGDP. 
* `-r` is the reference human sample (default _Mbuti_) used for the lineage assignment. The other option is _Yoruba_ but notice that it will not produce reliable results if _Yoruba_ is also the source of contamination. 


It is possible to constrain on the split time range with the option `-t`, in order to save some time although this gain is marginal since the method run already quite fast (less than a minute). This should be done only after a first estimate so that the constrain is within the plausible ranges, which in this example would be between 70 and 100 kya. 

```
cecast -p AncientCtrl -t 60,150 lin_den4.tabs
#t 	t_low 	t_high 	pop 	c 	c_low 	c_high 	logLike 	c_source	nseq	boot	comments 	
92	72	132	den	0.665	0.6334	0.685	-195.2035	AncientCtrl	16456	data
```


It is possible to change the source of human contamination with `-p`. There are samples from the SGDP (see help message for the full list).


Here it is an excample of using different source of contamination. The following bash command loop over the different sources and store the different estimates (more friendly) in results_c_source.tsv.    
```
OUT=results_c_source.tsv
SAMPLES='Dai French Han Mandenka Mbuti Papuan San Sardinian Karitiana'
cecast data/example_L30MQ25.tabs > ${OUT}
for s in ${SAMPLES} ; do cecast -p ${s} data/example_L30MQ25.tabs | sed 1d; done >> ${OUT}

column -t ${OUT}
#t  t_low  t_high  pop  c       c_low   c_high  logLike    c_source   nseq  chr   boot                  comments
90  75     90      vin  0.0664  0.0498  0.0952  -131.5178  Yoruba     1108  1:22  data:vin=57%,cha=43%  57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
90  75     90      vin  0.0659  0.0501  0.1061  -131.5083  Sardinian  1108  1:22  data:vin=69%,cha=31%  57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
90  75     90      vin  0.0664  0.0525  0.0978  -131.4173  French     1108  1:22  data:vin=56%,cha=44%  57_estimates_with_logLike<3.4:vin=53%,cha=14%,v-c=33%
90  84     90      cha  0.0707  0.0526  0.1682  -131.254   Han        1108  1:22  data:vin=44%,cha=56%  57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
90  75     90      vin  0.0671  0.0556  0.1018  -131.0246  Mandenka   1108  1:22  data:vin=50%,cha=50%  57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
90  75     90      vin  0.0571  0.0439  0.0997  -132.8482  Mbuti      1108  1:22  data:vin=63%,cha=37%  57_estimates_with_logLike<3.4:vin=53%,cha=14%,v-c=33%
90  84     90      cha  0.0746  0.0588  0.1256  -130.5146  Papuan     1108  1:22  data:vin=32%,cha=68%  56_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
90  75     90      vin  0.0659  0.0501  0.1061  -131.5083  Sardinian  1108  1:22  data:vin=54%,cha=46%  57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
90  75     90      vin  0.0659  0.0501  0.1259  -131.5083  Sardinian  1108  1:22  data:vin=49%,cha=51%  57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
90  75     90      vin  0.0691  0.0604  0.1037  -131.7344  Karitiana  1108  1:22  data:vin=42%,cha=58%  57_estimates_with_logLike<3.4:vin=53%,cha=14%,v-c=33%

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

