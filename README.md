# _cecast_
CoEstimating Contamination And Split Time


## Table of contents
* [General info](#general-info)
* [Necessary scripts](#necessary-scripts)
* [Necessary files](#necessary-files)
* [Installation](#installation)
* [Generate input](#generate-input)
* [Run estimate](#run-estimate)



## General info
_cecast_ is a maximum likelihood estimate of human contamination and split time in low-coverage archaic samples.
The likelihood is calculated using coalescent simulations that generates the expectations for the number of ancestral and derived sites for all possibile class of sites, using _M-Y-A-V-C-D_ high-coverage genomes, where _M_ is modern human Mbuti, _Y_ is modern human Yoruba, _A_ is Altai Neandertal, _V_ is Vindija Neandertal, _C_ is Chagyrskaya Neandertal and _D_ is Denisova. The simulations also include modern human, European which experienced 2.5% gene flow from Neandertals around 70 kya.  


## Necessary scripts

- **bam2cecast.py** prepare the input file from bam file. It can also generate a pileup with quite a few options at informative sites.
   (some standard python libraries are required, look within the script and install them.) The scripts also does other useful things.  

- **vcf2cecast.py** (optional) prepare the input file from vcf file. It also does something else. 

- **cecast.R** runs the estimate with R. 

- **install.R** will download the necessary R-packages and make sure that the all files will be properly load when running `cecast.R`. 

## Other necessary files 

- **cecast_funtions.R** functions to be loaded.

- **expected_counts_fine_grid.RData** the expectations from simulations to calculate the likelihood.

- **expected_counts_boot.RData** the expectations from simulations to calculate the likelihood using 100 bootstraps of the demography. 

- **expected_counts_fine_grid_archaics_afr_combined.RData** the expectations from simulations using ascertain SNPs captured according to the combined variations found in Africans and Archaics. 

- **expected_counts_fine_grid_archaics_var.RData** the expectations from simulations using ascertain SNPs captured according to the  variations found in Archaics. 

- **infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz** bed file with the informative sites and the infos necessary to create input files for cecast.
   The file should be indexed with `<tabix -p bed infosites*.gz>` before running the `bam2cecast.py`.

This is an example of the bed file with the informative sites (~9.9 millions) which is necessary to generate the input file, like the lineage assignment. Columns 1,2 and 3 are the chromosome, start and end coordinates, columns 4 and 5 the ancestral and derived alleles, column 6 specify the ancestry's type and the lineage, and the rest of the columns are the alleles in 10 modern humans.  

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

## Installation
First, create symbolic links to your executable path in order to call the scripts easily. I have them in my home 'bin' folder but you can do whatever you want. Notice that I link `cecast.R` to just `cecast` for an even easier calling and don't forget to make them executable. 
```
ln -s scripts/bam2cecast.py ~/bin/bam2cecast.py
ln -s scripts/cecast.R ~/bin/cecast
```

Second, run `install.R` within the directory where cecast is download/cloned. It will download the necessary files/scripts/R-packages and make sure they are in the proper sub-folders. If you have them already, it will do nothing. However, this will create a hidden jason file `~/.config_cecast.json` in your home directory, but if you want to place it somewhere else, you have to hack a bit `cecast.R` for loading the necessary files/scripts. Also notice that if the file `*json` is already present the installation will not do anything as well. Therefore, if something is messed up and you want to re-do the installation remove the file `*json`.

```
R < install.R --no-save  
```

## Generate input
I describe the pipeline by using an example of 100,000 sequences from Mezmaskaya1 bam file that can be downloaded from [here](https://bioinf.eva.mpg.de/SpAl/downloads/example.bam). 

```
cd data
infosites=infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz
bam2cecast.py -s ${infosites} -l 30 -L 100 -IRES -m 25 -b 10 example.bam > example_ires.tabs
bam2cecast.py -s ${infosites} -l 30 -L 100 -IREF -m 25 -b 10 example.bam > example_iref.tabs
bam2cecast.py -s ${infosites} -l 30 -L 100 -IREF -m 25 -b 10 -d 3,3 example.bam > example_iref_d3.tabs
```

The output *.tabs files above (or input for cecast) will be a lineage assignment table in blocks of 10Mb specified by the `-b` option. These blocks are used to produce bootstraps and calculate confidence intervals. 
The `-IRES` options: remove reads with indels (`-I`), randomly samples one read at each site (`-R`), excludes reads that are neither ancestral or derived in the _infosites_ file (`-E`) and does the strand-specific orientation sampling (`-S`). One can used `-F` instead of `-S` for a looser strand-specific orientation sampling. Also notice that this filters only work on single-strand libraries. 
Option `-m` is the mappping quality, which can be skept if the mappability by length filter [MapL](https://bioinf.eva.mpg.de/MapL/), desfribed in [de Filippo et al. 2018](https://rdcu.be/d0vIf), is applied before. 
Option `-d` consider only deaminated reads in the last _n_-terminal positions of the 5' and 3' ends. In the example above, `-d 3,3` consider reads with a C-to-T change in the last three-terminal positions of both ends, but other values also with different _n_ for the terminals (e.g. `-d 2,4`) could also be specified. 
Options `-l` and `-L` filter for read minimum and maximum read lengths, respectively. In the example above, `-L 100` can be omitted since it is the default. 

In order to speed up the analyses and if you want to process the bam file(s) many times, you can filter the bam for reads overlapping the informative sites as following:
```
samtools view -bh input.bam -L ${infosites} > input_infosites.bam
```

## Run estimate

Run cecast with default paramters. 

```
cecast data/example_ires.tabs | column -t
#t  t_low  t_high  pop  c       c_low   c_high  logLike    c_source  nseq  chr   boot                  comments
90  75     90      vin  0.0664  0.0498  0.0952  -131.5178  Yoruba    1108  1:22  data:vin=54%,cha=46%  57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
```

The columns of the output:
* _**t**_, _**t_low**_ and _**t_up**_ are the estimate, the lower and upper intervals of split time in kya using 4x29 as scaling coalescent unit (29 are the years per generation).
* _**pop**_ is the lineage or population from which the test sample splits.
* _**c**_, _**c_low**_ and _**c_up**_ are the estimate, the lower and upper intervals of human contamination.
* _**logLike**_ is the log of the maximum likelihood values which is the point estimate.
* _**c_source**_ is the source of contamination specified by the otion `-p`. 
* _**nseq**_ is the number of sequences/reads. 
* _**chr**_ is the chromosome(s) used. 
* _**boot**_ is the type of bootstrap (`-b` option) generated from uncertainity in data, history and both. Eventually, it reports the fraction of bootstraps supporting different split populations. 
* _**comments**_ reports whether there are other likelihhod values that do not pass the likelihood ratio test, i.e. when the difference between the max and these values is lower than 3.4. As for _boot_ column, it reports the fraction of different populations split.  

Run cecast for all examples generated above with a loop in bash and store them in one output file `results.tsv`. The output will have an extra first column which refers to the file/filters used. 
```
OUT=results.tsv
INDECES='ires iref iref_d3'
echo -e "#file\t$(cecast -H)" > ${OUT} # create the header of the output file
for i in ${INDECES} ; do 
   echo -e "${i}\t$(cecast data/example_${i}.tabs | sed 1d)"
done >> ${OUT}
```

Print the results.
```
column -t results.tsv
#file    t   t_low  t_high  pop  c       c_low   c_high  logLike    c_source  nseq  chr   boot                         comments
ires     90  84     90      cha  0.0681  0.0512  0.1512  -131.5273  Yoruba    1108  1:22  data:vin=45%,cha=55%         55_estimates_with_logLike<3.4:vin=51%,cha=15%,v-c=35%
iref     90  84     101     cha  0.058   0.0418  0.136   -140.3437  Yoruba    1454  1:22  data:vin=29%,cha=67%,v-c=4%  50_estimates_with_logLike<3.4:vin=44%,cha=16%,v-c=40%
iref_d3  68  58     71      vin  0.2089  0.1829  0.2968  -110.7259  Yoruba    397   1:22  data:vin=93%,cha=7%          73_estimates_with_logLike<3.4:vin=55%,cha=15%,v-c=30%
```
It is suprising that when using deaminated reads, contamination is much higher than when using all sequences. This is rather a due to stochasticity, and the small number of reads (<400). Notice the drastic effect in terms of number of sequences/reads of the strict strand filter, compared to the looser one, which in this case does not affect much the estimates. Furthermore keep in mind that you cannot compare the likelihoods across these three examples since they are highly dependent on the number of sequences: higher the number of sites, lower the likelihood.     

Here, follows a description of some options but look at the help of cecast to see all options.

* `-p` is the source of human contamination (dfault _Yoruba_). There are also other nine human samples from the SGDP. 
* `-r` is the reference human sample (default _Mbuti_) used for the lineage assignment. The other option is _Yoruba_ but notice that it will not produce reliable results if _Yoruba_ is also the source of contamination. 
* `-t` constrain on the split time range, in order to save some time although this gain is marginal since the method run already quite fast (less than a minute). This should be done only after a first estimate so that the constrain is within the plausible ranges, which in this example would be between 70 and 100 kya. (example `-t 70,100`). 


Here it is an excample of using different source of contamination. The following bash command loop over the different sources and store the different estimates (more friendly) in results_c_source.tsv.    
```
OUT=results_c_source.tsv
SAMPLES='Dai French Han Mandenka Mbuti Papuan San Sardinian Karitiana'
./cecast.R data/example_L30MQ25_sf_loose.tabs > ${OUT}
for s in ${SAMPLES} ; do ./cecast.R -p ${s} data/example_L30MQ25_sf_loose.tabs | sed 1d; done >> ${OUT}

column -t ${OUT}
#t  t_low  t_high  pop  c       c_low   c_high  logLike    c_source   nseq  chr   boot                          comments
90  88     101     vin  0.0529  0.0385  0.0817  -139.6876  Yoruba     1454  1:22  data:vin=56%,cha=34%,v-c=10%  51_estimates_with_logLike<3.4:vin=47%,cha=14%,v-c=39%
90  86     91      cha  0.0547  0.0292  0.1116  -139.6168  Dai        1454  1:22  data:vin=36%,cha=62%,v-c=2%   52_estimates_with_logLike<3.4:vin=46%,cha=15%,v-c=38%
90  84     99      vin  0.0529  0.0398  0.0868  -139.5285  French     1454  1:22  data:vin=45%,cha=44%,v-c=11%  52_estimates_with_logLike<3.4:vin=46%,cha=15%,v-c=38%
90  88     93      vin  0.0548  0.04    0.0862  -139.5552  Han        1454  1:22  data:vin=33%,cha=65%,v-c=2%   52_estimates_with_logLike<3.4:vin=46%,cha=15%,v-c=38%
90  81     98      vin  0.0537  0.0395  0.0883  -139.4142  Mandenka   1454  1:22  data:vin=50%,cha=38%,v-c=12%  52_estimates_with_logLike<3.4:vin=46%,cha=15%,v-c=38%
90  85     101     vin  0.0452  0.0313  0.0808  -140.6318  Mbuti      1454  1:22  data:vin=49%,cha=38%,v-c=13%  52_estimates_with_logLike<3.4:vin=46%,cha=13%,v-c=40%
90  86     93      cha  0.0593  0.0316  0.1165  -139.0503  Papuan     1454  1:22  data:vin=40%,cha=57%,v-c=3%   52_estimates_with_logLike<3.4:vin=46%,cha=15%,v-c=38%
90  81     93      vin  0.0527  0.0393  0.0701  -139.5258  San        1454  1:22  data:vin=48%,cha=45%,v-c=7%   52_estimates_with_logLike<3.4:vin=46%,cha=15%,v-c=38%
90  81     93      vin  0.0538  0.0407  0.0899  -139.3427  Sardinian  1454  1:22  data:vin=50%,cha=41%,v-c=9%   52_estimates_with_logLike<3.4:vin=46%,cha=15%,v-c=38%
90  86     93      cha  0.0548  0.0333  0.0918  -139.7738  Karitiana  1454  1:22  data:vin=38%,cha=57%,v-c=5%   52_estimates_with_logLike<3.4:vin=46%,cha=15%,v-c=38%
```

In this example, Papuan is the best source of contamination, but this is rather due to chance since there are only about 1500 reads. Furthermore, the log likelihoods are not that different. This is only one example, where I used just a few sequences to explain the proceduere, but I will show later on some examples where the source of contamination really matters. 

Notice instead that  