# _cecast_: CoEstimating Contamination And Split Time

The README.md file provides detailed instructions for using _cecast_, a maximum likelihood tool to estimate human contamination and split time in low-coverage archaic samples. It outlines the required scripts, necessary files, installation process, input generation, and execution steps. The document includes example workflows, explanations for different parameter options, and tips for interpreting the results.


## Table of contents
* [General info](#general-info)
* [Necessary scripts](#necessary-scripts)
* [Necessary files](#necessary-files)
* [Installation](#installation)
* [Generate input](#generate-input)
* [Run estimate](#run-estimate)



## General info
_cecast_ is a maximum likelihood tool designed to estimate human contamination and split times in low-coverage archaic samples. It calculates likelihoods using coalescent simulations, generating expectations for ancestral and derived sites across all possible site classes. These simulations incorporate high-coverage genomes, including Mbuti (m), Yoruba (y), Altai Neandertal (a), Vindija Neandertal (v), Chagyrskaya Neandertal (c), and Denisova (d). 


## Necessary scripts

- **`bam2cecast.py`**: Prepares input files from BAM files and can generate pileups with various options at informative sites. It requires some standard Python3 libraries (specified within the script). This script also includes additional useful functions.

- **`vcf2cecast.py`**: (Optional) Prepares input files from VCF files and performs other auxiliary tasks.

- **`cecast.R`**: Executes the estimation process using R. 

- **`install.R`**: Installs necessary R packages and ensures all files are properly loaded when running `cecast.R`.

## Other necessary files 

- **`cecast_funtions.R`**: A set of functions to be loaded for the analyses.

- **`expected_counts_fine_grid.RData`**: Contains simulated expectations for likelihood calculations.

- **`expected_counts_boot.RData`**: Includes expectations for likelihood calculations using 100 bootstrap replicates of the demography.

- **`expected_counts_fine_grid_archaics_afr_combined.RData`**: Simulated expectations using SNPs captured based on combined variations in Africans and Archaics.

- **`expected_counts_fine_grid_archaics_var.RData`**: Expectations using SNPs captured based on variations found in Archaics.

- **`infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz`**: A BED file containing informative sites and necessary details to create input files for cecast. Ensure the file is indexed with `<tabix -p bed infosites*.gz>` before running the `bam2cecast.py`.
   - This BED file (~9.9 million sites) contains informative positions required for generating input files, such as lineage assignments.
      ```bash
      zcat infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz | head -n 10
      #Chr  Start   End     Anc  Der  Atype_lin    Dai  Fre  Han  Man  Mbu  Pap  San  Sar  Yor  Kar
      1     754162  754163  C    T    CBGO_avcdmy  T    C    T    T    T    T    T    T    T    T
      1     755211  755212  T    G    CBGO_avcdmy  G    G    G    G    G    G    G    ?    G    G
      1     755219  755220  G    A    CBGO_avcdmy  A    A    ?    A    A    ?    A    ?    A    A
      1     755238  755239  G    T    CBG_avcdmy   T    ?    ?    T    T    ?    T    ?    T    T
      1     769885  769886  T    C    CBGO_avcdmy  ?    C    C    ?    C    C    ?    C    C    ?
      1     776806  776807  G    A    CBGO_avcdmy  A    A    ?    A    A    A    A    ?    A    A
      1     778112  778113  C    A    BGO_avcdmy   A    A    A    A    A    A    A    A    A    A
      1     779452  779453  G    A    CGO_avcdmy   A    A    A    A    A    A    A    A    A    ?
      1     780033  780034  C    A    CGO_avcdmy   A    A    A    A    A    A    A    A    A    A
      ```
      - *Columns 1-3*: Chromosome, start, and end coordinates.
      - *Columns 4-5*: Ancestral and derived alleles.
      - *Column 6*: Specifies the ancestry type and lineage. For example, *CBGO_avcdmy* means that Chimpanzee, Bonobo, Gorilla and Orangutan have the same ancestral allele, and that Altai, Vindija, Chagyrskaya, Denisova, Mbuti and Yoruba the derived allele. *BGO_acy* means that Bonobo, Gorilla and Orangutan have the same ancestral allele and Altai, Chagyrskaya, and Yoruba the eerived allele. 
      - *Remaining columns*: Represent alleles from 10 modern humans.   
   This structured format is essential for processing data and generating input files effectively.


## Installation
To simplify calling the scripts, create symbolic links to your executable path. For example, you can use your `~/bin` folder, but any location works. Optionally, link `cecast.R` as `cecast` for easier access. Ensure the scripts are executable.
```bash
ln -s scripts/bam2cecast.py ~/bin/bam2cecast.py
ln -s scripts/cecast.R ~/bin/cecast
chmod +x ~/bin/bam2cecast.py ~/bin/cecast
```

Second, run `install.R` in the directory where cecast was downloaded or cloned. This script will download the necessary files, scripts, and R packages, placing them into the appropriate subfolders. If the required components already exist, the script will do nothing.
```bash
R < install.R --no-save  
```
During installation, a hidden JSON file (`~/.config_cecast.json`) is created in your home directory. To place it elsewhere, modify cecast.R accordingly. If reinstallation is needed, delete the existing JSON file to allow the process to reset.

Keep in mind, that you might have to re-index the downloaded BED and BAM files, with `tabix` and `samtools` since the index file might be younger after the download. 

## Generate input
This example demonstrates the cecast pipeline using 100,000 sequences from the Mezmaskaya1 Neanderthal BAM file, downloadable [here](https://bioinf.eva.mpg.de/SpAl/downloads/example.bam). 

```bash
cd data
infosites=infosites_MYAVCD_manifesto_anc3out4CBGO_dist50bp_some_humans.bed.gz
bam2cecast.py -s ${infosites} -l 30 -L 100 -IRES -m 25 -b 10 example.bam > example_ires.tabs
bam2cecast.py -s ${infosites} -l 30 -L 100 -IREF -m 25 -b 10 example.bam > example_iref.tabs
bam2cecast.py -s ${infosites} -l 30 -L 100 -IREF -m 25 -b 10 -d 3,3 example.bam > example_iref_d3.tabs
```
### Description
The resulting `*.tabs` files (inputs for `cecast`) contain lineage assignment tables divided into 10 Mb blocks (specified by -b). These blocks support bootstrapping and confidence interval for later calculations. In the examples, there are only four blocks since the BAM file contains only 100,000 reads on chromosome 1. 

### Key Options:
  - `-I`: Removes reads with indels.
  - `-R`: Randomly samples one read per site.
  - `-E`: Excludes reads neither ancestral nor derived (based on the _infosites_ file).
  - `-S`: Performs strand-specific orientation sampling. Use `-F` instead of `-S` for looser sampling. *Note*: These filters apply only to single-strand libraries.
  - `-m`: Sets the minimum mapping quality. Can be omitted if prefiltered for mappability by length using [_MapL_](https://bioinf.eva.mpg.de/MapL/), desfribed in [de Filippo et al. 2018](https://rdcu.be/d0vIf).
  - `-d`: Filters deaminated reads based on terminal positions. For example, `-d 3,3` considers reads with C-to-T changes in the last three positions at both ends. You can customize the terminal lengths, e.g., `-d 2,4`.
  - `-l` and `-L`: Set the minimum and maximum read lengths, respectively. `-L 100` is optional, as it’s the default.

### Speeding Up the Process
For repeated processing of the same BAM file, filter it to include only reads overlapping informative sites:
```bash
samtools view -bh input.bam -L ${infosites} > input_infosites.bam
```

## Run estimate

Run _cecast_ with default paramters:
```bash
cecast data/example_ires.tabs | column -t
#t  t_low  t_high  pop  c       c_low   c_high  logLike    c_source  nseq  chr   boot                  comments
90  75     90      vin  0.0664  0.0498  0.0952  -131.5178  Yoruba    1108  1:22  data:vin=54%,cha=46%  57_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
```

### Output description
* _**t**_, _**t_low**_ and _**t_up**_: Estimate, lower and upper intervals of split time (in kya) using 4x29 as scaling coalescent unit (29 are the years per generation).
* _**pop**_: The lineage or population from which the test sample splits.
* _**c**_, _**c_low**_ and _**c_up**_: Estimate, lower and upper intervals of human contamination.
* _**logLike**_: Logarithm of the maximum likelihood values (point estimate).
* _**c_source**_: Source of contamination (specified by the `-p` option, Yoruba as default). 
* _**nseq**_: Number of sequences/reads. 
* _**chr**_: Chromosome(s) used. 
* _**boot**_: Type of bootstrap (`-b` option), representing uncertainity from data (as default), history or both. Also reports the fraction of bootstraps supporting different split populations.
* _**comments**_: Indicates other likelihood values that fail the likelihood ratio test (difference < 3.4) and provides bootstrap fractions of different populations split. 

### Run cecast for All Examples
To run _cecast_ for all the examples generated above and store the results in one output file `results.tsv`, use the following bash loop:
```bash
OUT=results.tsv
INDECES='ires iref iref_d3'
echo -e "#file\t$(cecast -H)" > ${OUT} # create the header of the output file
for i in ${INDECES} ; do 
   echo -e "${i}\t$(cecast data/example_${i}.tabs | sed 1d)";
done >> ${OUT}
# Print the results
column -t results.tsv
#file    t   t_low  t_high  pop  c       c_low   c_high  logLike    c_source  nseq  chr   boot                         comments
ires     90  84     90      cha  0.0681  0.0512  0.1512  -131.5273  Yoruba    1108  1:22  data:vin=45%,cha=55%         55_estimates_with_logLike<3.4:vin=51%,cha=15%,v-c=35%
iref     90  84     101     cha  0.058   0.0418  0.136   -140.3437  Yoruba    1454  1:22  data:vin=29%,cha=67%,v-c=4%  50_estimates_with_logLike<3.4:vin=44%,cha=16%,v-c=40%
iref_d3  68  58     71      vin  0.2089  0.1829  0.2968  -110.7259  Yoruba    397   1:22  data:vin=93%,cha=7%          73_estimates_with_logLike<3.4:vin=55%,cha=15%,v-c=30%
```

#### Some explanations
It is surprising that when using deaminated reads, contamination is much higher than when using all reads. This is likely due to stochasticity of the small number of reads (<400). Notice the drastic effect of the strict strand filter (`-S` option in `bam2cecast.py`) on the number of sequences/reads, compared to the looser one (`-F` option in `bam2cecast.py`), which in this case does not affect much the estimates. Keep in mind that you cannot compare the likelihoods across these three examples, as they are highly dependent on the number of sequences: the higher the number of sites, the lower the likelihood.

### _cecast_ Key Options:
   * `-p`: Source of human contamination (default: _Yoruba_). Other options include nine additional human samples. Example: `-p French` or `-p Fr` for French sample. 
   * `-r`: Reference human sample (default: _Mbuti_) used for lineage assignment. _Yoruba_ is another option, but it may not yield reliable results if it is also used as the contamination source. Example: `-c y`. 
   * `-t`: Constraint on the split time range to save time (though the method is already fast). This should be done only after a first estimate to ensure the constraint is within plausible ranges. Example: `-t 70,100`.
   * `-c`: Constraint on contamination as for split times (default: _0,1_). Esample: `-c 0,0.2`. 
   * `-o`: Type of outgroup used to determined ancestry (default: use at least 3 out of 4 outgroup genomes). To be more strict and confindent on the ancestry, specify `-o CBGO` but this reduce the amount of sites (by ~17% in the examples above). 
   * `-O`: Outputs all estimates and bootstraps if you want to look deeper at the results.   

 
### Example source of contamination.
The following bash command loops over the different sources and stores the different estimates (more comprehensible) in `results_c_source.tsv`:

```bash
OUT=results_c_source.tsv
SAMPLES='Dai French Han Mandenka Mbuti Papuan San Sardinian Karitiana'
cecast data/example_iref.tabs > ${OUT}
for s in ${SAMPLES};
   do cecast -p ${s} data/example_iref.tabs | sed 1d;
done >> ${OUT}

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

In this example, Papuan emerges as the best source of contamination, likely due to chance, given that only ~1,500 reads are involved. The log-likelihood differences are minimal. This serves as a simplified illustration, using a limited dataset, to explain the procedure. More extensive examples with larger datasets will demonstrate scenarios where the choice of contamination source significantly impacts results.

 