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

- **`bam2cecast.py`**: A script that processes BAM files to prepare input files for downstream analysis, including generating pileups at informative sites. It supports a range of options and features useful functions for filtering and formatting data. Dependencies include Python3 standard modules (requiring no additional installation) such as `os`, `copy`, `random`, `itertools`, `signal`, `os.path`, `argparse`, and `textwrap` and the package `pysam`. Ensure all required packages are installed for optimal functionality.

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
   - This BED file (~9.9 million sites) contains informative positions required for generating input files, such as lineage assignments. This structured format is essential for processing data and generating input files effectively.
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
      - *Columns 1-3*: Chromosome, start, and end coordinates for hg19/GCRh37.
      - *Columns 4-5*: Ancestral and derived alleles.
      - *Column 6*: Specifies the ancestry type and lineage. For example, *CBGO_avcdmy* means that Chimpanzee, Bonobo, Gorilla and Orangutan have the same ancestral allele, and that Altai, Vindija, Chagyrskaya, Denisova, Mbuti and Yoruba are all derived. *BGO_acy* means that Bonobo, Gorilla and Orangutan have the same ancestral allele and Altai, Chagyrskaya, and Yoruba the derived allele. 
      - *Remaining columns*: Represent alleles from 10 modern humans.   


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
This example demonstrates the _cecast_ pipeline using 100,000 sequences from the Mezmaskaya1 Neanderthal BAM file, downloadable [here](https://bioinf.eva.mpg.de/SpAl/downloads/example.bam). 

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
  - `-S`: Performs strand-specific orientation sampling: selects reverse-oriented sequences at C-variable sites and forward-oriented sequences at G-variable sites. 
  - `-F`: Applies a less stringent strand-specific sampling: selects reverse and forward orientations at C/T and G/A variable sites, respectively.  
  **Note**: `-S` and `-F` are for single-strand libraries only.
  - `-m`: Sets the minimum mapping quality. Can be omitted if prefiltered for mappability by length using [_MapL_](https://bioinf.eva.mpg.de/MapL/), desfribed in [de Filippo et al. 2018](https://rdcu.be/d0vIf).
  - `-d`: Select deaminated reads based on terminal positions (e.g., `-d 3,3` considers reads with C-to-T changes in the last three positions of both ends).
  - `-l`/`-L`: Set the minimum/maximum read lengths (`-L 100` is the default).

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
90  84     90      cha  0.0692  0.0552  0.1158  -119.4305  Mbuti     1066  1:22  data:vin=37%,cha=63%  56_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
```

### Output description
* _**t**_, _**t_low**_ and _**t_up**_: Estimate, lower and upper intervals of split time (in kya) using 4x29 as scaling coalescent unit (29 are the years per generation).
* _**pop**_: The lineage or population from which the test sample splits.
* _**c**_, _**c_low**_ and _**c_up**_: Estimate, lower and upper intervals of human contamination.
* _**logLike**_: Logarithm of the maximum likelihood values (point estimate).
* _**c_source**_: Source of contamination (specified by the `-p` option, Mbuti as default). 
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
#file    t   t_low  t_high  pop  c       c_low   c_high  logLike    c_source  nseq  chr   boot                          comments
ires     90  84     90      cha  0.0692  0.0527  0.1278  -119.4305  Mbuti     1066  1:22  data:vin=42%,cha=58%          56_estimates_with_logLike<3.4:vin=54%,cha=14%,v-c=32%
iref     90  86     100     cha  0.0523  0.0346  0.0844  -129.0306  Mbuti     1424  1:22  data:vin=12%,cha=71%,v-c=17%  50_estimates_with_logLike<3.4:vin=44%,cha=16%,v-c=40%
iref_d3  68  56     77      vin  0.199   0.1891  0.2663  -108.4103  Mbuti     382   1:22  data:vin=93%,cha=7%           74_estimates_with_logLike<3.4:vin=54%,cha=15%,v-c=31%
```

#### Some explanations
It is surprising that when using deaminated reads, contamination is much higher than when using all reads. This is likely due to stochasticity of the small number of reads (<400). Notice the drastic effect of the strict strand filter (`-S` option in `bam2cecast.py`) on the number of sequences/reads, compared to the looser one (`-F` option in `bam2cecast.py`), which in this case does not affect much the estimates. Keep in mind that you cannot compare the likelihoods across these three examples, as they are highly dependent on the number of sequences: the higher the number of sites, the lower the likelihood.

### _cecast_ Key Options:
   * `-p`: Source of human contamination (default: _Mbuti_). Other options include nine additional human samples. Example: `-p French` or `-p Fr` for French sample. 
   * `-r`: Reference human sample (default: _Yoruba_) used for lineage assignment. _Mbuti_ is another option, but it may not yield reliable results if it is also used as the contamination source. Example: `-r m`. 
   * `-t`: Constraint on the split time range to save time (though the method is already fast). This should be done only after a first estimate to ensure the constraint is within plausible ranges. Example: `-t 70,100`.
   * `-c`: Constraint on contamination as for split times (default: _0,1_). Esample: `-c 0,0.2`. 
   * `-o`: Type of outgroup used to determined ancestry (default: use at least 3 out of 4 outgroup genomes). To be more strict and confindent on the ancestry, specify `-o CBGO` but this reduce the amount of sites (by ~17% in the examples above). 
   * `-O`: Outputs all estimates and bootstraps if you want to look deeper at the results.   

 
### Example source of contamination.
The following bash command loops over the different sources and stores the different estimates (in a more comprehensible format) in `results_c_source.tsv`. Before running the loop, I execute the estimate using Yoruba as both the reference (`-r y`) and the source of contamination, which corresponds to the same individual.

```bash
OUT=results_c_source.tsv
# Change the reference to Mbuti since Yoruba is the contamination source
cecast -p Yoruba -r m data/example_iref.tabs > ${OUT}
SAMPLES='Mbuti Dai French Han Mandenka Papuan San Sardinian Karitiana'
# Loop over the different sources with Yoruba (default) as reference used in the lineage assignment.
for s in ${SAMPLES};
   do cecast -p ${s} data/example_iref.tabs | sed 1d;
done >> ${OUT}
# Print the results


column -t ${OUT}
#t  t_low  t_high  pop  c       c_low   c_high  logLike    c_source   nseq  chr   boot                          comments
90  84     101     cha  0.058   0.0367  0.1223  -140.3437  Yoruba     1454  1:22  data:vin=34%,cha=60%,v-c=6%   50_estimates_with_logLike<3.4:vin=44%,cha=16%,v-c=40%
90  86     101     cha  0.0523  0.0358  0.074   -129.0306  Mbuti      1424  1:22  data:vin=22%,cha=61%,v-c=17%  50_estimates_with_logLike<3.4:vin=44%,cha=16%,v-c=40%
90  84     101     cha  0.0522  0.0413  0.0928  -129.2257  Dai        1424  1:22  data:vin=14%,cha=71%,v-c=15%  50_estimates_with_logLike<3.4:vin=42%,cha=16%,v-c=42%
90  86     101     cha  0.0525  0.0398  0.0789  -128.9354  French     1424  1:22  data:vin=11%,cha=71%,v-c=18%  50_estimates_with_logLike<3.4:vin=42%,cha=16%,v-c=42%
90  86     101     cha  0.0525  0.0372  0.0955  -129.0684  Han        1424  1:22  data:vin=7%,cha=78%,v-c=15%   49_estimates_with_logLike<3.4:vin=43%,cha=16%,v-c=41%
90  86     101     cha  0.0521  0.0381  0.0687  -128.9366  Mandenka   1424  1:22  data:vin=13%,cha=75%,v-c=12%  49_estimates_with_logLike<3.4:vin=43%,cha=16%,v-c=41%
90  86     101     cha  0.0546  0.0353  0.0879  -128.9745  Papuan     1424  1:22  data:vin=14%,cha=73%,v-c=13%  49_estimates_with_logLike<3.4:vin=43%,cha=16%,v-c=41%
90  86     93      cha  0.0546  0.0376  0.079   -128.9185  San        1424  1:22  data:vin=22%,cha=72%,v-c=6%   51_estimates_with_logLike<3.4:vin=45%,cha=16%,v-c=39%
90  84     101     cha  0.0518  0.0348  0.0887  -128.9004  Sardinian  1424  1:22  data:vin=11%,cha=76%,v-c=13%  50_estimates_with_logLike<3.4:vin=44%,cha=16%,v-c=40%
90  86     101     cha  0.0525  0.0364  0.0931  -129.2951  Karitiana  1424  1:22  data:vin=17%,cha=68%,v-c=15%  50_estimates_with_logLike<3.4:vin=42%,cha=16%,v-c=42%
```

First, note that when using Yoruba as the source of contamination and Mbuti as the reference population (first row), the number of reads (column **nseq**) differs from the other estimates. As a result, the likelihoods are not directly comparable. Excluding this Yoruba estimate, Sardinian appears to be the best source of contamination. However, this result could be due to chance, given the limited dataset (~1,500 reads). While the log-likelihood differences are small, larger datasets would more clearly illustrate the effect of the contamination source selection. This example serves as a simplified demonstration of the procedure.


#### Effect of using the same sample as reference and source of contamination.
The following command generates estimates using the same population as both the reference and contamination source, specifically Yoruba and Mbuti, which are currently the only samples used as references for lineage assignment.

```bash
OUT=results_c_source_ref.tsv
echo -e "#Ref\t$(cecast -H)" > ${OUT}
# Mbuti as reference and source of contamination
echo -e "Mbuti\t$(cecast -r m data/example_iref.tabs | sed 1d)" >> ${OUT}
# Mbuti as reference and Yoruba as source of contamination
echo -e "Mbuti\t$(cecast -r m -p Yoruba data/example_iref.tabs | sed 1d)" >> ${OUT}
# Yoruba as reference and source of contamination
echo -e "Yoruba\t$(cecast -p Yoruba data/example_iref.tabs | sed 1d)" >> ${OUT}
# Yoruba as reference and Mbuti as source of contamination
echo -e "Yoruba\t$(cecast -r y -p Mbuti data/example_iref.tabs | sed 1d)" >> ${OUT}

column -t ${OUT}
#Ref    t   t_low  t_high  pop  c       c_low   c_high  logLike    c_source  nseq  chr   boot                         comments
Mbuti   90  86     101     cha  0.0499  0.0315  0.0808  -141.2587  Mbuti     1454  1:22  data:vin=26%,cha=70%,v-c=4%  50_estimates_with_logLike<3.4:vin=42%,cha=16%,v-c=42%
Yoruba  93  86     103     v-c  0.0305  0.024   0.0522  -131.4266  Yoruba    1424  1:22  data
```

Using the same sample as both the reference and the contamination source typically results in a lower contamination estimate. However, this is not entirely the case for Mbuti, where the contamination estimate is only 5% lower. This occurs because the reference and contamination samples, while from the same population, are not identical. In contrast, using Yoruba results in almost half the contamination estimate (0.030 vs. 0.058) since the same individual is used for both the reference and contamination source. Keep this in mind when using the `-p` and `-r` options. 


### Examples with more sequences from Mezmaskaya 1
The following examples use the same sample Mezmaskaya but with many more reads using this input `data/mez1_L34MQ1_iref.tabs`. This has been generated with `bam2cecast.py` using a much bigger BAM file, which I omitted from this repository for size reason. 

#### Less contamination when selecting for deaminated sequences
In the example before, selecting for deaminated sequences using `the example.bam` resulted in more contamination, but this is not the case when using much more data, as our expecations of aDNA. 

```bash
cat <(cecast data/mez1_iref.tabs ) <(cecast data/mez1_iref_d3.tabs | sed 1d) | column -t
#t  t_low  t_high  pop  c       c_low   c_high  logLike     c_source  nseq    chr   boot  comments
93  86     93      v-c  0.0233  0.0219  0.0261  -1423.9877  Mbuti     300897  1:22  data
93  85     93      v-c  0.0046  0.0028  0.0078  -690.9816   Mbuti     92120   1:22  data
```

#### Compare the likelihoods
The following bash command does the same as before with different contamination sources. 

```bash
OUT=results_c_source_mez1.tsv
# For Yoruba as source of contamination change the reference to Mbuti since Yoruba is the default reference.
cecast -p Yoruba -r m data/mez1_iref.tabs > ${OUT}
# Loop for the other sources
SAMPLES='Mbuti Dai French Han Mandenka Papuan San Sardinian Karitiana'
for s in ${SAMPLES};
   do cecast -p ${s} data/mez1_iref.tabs | sed 1d;
done >> ${OUT}

column -t ${OUT}
#t  t_low  t_high  pop  c       c_low   c_high  logLike     c_source   nseq    chr   boot  comments
93  86     93      v-c  0.0223  0.0205  0.0245  -1793.7847  Yoruba     305005  1:22  data
93  88     93      v-c  0.0233  0.0218  0.0261  -1423.9877  Mbuti      300897  1:22  data
93  88     93      v-c  0.0223  0.021   0.0248  -1425.9039  Dai        300897  1:22  data
93  88     93      v-c  0.0224  0.0208  0.0253  -1425.5827  French     300897  1:22  data
93  86     93      v-c  0.0225  0.0208  0.0248  -1425.3675  Han        300897  1:22  data
93  88     93      v-c  0.0215  0.0201  0.024   -1428.7057  Mandenka   300897  1:22  data
93  88     93      v-c  0.0228  0.0212  0.0254  -1423.2131  Papuan     300897  1:22  data
93  85     93      v-c  0.0245  0.0234  0.0278  -1423.8176  San        300897  1:22  data
93  86     93      v-c  0.0223  0.021   0.0252  -1426.1349  Sardinian  300897  1:22  data
93  86     93      v-c  0.0225  0.0207  0.0251  -1425.4836  Karitiana  300897  1:22  data
```


