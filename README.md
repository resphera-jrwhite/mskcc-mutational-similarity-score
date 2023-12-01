# Mutational Similarity Score (mskcc-mutational-similarity-score)
### Diaz Laboratory | MSKCC
This is a version-controlled codebase for the mutational similarity score project
under the Diaz Lab at Memorial Sloan Kettering Cancer Center. 

### Introduction
The purpose of the code is to utilize shared mutation observations between two input samples to estimate the 
likelihood of sharing the same origin. The primary script `A01-mutational-profile-similarity.pl` can be run 
from the command line by providing two single-column input text files with the mutations observed for the
two samples to be compared. The input file mutation format is based on [TCGA MAF formatted coordinates](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format): `<gene name>`\_`<chr#>`\_`<beginning coordinate>`\-`<ending coordinate>`\_`<ref allele>`\_`<alt allele>` for example: `BRAF_chr7_140453136-140453136_A_T`. 

You will also need to select a cancer set and panel type to serve as the background database. See the `Usage` sections below for more information.


### Repository Structure
* Top-level directory: `mskcc-mutational-similarity-score/`
```
./mskcc-mutational-similarity-score/
├── A01-mutational-profile-similarity.pl     <-- main code to run similarity assessment
├── B01-example_data_run.sh                  <-- example runs using samples from Chaudhuri et al. 2017
├── README.md                                <-- this readme file
├── example-msi-mss-results                  <-- MPS results summary of MSS and MSI CRC / UCEC samples
└── supporting-materials         
    ├── background_dbs                       <-- reference lists of mutations (.gz) and sample types (.txt)
    │   ├── COADMSS.samples.txt
    │   ├── ...
    │   └── stage2.output.PS64.txt.gz
    └── example_data                         <-- example data from Chaudhuri et al. 2017 (PMID: 28899864)
        ├── LUP14_TimePt0Tumor.txt
        ├── ...
        └── LUP20_TimePt3.txt
```

### Requirements
This code requires Perl (v5.30) and a few other modules to be installed in order to run (e.g. via cpan):
* File::Spec
* Data::Dumper
* List::Util
* Getopt::Std

### Prior to use: uncompress gz files
```
cd mskcc-mutational-similarity-score
gunzip supporting-materials/background_dbs/*gz
```

### Primary Usage
```
Usage:  ./A01-mutational-profile-similarity.pl
        -i file with Sample1 mutations (mutations formatted as: BRAF_chr7_140453136-140453136_A_T)
        -j file with Sample2 mutations (mutations formatted as: BRAF_chr7_140453136-140453136_A_T)
        -t cancer type (see below for options)
        -p panel ID (ASX, PS64, CS125, CS447, MSKIMPACT341, MSKIMPACT410, MSKIMPACT468, CAPPSeq, WEX, Sample1)
        -e list of TCGA ids to exclude from empirical analysis (for internal validation)

.OPTIONS.
  Available cancer types: ALLTCGA,     ACC,    BLCA,  BRCA,  CESC,  CHOL,  COAD,  DLBC,
                             ESCA,     GBM,    HNSC,  KICH,  KIRC,  KIRP,   LGG,  LIHC,
                             LUAD,    LUSC,      OV,  PAAD,  PCPG,  PRAD,  READ,  SARC,
                             SKCM,    STAD,    TGCT,  THCA,  THYM,  UCEC,   UCS,   UVM, COADREAD,

                          COADMSS,   UCECMSS, COADREADMSS,
                          DIAZMSI, MOSAICMSI, DIAZMSIPOLE

  Available panel IDs:    ASX, PS64, CS125, CS447, MSKIMPACT341, MSKIMPACT410,
                          MSKIMPACT468, CAPPSeq, WEX, Sample1

.DESCRIPTION.
  This code is designed to use recurrent mutations mined from TCGA (MC3)
  to predict whether two tumor samples are related in nature. The code
  computes a probability that two random tumor samples overlap by chance.
  The code adjusts for mutations that are highly recurrent as well, e.g.
  the BRAF V600E mutation in Thyroid cancer is present in approximately 60%
  of patients, so it's much more likely to occur by chance than a different
  more unique mutation.

  Mutation coordinates are in HG19 space.

  If the user chooses "Sample1" as the input panel ID, then the ROIs are limited to those specific mutations
  listed in the Sample1 input file (simulates ddPCR results).
```

### Example Usage (outputs to STDOUT)
```
./A01-mutational-profile-similarity.pl \
    -i supporting-materials/example_data/LUP14_TimePt0Tumor.txt \
    -j supporting-materials/example_data/LUP14_TimePt1.txt \
    -p ASX \
    -t LUAD
```

### Example Output 
```
STATUS: Computing background frequencies of each Shared mutation...
STATUS: Empirical assessment comparing Sample1 to negative controls...

EMPIRICAL RESULTS: 0 out of 462 greater than or equal to observed value of 2 Shared mutations...

BEGIN REPORT:
Input Cancer Type: LUAD (Lung_adenocarcinoma)
Input Panel ID: ASX
Sample1: supporting-materials/example_data/LUP14_TimePt0Tumor.txt
Sample2: supporting-materials/example_data/LUP14_TimePt1.txt
Negative control samples: 462
Total mutations in Sample1: 5
Total mutations in Sample2: 2
Total mutations Shared: 2

SHARED MUTATION STATISTICS:
Shared_Mutation	General_Frequency_in_LUAD(%)	95%CI
EGFR_chr7_55259515-55259515_T_G	4.545%  	(2.906%,6.977%)
NF1_chr17_29553484-29553484_C_T	0.216%  	(0.011%,1.394%)

NEGATIVE CONTROL EMPIRICAL RESULTS:
TYPE1 ERROR ESTIMATE: 0.216% 95%CI:(0.011%,1.394%)
OVERALL PROBABILITY OF RELATEDNESS: 99.784% 95%CI:(98.606%,99.989%)
```

### Example Multisample Runs
```
./B01-example_data_run.sh
```


## [Issues and Requests](https://github.com/resphera-jrwhite/mskcc-mutational-similarity-score/issues)
For stakeholders interested in requests for data transfers, methods developement, analyses, documentation, etc, please visit the [issues](https://github.com/resphera-jrwhite/mskcc-mutational-similarity-score/issues) tab to submit a ticket. Please select an associated label, mention users in the description and assign to a username when possible.

## Authors

* **James Robert White PhD** - *Initial work* - [GitHub](https://github.com/resphera-jrwhite)


