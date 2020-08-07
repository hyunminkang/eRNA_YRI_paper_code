# eRNA_YRI_paper_code

## Overview

This repository contains the source code and example data to reproduce
the results described in the manuscript entitled
"Population-scale study of eRNA transcription reveals bipartite
functional enhancer architecture"

This repository contains the following two folders:
* `TRE_identification` folder contains the code and examples to identify
  transctibe regulatory elements (TREs).
* `variant_alignment` folder contains the code used to account for
  mappability and reference allele bias using allele-specific
  mapping. 
  
## Prerequisites

To use this repository, the following software tools need to be
installed and should be accessible from your default path.
The software version indicates the version used in our experiments.
Other versions may work, but the same version may be required for
a complete reproducbility.

* `bedtools` (v2.28)
* `bowtie` (v1.2.2)
* `tabix` (v1.7-2)
* `R` (v3.6.1)
* `samtools`(v1.9)
* `git` (v2.24.3)
* `wget` (v1.20.1)

## Getting Started

### Clone the repository

First, clone the repository with the following commands

```
## clone the repository
git clone https://github.com/hyunminkang/eRNA_YRI_paper_code.git
## change the working directory
cd eRNA_YRI_paper_code
```

### Download the example 

Second, download the example data (1.2GB) from our FTP site using the
following commands

```
## your current working directory must be eRNA_YRI_paper_code
sh download.sh
```

The download may take several minutes, so please be patient.

### Review README and run the codes

For each subdirectory, there is a `README.txt` file. Review the file
information carefully to understand the expected behavior.

To run `TRE_identification`, use the following commands:

```
## your current working directory must be eRNA_YRI_paper_code
cd TRE_identiciation
sh sh/make.PROcap.eTSS.sh 
cd .. ## move to the original directory
```

To run `variant_alignment`, use the following commands:

```
## your current working directory must be eRNA_YRI_paper_code
cd variant_alignment
sh sh/make.sh
cd .. ## move to the original directory
```

You may modify the variables to apply the scripts to other samples

## Reference

To cite our repository, please cite the following information:

  Kristjánsdóttir K, Kwak Y, Tippens ND, Lis JT, Kang HM, Kwak H. Population-scale study of eRNA transcription reveals bipartite functional enhancer architecture. bioRxiv. 426908.
