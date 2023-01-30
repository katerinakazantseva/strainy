[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

# metaPhase

metaPhase is a graph-based phasing algorithm, that takes a de novo assembly graph (in gfa format) and simplifies it by combining phasing information and graph structure.


<img width="694" alt="Screenshot 2023-01-30 at 16 47 16" src="https://user-images.githubusercontent.com/82141791/215481164-2b23544f-589d-4cd2-83f9-a6668ecb8ca6.png">


**metaPhase.py phase** - performs reads clustering according to SNP positions using community detection approach

**metaPhase.py transfom** - transforms assembly graph 

## Installation
### Requirements
```
python3
pip
samtools
```

### Installation

**Ubuntu 20.04:**
```
$ sudo apt update && sudo apt install libbz2-dev liblzma-dev graphviz graphviz-dev bcftools
$ pip install -r requirements.txt
```

**MacOS 12:**
```
$ brew install bcftools
$ pip install -r requirements.txt
```

## Input requirements

metaPhase takes as input gfa graph (can be produced with [**metaFlye**](https://github.com/fenderglass/Flye) or minigraph), 
fasta file and BAM (reads alignned to fasta reference). Also it supports hifi and nanopore modes.

Usage:
```
metaphase.py [-h] [-s SNP] [-t THREADS] -o OUTPUT -b BAM -g GFA -f FA -m {hifi,nano} stage

positional arguments:
  stage                 stage to run: either phase or transform

optional arguments:
  -h, --help            show this help message and exit
  -s SNP, --snp SNP     vcf file
  -t THREADS, --threads THREADS
                        number of threads

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output dir
  -b BAM, --bam BAM     bam file
  -g GFA, --gfa GFA     gfa file
  -f FA, --fa FA        fa file
  -m {hifi,nano}, --mode {hifi,nano}
```

Please specify Flye path in `params.py`:
- `flye` - path to the installed Flye executable

It is not recommended to other parameters.

## Run and outputs

```
python3 ./metaPhase.py phase -o output_dir -b bam_file -g gfa_graph -f fasta file -m mode -t threads
```
Phasing stage clusters reads and produce csv files with read names and corresponding cluster names and BAM file wich visualise reads clustering

<img width="853" alt="Screenshot 2023-01-30 at 16 36 31" src="https://user-images.githubusercontent.com/82141791/215479038-2e0cbba0-0b90-4a12-b84e-5d1965ca193b.png">
<img width="444" alt="Screenshot 2023-01-30 at 16 38 55" src="https://user-images.githubusercontent.com/82141791/215479426-890859d9-75cb-45ad-86ac-f4c41a05ee80.png">



```
python3 ./metaPhase.py transform -o output_dir -b bam_file -g gfa_graph -f fasta file -m mode -t threads
```
Transform stage transform and simplify initial assembly graph, produce  final gfa file transformed_after_simplification_merged.gfa
<img width="903" alt="Screenshot 2023-01-30 at 16 45 20" src="https://user-images.githubusercontent.com/82141791/215480788-3b895736-c43e-43db-a820-6f46c3216a81.png">

## Credits

metaPhase was originally developed at at [**Kolmogorov lab at NCI**](https://ccr.cancer.gov/staff-directory/mikhail-kolmogorov)  


## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg



