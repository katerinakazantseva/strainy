[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

# stRainy

stRainy is a graph-based phasing algorithm, that takes a de novo assembly graph (in gfa format) and simplifies it by combining phasing information and graph structure.


<img width="694" alt="Screenshot 2023-01-30 at 16 47 16" src="https://user-images.githubusercontent.com/82141791/215481164-2b23544f-589d-4cd2-83f9-a6668ecb8ca6.png">


**strainy.py phase** - performs reads clustering according to SNP positions using community detection approach

**strainy.py transfom** - transforms assembly graph 

## Conda Installation


The recommended way of installing is though [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```
git clone https://github.com/katerinakazantseva/stRainy
cd stRainy
git submodule update --init
make -C submodules/Flye
conda env create -f environment.yml -n strainy
```

Note that if you use an M1 conda installation, you should run `conda config --add subdirs osx-64` before installation. 
Find details [**here**](https://github.com/conda/conda/issues/11216)

Once installed, you will need to activate the conda environemnt prior running:

```
conda activate strainy
./strainy.py
```

## Quick usage example

After successful installation, you should be able to run:

```
conda activate strainy
./strainy.py phase  -o out_strainy -b test_set/toy.bam -g test_set/toy.gfa -t 4 -m hifi 
./strainy.py transfrom  -o out_strainy -b test_set/toy.bam -g test_set/toy.gfa -t 4 -m hifi 
```

## Input requirements

stRainy can support both hifi or nanopore models. The inputs to metaPhase are:
1. GFA file (can be produced with [**metaFlye**](https://github.com/fenderglass/Flye) or minigraph) 
and
2. BAM file (reads aligned to the fasta reference generated from the GFA file).


How to get fasta from gfa and perform alignment (assuming ONT reads):
```
awk '/^S/{print ">"$2"\n"$3}’ assembly_graph.gfa > assembly_graph.fasta
minimap2 -ax map-ont assembly_graph.gfa reads.fastq | samtools sort -@4 -t 8 > assembly_graph.bam
samtools index assmebly_grpah.bam
```

Usage:
```
strainy.py [-h] [-s SNP] [-t THREADS] [-f FASTA] -o OUTPUT -b BAM -g GFA -m {hifi,nano} stage

positional arguments:
  stage                 stage to run: either phase or transform

optional arguments:
  -h, --help            show this help message and exit
  -s SNP, --snp SNP     vcf file
  -t THREADS, --threads THREADS
                        number of threads
  -f FASTA, --fasta FASTA
                        fasta file

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output dir
  -b BAM, --bam BAM     bam file
  -g GFA, --gfa GFA     gfa file
  -m {hifi,nano}, --mode {hifi,nano}
```

## Run and outputs

```
./strainy.py phase -o output_dir -b bam_file -g gfa_graph -m mode -t threads
```
Phasing stage clusters reads and produce csv files with read names and corresponding cluster names and BAM file wich visualise reads clustering

<img width="500" alt="Screenshot 2023-01-30 at 17 01 47" src="https://user-images.githubusercontent.com/82141791/215484889-6a032cc0-9c90-4a26-9689-7d5cb41a2ab5.png">



```
./strainy.py transform -o output_dir -b bam_file -g gfa_graph -m mode -t threads
```
Transform stage transform and simplify initial assembly graph, produce  final gfa file transformed_after_simplification_merged.gfa
<img width="500" alt="Screenshot 2023-01-30 at 16 45 20" src="https://user-images.githubusercontent.com/82141791/215480788-3b895736-c43e-43db-a820-6f46c3216a81.png">

## Asknowledgements

Consesus function of stRainy is [**Flye**](https://github.com/fenderglass/Flye)

Community detection algorithm is [**Karate club**](https://github.com/benedekrozemberczki/KarateClub/blob/master/docs/source/notes/introduction.rst)


## Authors

stRainy was originally developed at at [**Kolmogorov lab at NCI**](https://ccr.cancer.gov/staff-directory/mikhail-kolmogorov)  

Code contributors:
- Ekaterina Kazantseva
- Ataberk Donmez
- Mikhail Kolmogorov


## Citation

TBD

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg



