[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

# stRainy

stRainy is a graph-based phasing algorithm, that takes a de novo assembly graph (in gfa format) and simplifies it by combining phasing information and graph structure.

<p align="center">
<img width="694" alt="Screenshot 2023-01-30 at 16 47 16" src="https://user-images.githubusercontent.com/82141791/215481164-2b23544f-589d-4cd2-83f9-a6668ecb8ca6.png">
</p>

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
./strainy.py transform  -o out_strainy -b test_set/toy.bam -g test_set/toy.gfa -t 4 -m hifi 
```

## Limitations

stRainy is under active development! The current version is optimized for a relatively simple 
bacterial communities (one or a few bacterial species, 2-5 strains each). Extending
stRainy to larger metagenomes is a work in progress.

## Input requirements

stRainy supports PacBio HiFi and Nanopore (Guppy5+) sequencing. 

The inputs to stRainy are:
1. GFA file (can be produced with [**metaFlye**](https://github.com/fenderglass/Flye) or minigraph) 
and
2. BAM file (reads aligned to the fasta reference generated from the GFA file).

How to get fasta from gfa and perform alignment (assuming ONT reads):
```
awk '/^S/{print ">"$2"\n"$3}' assembly_graph.gfa > assembly_graph.fasta
minimap2 -ax map-ont assembly_graph.fasta reads.fastq | samtools sort -@4 -t 8 > assembly_graph.bam
samtools index assembly_graph.bam
```

## Improving de novo metagenomic assmbelies

We have tested stRainy using metaFlye metagenoimic assembly graphs as input. The recommended
set of parameters is `--meta --keep-haplotypes -i 0`. 

Note that `-i 0` disables metaFlye's polishing procedure, which we found to improve read assignemnt
to bubble branches during `minimap2` realignment. `--keep-haplotypes` retains structural
variations between strains on the assmebly graph.

## Usage and outputs

**strainy.py phase** - performs reads clustering according to SNP positions using community detection approach

**strainy.py transfom** - transforms assembly graph 

```
./strainy.py phase -o output_dir -b bam_file -g gfa_graph -m mode -t threads
```

Phase stage clusters reads and produces csv files with read names, corresponding cluster names and a BAM file. The BAM file visualises the clustering of the reads

<p align="center">
<img width="500" alt="Screenshot 2023-01-30 at 17 01 47" src="https://user-images.githubusercontent.com/82141791/215484889-6a032cc0-9c90-4a26-9689-7d5cb41a2ab5.png">
</p>

```
./strainy.py transform -o output_dir -b bam_file -g gfa_graph -m mode -t threads
```

Transform stage transforms and simplifies initial assembly graph, producing the final gfa file: `transformed_after_simplification_merged.gfa`

<p align="center">
<img width="500" alt="Screenshot 2023-01-30 at 16 45 20" src="https://user-images.githubusercontent.com/82141791/215480788-3b895736-c43e-43db-a820-6f46c3216a81.png">
</p>

## Parameters desciption

```
usage: strainy.py [-h] -o OUTPUT -b BAM -g GFA -m {hifi,nano} [-s SNP] [-t THREADS] [-f FASTA] [-splu] [-fq FASTQ]
                  [-splen UNITIG_SPLIT_LENGTH]
                  {phase,transform}

positional arguments:
  {phase,transform}     stage to run: either phase or transform

options:
  -h, --help            show this help message and exit
  -s SNP, --snp SNP     vcf file (default: None)
  -t THREADS, --threads THREADS
                        number of threads to use (default: 4)
  -f FASTA, --fasta FASTA
                        fasta file (default: None)
  -splu, --split-long-unitigs
                        Split unitigs with long sequences into multiple smaller unitigs for faster processing, requires -fq to be
                        provided (default: False)
  -fq FASTQ, --fastq FASTQ
                        fastq file containing reads to perform alignment, only used if -splu is set (default: False)
  -splen UNITIG_SPLIT_LENGTH, --unitig-split-length UNITIG_SPLIT_LENGTH
                        The length (in kb) which the unitigs that are longer will be split, only used if -splu is set (default: 50)

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output dir (default: None)
  -b BAM, --bam BAM     bam file (default: None)
  -g GFA, --gfa GFA     gfa file (default: None)
  -m {hifi,nano}, --mode {hifi,nano}
```

## Acknowledgements

Consensus function of stRainy is [**Flye**](https://github.com/fenderglass/Flye)

Community detection algorithm is [**Karate club**](https://github.com/benedekrozemberczki/KarateClub/blob/master/docs/source/notes/introduction.rst)

## Credits

stRainy was originally developed at at [**Kolmogorov lab at NCI**](https://ccr.cancer.gov/staff-directory/mikhail-kolmogorov)  

Code contributors:
- Ekaterina Kazantseva
- Ataberk Donmez
- Mikhail Kolmogorov

## Citation

Ekaterina Kazantseva, Ataberk Donmez, Mihai Pop, Mikhail Kolmogorov.
"stRainy: assembly-based metagenomic strain phasing using long reads"
bioRxiv 2023, https://doi.org/10.1101/2023.01.31.526521

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
