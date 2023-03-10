[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

# stRainy

stRainy is a graph-based phasing algorithm, that takes a de novo assembly graph (in gfa format) and simplifies it by combining phasing information and graph structure.

<p align="center">
<img width="694" alt="Screenshot 2023-01-30 at 16 47 16" src="https://user-images.githubusercontent.com/82141791/215481164-2b23544f-589d-4cd2-83f9-a6668ecb8ca6.png">
</p>

## Conda Installation

The recommended way of installing is through [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```
git clone https://github.com/katerinakazantseva/stRainy
cd stRainy
git submodule update --init
make -C submodules/Flye
conda env create -f environment.yml -n strainy
```

Note that if you use an M1 conda installation, you should run `conda config --add subdirs osx-64` before installation. 
Find details [**here**](https://github.com/conda/conda/issues/11216)

Once installed, you will need to activate the conda environment prior running:

```
conda activate strainy
./strainy.py -h
```

## Quick usage example

After successful installation, you should be able to run:

```
conda activate strainy
./strainy.py -g test_set/toy.gfa -q test_set/toy.fastq.gz -o out_strainy -m hifi 
```

## Limitations

stRainy is under active development! The current version is optimized for a relatively simple 
bacterial communities (one or a few bacterial species, 2-5 strains each). Extending
stRainy to larger metagenomes is a work in progress.

## Input requirements

stRainy supports PacBio HiFi and Nanopore (Guppy5+) sequencing. 

The two main inputs to stRainy are:
1. GFA file (can be produced with [**metaFlye**](https://github.com/fenderglass/Flye) or minigraph) 
and
2. FASTQ file (containing reads to be aligned to the fasta reference generated from the GFA file).

## Improving de novo metagenomic assmbelies

We have tested stRainy using metaFlye metagenoimic assembly graphs as input. The recommended
set of parameters is `--meta --keep-haplotypes -i 0`. 

Note that `-i 0` disables metaFlye's polishing procedure, which we found to improve read assignemnt
to bubble branches during `minimap2` realignment. `--keep-haplotypes` retains structural
variations between strains on the assmebly graph.

## Usage and outputs
stRainy has 2 stages: **phase** and **transform**. With the command below, stRainy will phase and transform by default. Please see Parameter Description section for the full list of available arguments:

```
./strainy.py -g [gfa_file] -q [fastq_file] -m [mode] -o [output_dir]
```  

 **1. phase** stage performs read clustering according to SNP positions using community detection approach and produces csv files with read names, corresponding cluster names and a BAM file. The BAM file visualises the clustering of the reads.

<p align="center">
<img width="500" alt="Screenshot 2023-01-30 at 17 01 47" src="https://user-images.githubusercontent.com/82141791/215484889-6a032cc0-9c90-4a26-9689-7d5cb41a2ab5.png">
</p>

<br>

**2. transform** stage transforms and simplifies the initial assembly graph, producing the final gfa file: `strainy_final.gfa`

<p align="center">
<img width="500" alt="Screenshot 2023-01-30 at 16 45 20" src="https://user-images.githubusercontent.com/82141791/215480788-3b895736-c43e-43db-a820-6f46c3216a81.png">
</p>

## Parameter desciption

```
usage: strainy.py [-h] -o OUTPUT -g GFA -m {hifi,nano} -q FASTQ [-stage {phase,transform,e2e}] [-s SNP] [-t THREADS] [-f FASTA] [-b BAM] [--unitig-split-length UNITIG_SPLIT_LENGTH]

options:
  -h, --help            show this help message and exit
  -stage {phase,transform,e2e}
                        stage to run: either phase, transform or e2e (phase + transform) (default: e2e)
  -s SNP, --snp SNP     vcf file (default: None)
  -t THREADS, --threads THREADS
                        number of threads to use (default: 4)
  -f FASTA, --fasta FASTA
                        fasta file (default: None)
  -b BAM, --bam BAM     bam file (default: None)
  --unitig-split-length UNITIG_SPLIT_LENGTH
                        The length (in kb) which the unitigs that are longer will be split, set 0 to disable (default: 50)

Required named arguments:
  -o OUTPUT, --output OUTPUT
                        directory that will contain the output files (default: None)
  -g GFA, --gfa GFA     gfa file (default: None)
  -m {hifi,nano}, --mode {hifi,nano}
                        type of reads (default: None)
  -q FASTQ, --fastq FASTQ
                        fastq file containing reads to perform alignment, used to create a .bam file (default: None)
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
