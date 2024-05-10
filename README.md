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
1. GFA file: A de novo metagenomic assembly that can be produced with [**metaFlye**](https://github.com/fenderglass/Flye) or minigraph).
For metaFlye parameters please check **Improving de novo metagenomic assemblies** below.
You can also try to format fasta into a gfa file in the case of a linear genome (experimentally,not tested)
3. FASTQ file (containing reads to be aligned to the fasta reference generated from the GFA file).

## Improving de novo metagenomic assemblies

We have tested stRainy using metaFlye metagenoimic assembly graphs as input. The recommended
set of parameters is `--meta --keep-haplotypes --no-alt-contigs -i 0`. 

Note that `-i 0` disables metaFlye's polishing procedure, which we found to improve read assignemnt
to bubble branches during `minimap2` realignment. `--keep-haplotypes` retains structural
variations between strains on the assmebly graph. `--no-alt-contigs` disables the output of
"alternative" contigs, that can later confuse the read aligner.

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

## Parameter description

| Key  | Description |
| ------------- | ------------- |
|-o	|directory that will contain the output files (default: None)|
|-g	|assembly gfa file (may be produced with metaFlye or minigraph)|
|-q	|FASTQ file ( PacBio HiFi or  Nanopore sequencing)|
|-m	|type of reads {hifi,nano}|
| -snp (Optional)	|you can provide your own vcf file, with variants of the desired AF. Otherwise stRainy will use the built-in pileup-based caller(default: None)|
| -a (Optional)	|ALLELE_FREQUENCY: allele frequency trhreshhold for built-in pileup-based caller. Will only work if snp=None (default: None)|
| -d (Optional)|	CLUSTER_DIVERGENCE: Cluster sensitivyty-maximum number of total mismatches allowed in the cluster per 1kbp. Should be selected depending on SNP rates and their accuracy. Increasing the parameter can reduce high fragmentation but at the same time reduce clustering accuracy (default: None)|
| -unitig-split-length (Optional)	|stRainy splits the unitigs of the original graph to improve performance so that they do not exceed the threshold (defailt: 50K)|
|-b (Optional)	| BAM: If you want, you can use a prebuilt bam file (aligned reads per assembly), but this will only work in case of --unitig-split-length 0 or (default: None)|
|-s (Optional)	| stage to run: either phase, transform or e2e (phase + transform) (default: e2e)|
|--min-unitig-coverage/--max-unitig-coverage (Optional)	|thresholds to determine which unitigss to phase|
|-threads (Optional)	|number of threads to use (default: 4)|
|-debug (Optional) |	debug mode enable |

## Output description




| Path  | Description | Example |
| ------------- | ------------- | ------------- | 
| strain_contigs.gfa | phased graph (before simplifying links and merging contigs) | <img width="170" alt="image" src="https://github.com/katerinakazantseva/stRainy/assets/82141791/71347e47-d18d-4d35-a3c8-3f36561ae970"> 
| strain_unitigs.gfa	| phased graph (after simplifying links and merging contigs) | <img width="170" alt="image" src="https://github.com/katerinakazantseva/stRainy/assets/82141791/6c254095-d7f7-4104-bc94-f8b1a15a8c07"> 
| strain_variants.vcf |	vcf produced by stRainy build-in caller if not provided by user  | 
|alignment_phased.bam	| alignment (input reads to the input gfa) if not provided by user | <img width="170" alt="image" src="https://github.com/katerinakazantseva/stRainy/assets/82141791/5f67678c-969d-4121-a694-dfddf694c0ce">
| multiplicity_stats.txt      | output statistics file (multiplicity and strain divergence info) | 
|phased_unitig_info_table.csv |  output statistics file (Length,Coverage, SNP rate) for phased unitigs | 
|reference_unitig_info_table.csv	|output statistics file (Length,Coverage, SNP rate) for reference unitigs | 
|preprocessing_data/ |	directory with files produced by stRainy (bam, fasta, etc) if not provided by user, also contain splitted gfa/fa/bam if split long unitigs enabled ||
|intermediate/10_fine_clusters.gfa | gfa with phased contigs only (same strain_unitigs.gfa) ||
|intermediate/20_extended_haplotypes.gfa |	gfa with phased and connected  contigs only (same strain_contigs.gfa)||
|intermediate/30_links_simplification.gfa | gfa after link simplification||
|intermediate/filtered_out_unitigs.txt | list of unitigs have not been processed (i.e. filtered by coverage) ||
|intermediate/bam/ | bam files with color tag based on cluster ID (availible only in debug mode) |<img width="170" alt="image" src="https://github.com/katerinakazantseva/stRainy/assets/82141791/565c08b4-da2b-4923-a4d2-5f3c8ed943df">
|intermediate/clusters/ |	csv files containing read names and corresponding cluster IDs |<img width="240" alt="image" src="https://github.com/katerinakazantseva/stRainy/assets/82141791/3665d84e-7486-405f-af79-de90874d5627">
|log_phase | log directory (phasing stage) ||
|log_transform	| log directory (transform stage) ||


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
