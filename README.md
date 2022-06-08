[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

# metaPhase

Graph-based phasing algorithm, that takes a de novo assembly graph (in gfa format) and simplifies it by combining phasing information and graph structure.

**metaPhase.phase** - perform reads clustering according to SNP positions using community detection approach

**metaPhase.transfom** - transform assembly graph 

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

## Parameters

Please specify parameters in `params.py`:
- `bam` - full path to alignment file in bam format
- `gfa` - full path to assembly gfa (produced with [**metaFlye**](https://github.com/fenderglass/Flye))
- `gfa_transformed` - full path to store transformed assembly gfa
- `gaf_file` - full path to alignment file in gaf format (produced with [**GraphAligner**](https://github.com/maickrau/GraphAligner))
- `snp` - full path to alignment file in vcf format (Optional). 

It is not recommended to change parameters below:
- `R` - reads max mismatch count (default=1)
- `I` - reads min intersection  (default=100)
- `AF` - allele frequency (default=0.1)
- `clipp` - filter alignment. max clipping lenght (default=100)
- `min_mapping_quality` - filter alignment. min mapping quality (default=20)
- `min_base_quality` - filter alignment. min base quality (default=0)
- `min_al_len` - filter alignment. min alignment lenght (default=1000)
- `de_max` - filter alignment. max identity (default=0.05)

## Run

```
python3 ./metaPhase.phase.py
python3 ./metaPhase.transform.py
```
## Credits

metaPhase was developed under the supervision of Mikhail Kolmogorov https://github.com/fenderglass

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg



