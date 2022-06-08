[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

# metaPhase

Graph-based phasing algorithm, that takes a de novo assembly graph (in gfa format) and simplifies it by combining phasing information and graph structure.

## Installation
### Requirements
```
python3
pip
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
python3 ./main.py
```

## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
