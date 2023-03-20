# Structure analysis for Narnavirus

This repository contains structure analysis for Narnavirus 20S/23S and related viruses.



## Dependencies

The analysis is done with ```Python3.9```, ```mfold-3.6``` (Linux distribution)and ```ViennaRNA-2.5.0``` (Linux distribution). 

python environment

```bash
conda create --name Narna_structure --file requirements.txt
conda activate Narna_structure
```

Mfold-3.6 needs to be installed locally

```bash
# run configuration for mfold-3.6
cd tools/mfold-3.6
chmod 755 tools/mfold-3.6/configure
./configure

# compile and install
make
make install

# if libgfortran.so.4 is not found, install libgfortran4
#apt-get install libgfortran4
```

ViennaRNA-2.5.0 doesn't need to be installed. The part of program used in the analysis was contained in ```tools/ViennaRNA-2.5.0```



## Dataset

Data used for analysis are contained in data/

**data/Narna_20S_refGenome.fasta**: reference genome for Saccharomyces 20S RNA narnavirus. Obtained from NCBI https://www.ncbi.nlm.nih.gov/nuccore/NC_004051.1

**data/Narna_23S_refGenome.fasta**: reference genome for Saccharomyces 23S RNA narnavirus. Obtained from NCBI https://www.ncbi.nlm.nih.gov/nuccore/NC_004050.1

**data/Narnaviridae_refGenome.fasta**: reference genomes for Narnaviridae.

**data/Mitoviridae_refGenome.fasta**: reference genomes for Mitoviridae.

The viruses reference genome and metadata are downloaded from NCBI virus database (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Genome&VirusLineage_ss=Viruses,%20taxid:10239). To extract references for Narnaviridae and Mitoviridae.

```bash
# Narnaviridae
python get_seqList.py --filename ../data/viralGenomeGenbank.txt --family Narnaviridae --ref > ../data/Narnaviridae_list.txt
python get_seq.py --filename ../data/viralGenomeGenbank.fasta --seqlist ../data/Narnaviridae_list.txt > ../data/Narnaviridae_refGenome.fasta

# Mitoviridae
python get_seqList.py --filename ../data/viralGenomeGenbank.txt --family Mitoviridae --ref > ../data/Mitoviridae_list.txt
python get_seq.py --filename ../data/viralGenomeGenbank.fasta --seqlist ../data/Mitoviridae_list.txt > ../data/Mitoviridae_refGenome.fasta
```



## Structure segment

The structure segment section aims to segment the sequence into regions with independent secondary structures.

- Predict RNA secondary structure by Mfold-3.6 (folding temperature T = 30&deg;C)
- Average between the ensemble of optimal and suboptimal structures
- Ignore long range interactions beyond 200 nucleotides apart
- Identify independent structural regions

```bash
python segment_structure.py --inputFasta ../data/Narna_23S_refGenome.fasta --mfoldOut --outdir ../output/
```

Outputs:

- **output/Narna_23S_refGenome_dotplot.png**: Dotplot of RNA structure, with stems marked by blue rectangles
- **output/Narna_23S_refGenome_regions.png**: Structural regions marked on the dotplot as red rectangle
- **output/Narna_23S_refGenome/mfold_***: Mfold output files (if ```--mfoldOut``` is given to keep the mfold outputs)

 

## Structural comparison with permutated sequences

The permutated sequences all shows larger folding energy than reference genome, indicating that the wild-type virus has envolved into a stable structure.

```bash
python compare_permutate.py --inputFasta ../data/Narna_23S_refGenome.fasta --outdir ../output/
Rscript compare_permutate_plot.R ../output/Narna_23S_refGenome.dat Narna23S
```

Outputs:

- **output/Narna_23S_refGenome.dat**: Recoding of folding energy for wildtype and permutated sequences. 
- **output/Narna_23S_energy_separate.pdf**: Boxplot of folding energy



## MFE comparison with viruses from similar family

The trend that wild-type viruses have lower folding energy than permutated sequence is not only restricted to 23S, but also viruses from similar family (Narnaviridae and Mitoviridae), regardless of GC contents.

```bash
python compare_gcContent.py
Rscript compare_gcContent_plot.R ../output/Narna_energy.txt ../output/Mito_energy.txt
```

Outputs:

- **output/Narna_energy.txt, output/Mito_energy.txt**: Folding energy statistics for genomes in corresponding viridae. Columns corresponds to:
  - Viridae
  - Reference genome name (accession | description)
  - Genome length
  - GC content
  - Mean of wild-type folding energy (optimal and suboptimal)
  - Standard error of wild-type folding energy
  - Numbers of structure considered (wild-type)
  - Minimum folding energy (wild-type)
  - Mean of permutated folding energy
  - Standard error of permutated folding energy
  - Numbers of structure considered (permutation)
  - Minimum folding energy (permutation)

- **output/energy2gc.pdf**: Showing the relationship between folding energy and GC contents



## Contact map

This section generates a contact map of a given sequence for structure visualization. RNAfold is used for data prediction because it could generate the base pairing probability between 2 bases.

```bash
python plot_contactMap.py --inputFasta ../data/Narna_23S_refGenome.fasta --rnafoldOut --outdir ../output/
```

Outputs:

- **output/Narna_23S_refGenome.png**: Contact map, with color represents paring probability (0-1)
- **output/Narna_23S_refGenome/rnafold_***: RNAfold output files (if ```--rnafoldOut``` is given to keep the mfold outputs)



