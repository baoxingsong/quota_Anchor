# quota_Anchor&nbsp;&nbsp;[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/quota_anchor/README.html)[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<div align="center">

[**English**](./README.md) | [**中文简体**](./README_zh.md)

</div>

---
<details open>
 <summary><strong>Table of Contents</strong></summary>

<!-- TOC -->
- [quota\_Anchor  ](#quota_anchor)
  - [Installation](#installation)
  - [Usage](#usage)
    - [Help information](#help-information)
    - [Example of synteny analysis between maize and sorghum](#example-of-synteny-analysis-between-maize-and-sorghum)
      - [Preparation of genome and annotation file](#preparation-of-genome-and-annotation-file)
      - [Generate the longest protein sequence files](#generate-the-longest-protein-sequence-files)
      - [Generate the chromosome length files from fai and gff file](#generate-the-chromosome-length-files-from-fai-and-gff-file)
      - [Generate the table files that will be used as the input file for synteny analysis](#generate-the-table-files-that-will-be-used-as-the-input-file-for-synteny-analysis)
      - [Performing synteny analysis](#performing-synteny-analysis)
      - [Generate the longest cds sequence file](#generate-the-longest-cds-sequence-file)
      - [Calculate synonymous and non-synonymous substitution rates for syntenic pairs](#calculate-synonymous-and-non-synonymous-substitution-rates-for-syntenic-pairs)
    - [Homologous pairs and syntenic pairs visualization](#homologous-pairs-and-syntenic-pairs-visualization)
      - [Dotplot visualiztion](#dotplot-visualiztion)
      - [Circos visualiztion](#circos-visualiztion)
      - [Chromosome line style visualization](#chromosome-line-style-visualization)
    - [Maize gene/gene pairs classification](#maize-genegene-pairs-classification)
<!-- /TOC -->
</details>
Here are the documents to conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave.

## Installation

You can simply install the software via conda:

```command
conda create -n quota_Anchor bioconda::quota_anchor
```

## Usage

### Help information

```command
quota_Anchor -h
```

```text
usage: quota_Anchor [-h] [-v] {longest_pep,longest_cds,pre_col,col,get_chr_length,dotplot,circle,line,line_proali,ks,class_gene,kde,kf,trios,correct} ...

Conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Gene collinearity analysis:
  {longest_pep,longest_cds,pre_col,col,get_chr_length,dotplot,circle,line,line_proali,ks,class_gene,kde,kf,trios,correct}
    longest_pep         Generate the longest protein sequence file from gffread result.
    longest_cds         Generate the longest cds sequence file from gffread result.
    pre_col             Generate input file of synteny analysis(table file).
    col                 Generate gene collinearity result file.
    get_chr_length      Generate chromosome length and total number of genes information from fai file.
    dotplot             Generate collinearity dotplot or blast dotplot.
    circle              Collinearity result visualization(circos).
    line                Collinearity result visualization(line style).
    line_proali         Anchors file from AnchorWave proali to visualization(line style).
    ks                  Syntenic gene pair synonymous/non-synonymous substitution rate using yn00.
    class_gene          Genes/Pairs were classified as tandem, proximal, transposed, wgd/segmental, dispersed, and singleton.
    kde                 Focal species all syntenic pairs ks / block ks median histogram / gaussian kde curve.
    kf                  ks fitting plot.
    trios               generate trios by nwk tree.
    correct             Divergent ks peaks and correction.
```

### Example of synteny analysis between maize and sorghum

Here is an example to identify syntenic genes between maize and sorghum. The maize lineage has undergone a whole genome duplication (WGD) since its divergence with sorghum, but subsequent chromosomal fusions resulted in these species having the same chromosome number (n = 10). AnchorWave can allow up to two collinear paths for each sorghum anchor while one collinear path for each maize anchor.

#### Preparation of genome and annotation file

The current working directory contains genome files in fasta format and genome annotation files in gff format.

```bash
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gff3/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3.gz
gunzip *gz
```

For convenience, rename the file as follows:

```text
├── maize.fa
├── maize.gff3
├── sorghum.fa
└── sorghum.gff3
```

#### Generate the longest protein sequence files

The process, implemented in the function `quota_Anchor longest_pep`, consists of two main steps:

1. Protein sequences are extracted from genomes and annotation files based on genetic code tables.
2. For each gene, the longest protein sequence was identified and extracted to ensure the most complete characterization for further analysis.

```command
quota_Anchor longest_pep -f sorghum.fa,maize.fa -g sorghum.gff3,maize.gff3 -p sb.p.fa,zm.p.fa -l sorghum.protein.fa,maize.protein.fa -t 2 --overwrite -merge merged.pep.fa
```

#### Generate the chromosome length files from fai and gff file

The length files are required as input for generating table files, which are subsequently used for synteny analysis and plotting.

```command
quota_Anchor get_chr_length -f sorghum.fa.fai,maize.fa.fai -g sorghum.gff3,maize.gff3 -s 0-9:chr -o sorghum.length.txt,maize.length.txt --overwrite
```

#### Generate the table files that will be used as the input file for synteny analysis

1. Use DIAMOND/BLASTP/BLASTN for sequence alignment.
2. Combine the BLAST results and GFF file information into a single table file.

```command
quota_Anchor pre_col -a diamond -rs sorghum.protein.fa -qs maize.protein.fa -db sorghum.database.diamond -mts 20 -e 1e-10 -b sorghum.maize.diamond -rg sorghum.gff3 -qg maize.gff3 -o sb_zm.table -bs 100 -al 0 -rl sorghum.length.txt -ql maize.length.txt --overwrite
```

#### Performing synteny analysis

1. Generate collinearity result and specify `-r -q` parameter.

    ```command
    quota_Anchor col -i sb_zm.table -o sb_zm.collinearity -r 2 -q 1 -s 0 --overwrite 
    ```

2. Get `all` collinearity result and `remove` relative inversion gene pairs.

    ```command
    quota_Anchor col -i sb_zm.table -o sb_zm.collinearity -s 1 -a 1 --overwrite
    ```

3. Get `all` collinearity result and `retain` relative inversion gene pairs.

    ```command
    quota_Anchor col -i sb_zm.table -o sb_zm.collinearity -s 0 -a 1 --overwrite
    ```

#### Generate the longest cds sequence file

The process, implemented in the function `quota_Anchor longest_cds`, consists of two main steps:

1. Extract cds sequences from genome files and annotation files.
2. Identify and extract the longest CDS sequence for each gene.

```command
quota_Anchor longest_cds -f sorghum.fa,maize.fa -g sorghum.gff3,maize.gff3 -p sb.cds.fa,zm.cds.fa -l sorghum.cds.fa,maize.cds.fa -t 2 --overwrite -merge merged.cds.fa
```

#### Calculate synonymous and non-synonymous substitution rates for syntenic pairs

```command
quota_Anchor ks -i sb_zm.collinearity -a muscle -p merged.pep.fa -d merged.cds.fa  -o sb_zm.ks -t 16 --overwrite 
```

### Homologous pairs and syntenic pairs visualization

#### Dotplot visualiztion

1. Homologous gene pairs visualization using identity as a legend.

    ```command
    quota_Anchor dotplot -i sb_zm.table  -o sb_zm.table.identity.png -r sorghum.length.txt -q maize.length.txt -t order -r_label "Sorghum bicolor" -q_label "Zea mays" -w 1500 -e 1200 -use_identity --overwrite 
    ```

    <p align="center">
    <img src="./quota_anchor/plots/sb_zm.table.png" alt="sb_zm.identity.table" width="800px" background-color="#ffffff" />
    </p>

2. Syntenic pairs visualization using identity as a legend.

    ```command
    quota_Anchor dotplot -i sb_zm.collinearity  -o sb_zm.collinearity.identity.png -r sorghum.length.txt -q maize.length.txt -t order -r_label "Sorghum bicolor" -q_label "Zea mays" -w 1500 -e 1200 -use_identity --overwrite
    ```

    <p align="center">
    <img src="./quota_anchor/plots/sb_zm.collinearity.identity.png" alt="sb_zm.identity.collinearity" width="800px" background-color="#ffffff" />
    </p>

3. Syntenic pairs visualization using ks value as a legend.

    ```command
    quota_Anchor dotplot -i sb_zm.collinearity  -o sb_zm.collinearity.ks.png -r sorghum.length.txt -q maize.length.txt -t order -r_label "Sorghum bicolor" -q_label "Zea mays" -w 1500 -e 1200 -ks sb_zm.ks --overwrite
    ```

    <p align="center">
    <img src="./quota_anchor/plots/sb_zm.collinearity.ks.png" alt= sb_zm.ks.collinearity. png width="800px" background-color="#ffffff" />
    </p>

4. Syntenic pairs visualization by R code.

    This file of `sb_zm.collinearity` could be visualized via the following R code:

    ```R
    library(ggplot2)
    changetoM <- function ( position ){
      position=position/1000000;
      paste(position, "M", sep="")
    }

    data = read.table("sb_zm.collinearity", header=T)
    data = data[which(data$refChr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
    data = data[which(data$queryChr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")),]
    data$refChr = factor(data$refChr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
    data$queryChr = factor(data$queryChr, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"))

    plot = ggplot(data=data, aes(x=queryStart, y=referenceStart))+geom_point(size=0.5, aes(color=strand))+facet_grid(refChr~queryChr, scales="free", space="free" )+ 
      theme_grey(base_size = 30) +
      labs(x="maize", y="sorghum")+scale_x_continuous(labels=changetoM) + scale_y_continuous(labels=changetoM) +
      theme(axis.line = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
            axis.text.y = element_text( colour = "black"),
            legend.position='none',
            axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black") )

    png("sorghum.maize.colinearity.png" , width=2000, height=1500)
    plot
    dev.off()
    ```

    <p align="center">
    <img src="./quota_anchor/plots/sorghum.maize.collinearity.png" alt= sb_zm.dotplot.collinearity width="800px" background-color="#ffffff" />
    </p>

#### Circos visualiztion

Inter-species

```command
quota_Anchor circle -i sb_zm.collinearity -o sb_zm.circle.png -q maize.length.txt -r sorghum.length.txt -rn "Sorghum bicolor" -qn "Zea mays" -cf 9 -sf 9 -rm chr,Chr -fs 14,14 --overwrite
```

<p align="center">
<img src="./quota_anchor/plots/sb_zm.circle.png" alt= sb_zm.circle.collinearity width="800px" background-color="#ffffff" />
</p>

Intra-species

```command
quota_Anchor circle -i sb_sb.collinearity -o sb_sb.circle.png --overwrite -r ../sorghum.length.txt -q ../sorghum.length.txt -rn "sorghum" -qn "sorghum" 
```

<p align="center">
<img src="./quota_anchor/plots/sb_sb.circle.png" alt= sb_sb.circle.collinearity width="800px" background-color="#ffffff" />
</p>

#### Chromosome line style visualization

1. Visualization of syntenic pairs of two species
  
    ```command
    quota_Anchor line -i sb_zm.collinearity -o sb_zm.line.png -l sorghum.length.txt,maize.length.txt -n "Sorghum bicolor,Zea mays" --overwrite  
    ```

    <p align="center">
    <img src="./quota_anchor/plots/sb_zm.line.png" alt= sb_zm.line.collinearity width="800px" background-color="#ffffff" />
    </p>

2. Multi-species syntenic pairs visualization

    ```command
    quota_Anchor line -i os_sb.collinearity,sb_sv.collinearity -o os_sb_sv.line.png -l os.length.txt,sb.length.txt,sv.length.txt -n "Oryza sativa, Sorghum bicolor,Zea mays" -rm "chr,Chr" -cf 7 -sf 10 -fs 14,14 --overwrite
    ```

    <p align="center">
    <img src="./quota_anchor/plots/os_sb_sv.line.png" alt= os_sb_sv.line.collinearity width="800px" background-color="#ffffff" />
    </p>

### Maize gene/gene pairs classification

1. Synteny Analysis of Maize and Maize

    ```command
    quota_Anchor pre_col -a diamond -rs maize.protein.fa -qs maize.protein.fa -db maize.database.diamond -mts 20 -e 1e-10 -b maize.maize.diamond -rg maize.gff3 -qg maize.gff3 -o zm_zm.table -bs 100 -al 0 -rl maize.length.txt -ql maize.length.txt --overwrite
    quota_Anchor col -i zm_zm.table -o zm_zm.collinearity -s 0 -m 500 -W 5 -E -0.005 -D 25 -a 1 --overwrite
    ```

2. Download Musa balbisiana data and rename the filename

    ```bash
    wget https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_004837865.1/download\?include_annotation_type\=GENOME_FASTA\&include_annotation_type\=GENOME_GFF
    ```

    ```text
    ├── banana.B.fa
    └── banana.B.gff
    ```

3. Synteny Analysis of Banana.B and Maize

    ```command
    quota_Anchor longest_pep -f banana.B.fa -g banana.B.gff -p B.p.pep -l banana.B.pep -t 1 --overwrite
    quota_Anchor get_chr_length -f banana.B.fa.fai -g banana.B.gff -s CM01 -o banana.B.length.txt --overwrite
    quota_Anchor pre_col -a diamond -rs banana.B.pep -qs maize.protein.fa -db banana.B.database.diamond -mts 20 -e 1e-10 -b banana.B.maize.diamond -rg banana.B.gff -qg maize.gff3 -o bananaB_zm.table -bs 100 -al 0 -rl banana.B.length.txt -ql maize.length.txt --overwrite
    quota_Anchor col -i bananaB_zm.table -o bananaB_zm.collinearity -s 0 --overwrite -D 25 -a 1 
    ```

4. Classifying maize genes/gene pairs
    Unique mode

    ```command
    quota_Anchor class_gene -b maize.maize.diamond -g maize.gff3 -q zm_zm.collinearity -qr bananaB_zm.collinearity -o maize_classify_dir -p maize -s 1 -d 10 --overwrite -u
    ```

    <p align="center">
    <img src="./quota_anchor/plots/maize-unique.stats.gene.png" alt= maize-unique.stats.gene.png width="800px" background-color="#ffffff" />
    </p>

    <p align="center">
    <img src="./quota_anchor/plots/maize-unique.stats.pair.png" alt= maize-unique.stats.gene.png width="800px" background-color="#ffffff" />
    </p>

    Non-unique mode

    ```command
    quota_Anchor class_gene -b maize.maize.diamond -g maize.gff3 -q zm_zm.collinearity -qr bananaB_zm.collinearity -o maize_classify_dir -p maize -s 1 -d 10 --overwrite
    ```

    <p align="center">
    <img src="./quota_anchor/plots/maize.stats.gene.png" alt= maize-unique.stats.gene.png width="800px" background-color="#ffffff" />
    </p>

    <p align="center">
    <img src="./quota_anchor/plots/maize.stats.pair.png" alt= maize-unique.stats.gene.png width="800px" background-color="#ffffff" />
    </p>
