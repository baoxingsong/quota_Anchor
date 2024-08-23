# quota_Anchor
Here are the scripts and documents to conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave. We currently provide three visualization methods for syntenic results.
## Installation
You can simple by the following command get this software in a independent conda envirment. This is a beta version, so we haven't uploaded it to bioconda yet.
```
conda install xiaodli::quota_anchor
```
## Usage
### Help info
```
quota_Anchor -h
```
```
usage: quota_Anchor [-h] [-v] {pre_col,col,get_chr_length,dotplot,circle,line_2,line_proali} ...

Conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave:
options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

gene collinearity analysis:
  {pre_col,col,get_chr_length,dotplot,circle,line_2,line_proali}
    pre_col             Get longest protein file from gffread result and input file for gene collinearity analysis
    col                 Get gene collinearity result file
    get_chr_length      Get chromosome length and name info from fai file
    dotplot             Collinearity result visualization
    circle              Collinearity result visualization
    line_2              Collinearity result visualization
    line_proali         Anchors file from AnchorWave proali visualization
```
## Example
Here is an example to identify syntenic genes between maize and sorghum. The maize lineage has undergone a whole genome duplication (WGD) since its divergence with sorghum, but subsequent chromosomal fusions resulted in these species having the same chromosome number (n = 10). AnchorWave can allow up to two collinear paths for each sorghum anchor while one collinear path for each maize anchor.
### Make some folders rather than a folder may be more clearer
Working directory structure are as follows.
```
├── length_file
│   └── get_length.conf
├── raw_data
│   ├── Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3
│   ├── Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
│   ├── Zm-B73-REFERENCE-NAM-5.0.fa
│   └── Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
└── sb_zm
    └── config_file
        ├── circle.conf
        ├── collinearity.conf
        ├── dotplot.conf
        ├── line.conf
        └── pre_collinearity.conf
```
### Genome and annotation data preparation
Put the genome and gff file into the sb_zm/raw_data directory.
```
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gff3/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3.gz
gunzip *gz
```
### Modify config file and running pre_collinearity(maize vs sorghum)
This includes four steps(implemented in "quota_Anchor pre_col")
1. Extract and translate protein sequences from genome sequences and annotations
2. Identify and extract the longest protein sequence encoded by each gene
3. Protein sequence alignment using DIAMOND or conduct protein sequence alignment using BLASTp
4. Put the gene strand information and the blast result into a single file

Put the following information into the sb_zm/config_file/pre_collinearity.conf file
```
[gffread]
ref_genome_seq = ../raw_data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
ref_gff_file = ../raw_data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3
output_ref_pep_seq = sb.p.fa
query_genome_seq = ../raw_data/Zm-B73-REFERENCE-NAM-5.0.fa
query_gff_file = ../raw_data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
output_query_pep_seq = zm.p.fa
# The next line is the description of the S parameter of gffread(https://github.com/gpertea/gffread), you need to set True in general.
# -S    for -y option, use '*' instead of '.' as stop codon translation
use_S_parameter = True

[longest_pep]
out_ref_longest_pep_name = sorghum.protein.fa
out_query_longest_pep_name = maize.protein.fa

[align]
align=  diamond

[diamond]
# use ref protein seq construct database
database_name = sorghum.db
output_blast_result = sorghum.maize.diamond.blastp
max_target_seqs = 20
evalue = 1e-10

[blastp]
database_name = sorghum.blastp.db
dtype = prot
output_blast_result = sorghum.maize.blastp
evalue = 1e-10
max_target_seqs = 20
thread = 6
outfmt = 6

[combineBlastAndStrand]
out_file = sb_zm.table
bitscore = 100 
align_length = 0  
```
You can run this command in the sb_zm directory.
```
quota_Anchor pre_col -c ./config_file/pre_collinearity.conf
```
### Collinearity analysis(maize vs sorghum)
Put the following information into the sb_zm/config_file/collinearity.conf file and running colllinearity analysis.
```
[AnchorWave]
# The R value indicates the maximum number of occurrences of a gene in the collinearity file, and Q means the same as R.
# For maize and sorghum, maize has undergone an additional whole-genome duplication compared to sorghum.
# If sorghum is used as a reference, you can set R to 2 and Q to 1.
R = 2
Q = 1
maximum_gap_size = 25
delete_tandem = 0
tandem_dis = 5
input_file_name = sb_zm.table
output_coll_name = sb_zm.table.collinearity
```

```
quota_Anchor col -c ./config_file/collinearity.conf
```
### Visualzing by R code
The `sb_zm.table` could be visualized via the following R code:
```
library(ggplot2)
changetoM <- function ( position ){
  position=position/1000000;
  paste(position, "M", sep="")
}
data =read.table("sb.zm.table")
data$strand = data$V6==data$V12
data[which(data$strand),]$strand = "+"
data[which(data$strand==FALSE),]$strand = "-"

data = data[which(data$V8 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data = data[which(data$V2 %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")),]
data$V8 = factor(data$V8, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
data$V2 = factor(data$V2, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"))

plot = ggplot(data=data, aes(x=V10, y=V4))+geom_point(size=0.5, aes(color=strand))+facet_grid(V2~V8, scales="free", space="free" )+ theme_grey(base_size = 30) +
    labs(x="sorghum", y="maize")+scale_x_continuous(labels=changetoM) + scale_y_continuous(labels=changetoM) +
    theme(axis.line = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          axis.text.y = element_text( colour = "black"),
          legend.position='none',
          axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black") )
png("sorghum.maize.table.png" , width=2000, height=1500)
plot
dev.off()
```
<p align="center">
<img src="./quota_anchor/plots/sorghum.maize.table.png" width="800px" background-color="#ffffff" />
</p>

This file of `sb_zm.table.colinearity` could be visualized via the following R code:
```
library(ggplot2)
changetoM <- function ( position ){
  position=position/1000000;
  paste(position, "M", sep="")
}

data = read.table("sb_zm.table.collinearity", header=T)
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
<img src="./quota_anchor/plots/sorghum.maize.colinearity.png" width="800px" background-color="#ffffff" />
</p>

### Get chromosome length info
Put the following information into the sb_zm/length_file/get_length.conf file
```
# In the process of quotaAnchor pre_col, you can get fai file. 
# By fai file and raw GFF file , you can get length information.
# The maize length information example file are as follows.

#chr     length  total_gene
#chr1    308452471       5892
#chr2    243675191       4751
#chr3    238017767       4103
#chr4    250330460       4093
#chr5    226353449       4485
#chr6    181357234       3412
#chr7    185808916       3070
#chr8    182411202       3536
#chr9    163004744       2988
#chr10   152435371       2705
# select_fai_chr_startswith parameter, 
# number: software select chromosome name start with number.
# chr: software select chromosome name start with chr string.
# Chr: software select chromosome name start with Chr string.
[length]
fai_file = ../raw_data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.fai, ../raw_data/Zm-B73-REFERENCE-NAM-5.0.fa.fai
gff_file = ../raw_data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3, ../raw_data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
# By default, the first column of the lines starting with chr or Chr or CHR in the fai file are extracted for plotting.
select_fai_chr_startswith = number,CHR,chr,Chr:number,CHR,chr,Chr
length_file = sb_length.txt, zm_length.txt
```

### Visualzing by quota_Anchor
Put the following information into the sb_zm/length_file/circle.conf file
```
[circle]
collinearity = sb_zm.table.collinearity
ref_length = ../length_file/sb_length.txt
query_length = ../length_file/zm_length.txt
ref_prefix = sb-
query_prefix = zm-
font_size = 7
savefig = sb_zm.circle.png
```
```
quota_Anchor circle -c circle.conf
```
<p align="center">
<img src="./quota_anchor/plots/sb_zm.circle.png" width="800px" background-color="#ffffff" />
</p>

### Visualzing by quota_Anchor
Put the following information into the sb_zm/length_file/line.conf file
```
[line]
collinearity = sb_zm.table.collinearity
length_file = ../length_file/sb_length.txt, ../length_file/zm_length.txt
prefix = Sorghum, Maize
remove_chromosome_prefix = chr,CHR,Chr
text_font_size = 8
savefig = sb_zm.line.png
```
```
quota_Anchor line -c line.conf
```
<p align="center">
<img src="./quota_anchor/plots/sb_zm.line.png" width="800px" background-color="#ffffff" />
</p>
