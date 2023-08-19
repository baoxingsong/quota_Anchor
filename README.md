# quota_Anchor
Here are the scripts and documents to conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave.
## Dependencies
AnchorWave \
[gffread](https://github.com/gpertea/gffread) \
python3 \
[NumPy](https://numpy.org/) \
BLAST or [DIAMOND](https://github.com/bbuchfink/diamond)
## Example
### Genome and annotation data preparation
Here is an example to identify syntenic genes between maize and sorghum. The maize lineage has undergone a whole genome duplication (WGD) since its divergence with sorghum, but subsequent chromosomal fusions resulted in these species having the same chromosome number (n = 10). AnchorWave can allow up to two collinear paths for each sorghum anchor while one collinear path for each maize anchor.

```
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gff3/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3.gz
gunzip *gz
```

### Extract and translate protein sequences from genome sequences and annotations
```
gffread -g Zm-B73-REFERENCE-NAM-5.0.fa -y zea.p.fa Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
gffread -g Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -y sb.p.fa Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3
```

### Identify and extract the longest protein sequence encoded by each gene
```
python3 longestPeps.py -g Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -f Zm-B73-REFERENCE-NAM-5.0.fa -p zea.p.fa -o maize.protein.fa
python3 longestPeps.py -g Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -f Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -p sb.p.fa -o sorghum.protein.fa
```


### Protein sequence alignment
Conduct protein sequence alignment using DIAMOND
```
sed -i -e 's/\./*/g' maize.protein.fa
sed -i -e 's/\./*/g' sorghum.protein.fa
diamond makedb --in maize.protein.fa --db maize
diamond blastp --db maize -q sorghum.protein.fa -k 5 -e 1e-10 -o sorghum.maize.blastp

```
Or conduct protein sequence alignment using BLASTp
```
./blast-2.2.26/bin/formatdb -p T -i maize.protein.fa -n maize.protein
./blast-2.2.26/bin/blastall -i sorghum.protein.fa -d maize.protein -p blastp -e 1e-10 -b 5 -v 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore " -o sorghum.maize.blastp
```

### Put the gene strand information and the blast result into a single file
```
python3 combineBlastAndStrandInformation.py -r Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -q Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -b sorghum.maize.blastp -o sorghum.maize.table
```

This table could be visualized via the following R code:
```
library(ggplot2)
changetoM <- function ( position ){
  position=position/1000000;
  paste(position, "M", sep="")
}
data =read.table("sorghum.maize.table")
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
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          axis.text.y = element_text( colour = "black"),
          legend.position='none',
          axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black") )
png("sorghum.maize.table.png" , width=2000, height=1500)
plot
dev.off()
```
<p align="center">
<img src="./plots/sorghum.maize.table.png" width="800px" background-color="#ffffff" />
</p>

### Conduct strand and WGD aware syntenic gene identification using AnchorWave

```
anchorwave pro -i sorghum.maize.table -n sorghum.maize.colinearity -R 1 -Q 2
```

This file of `sorghum.maize.colinearity` could be visualized via the following R code:
```
library(ggplot2)
changetoM <- function ( position ){
  position=position/1000000;
  paste(position, "M", sep="")
}

data = read.table("sorghum.maize.colinearity", header=T)
data = data[which(data$queryChr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data = data[which(data$refChr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")),]
data$queryChr = factor(data$queryChr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
data$refChr = factor(data$refChr, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"))

plot = ggplot(data=data, aes(x=queryStart, y=referenceStart))+geom_point(size=0.5, aes(color=strand))+facet_grid(refChr~queryChr, scales="free", space="free" )+ 
	theme_grey(base_size = 30) +
    labs(x="sorghum", y="maize")+scale_x_continuous(labels=changetoM) + scale_y_continuous(labels=changetoM) +
    theme(axis.line = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          axis.text.y = element_text( colour = "black"),
          legend.position='none',
          axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black") )

png("sorghum.maize.colinearity.png" , width=2000, height=1500)
plot
dev.off()
```
<p align="center">
<img src="./plots/sorghum.maize.colinearity.png" width="800px" background-color="#ffffff" />
</p>
