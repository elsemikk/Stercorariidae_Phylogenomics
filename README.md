# Stercorarius 
## Phylogenomics of Stercorariidae

## About this project
This project examines the phylogenomics of the Stercorariidae, an enigmatic family of seven seabirds, the skuas and jaegers. The goals of this project are to understand the role of introgression in the evolution of the Stercorariids, and evaluate the hypothesis that the Pomarine Jaeger (*S. pomarinus*) or the Great Skua (*S. skua skua*) are species of hybrid origin.

## About this repository
This repository contains all of the code and instructions needed to replicate the analyses of this project. We have organized the repository into folders corresponding to different steps of the project from data cleaning to final results. The entire pipeline could be rerun starting from the raw sequencing reads that will be available online from the NCBI SRA; alternatively, we provide intermediate files and output for each step where file sizes allow.

To rerun analyses on a different machine, it will be necessary to change the paths to the working directory (which is this repository, a folder named Stercorarius) as well as any paths to installed programs. The folders have been zipped in order to upload them to Dryad, and larger files have been compressed with gzip to reduce the file size of the repository. If the analyses are to be replicated starting from these intermediate files, they would need to be uncompressed, which can be done with the "gunzip" command.  

## Dataset  
We sequenced the whole genomes of each of the seven currently-recognized species of Stercorariidae.  

Throughout this repository, the samples are referred to by their short sample code names. Here is the correspondance between the code names and their species ID and sample names as recorded in Table S1 of the Supplementary Materials:     
LTJA_MKP990: *S longicaudus* (MKP990 (ROM Birds 156112))  
PAJA_B20730: *S parasiticus* (B-20730 (BROM573-07))  
PAJA_USNM606730: *S parasiticus* (Genbank SAMN12253778 (B10K-DU-001-20))  
POJA_MKP1559: *S pomarinus* (MKP1559 (ROM Birds 157062))  
POJA_4: *S pomarinus* (POJA4) (This specimen is pending an accession number at the Royal Ontario Museum, so I refer to it as POJA4)  
POJA_IB2659: *S pomarinus* (1B-2659 (ROM Birds 159823))  
GRSK_MKP1592: *S skua* (MKP1592)  
GRSK_MKP1593: *S skua* (MKP1593)  
CISK2: *S antarcticus lonnbergi* (C54 (Chatham Island Skua 2))  
CISK55: *S antarcticus lonnbergi* (C72 (Chatham Island Skua 5))  
CISK_3: *S antarcticus lonnbergi* (C55 (Chatham Island Skua 3))  
CHSK_MKP2451: *S chilensis* (MKP2451  (ROM Birds 158363))  
ANSK7: *S maccormicki* (E67 (Antarctic Skua 7))  
ANSK01: *S maccormicki* (E23 (Antarctic Skua 10))  
ANSK8: *S maccormicki* (E68 (Antarctic Skua 8))  
Alca_torda: *Alca torda* (Alcidae, outgroup)  
Uria_lomvia: *Uria lomvia* (Alcidae, outgroup)  

## Software  
These are the version numbers of software I have installed:  
* RagTag v2.0.0  
* sratoolkit prefetch v2.9.6-1-ubuntu64  
* FastQC v0.11.8  
* MultiQC v1.7  
* Trimmomatic v0.39  
* fastp v0.23.2  
* GNU parallel v20161222  
* NOVOPlasty v4.2.1  
* bowtie2 v2.4.4-linux-x86_64  
* samtools v1.14  
* qualimap v2.2.1  
* Rscript R scripting front-end version 3.6.3 (2020-02-29)  
* *R* libraries:  
    * ggplot2 v3.3.5  
    * dplyr v1.0.7  
* Mesquite v3.40 build 877  
* AliView v1.26  
* mafft v7.490  
* IQ-TREE v2.1.3  
* NCBI blast+ v2.12.0+  
* bcftools v1.14  
* GenomeTools v1.5.10  
* VCFtools v0.1.17  
* PLINK v1.90b6.24 64-bit (6 Jun 2021)   
* Tracer_v1.7.1  
* BEAST v2.6.7 (with StarBeast3 v1.0.4)  
* bedtools v2.30.0  
* Emboss v6.6.0  
* shapeit v2.904.3.10.0-693.11.6.el7.x86_64  
* seqkit v0.11.0  
* angsd v0.930 (htslib: 1.9) build(May 28 2021 01:24:37)  
* Rohan (no version number, accessed June 6 2020)  
* Dsuite v0.4  


# Contents
This repository is organized into three main sections: 1) Dataset, 2) Phylogeny, and 3) Introgression. Each of these sections contains multiple numbered steps for each analysis. Each step has a markdown (`.md`) file containing all the code and instructions detailing how each step was run. It also contains a folder by the same name, containing various intermediate files and raw results from that step. A few analyses were done in R; these analyses have an `.Rmd` file with the same name as their `.md` file containing all the R code and output from that step. The `.Rmd` files have an `html` file by the same name, containing output from the R code. This project used two different datasets: mapped to *Alca torda* or to *Stercorarius parasiticus*. Steps that use the *Stercorarius* reference genome datasets have the suffic "c", and steps that use the *Alca* reference genome dataset have the suffix "d".  

The tree files for each of the main phylogenies in the paper are also together in the folder `Phylogenetic_tree_files.zip` to be more easily accessed.  

Here is an overview of how the different markdown files/folders map to each analysis in the paper:  

## Phylogenies:  

### mtDNA ML phylogeny  
1.0_datasets -> 1.1c_read_trimming -> 1.2_mitogenome_assembly -> 2.1.1_mtDNA_ML_phylogeny  
**Final result: Stercorariidae_mtDNA.tre**  

### mtDNA BEAST2 Bayesian phylogeny  
1.0_datasets -> 1.1c_read_trimming -> 1.2_mitogenome_assembly -> 2.1.2c_mtDNA_BEAST_phylogeny   
**Final result: Stercorariidae_MitogenomesCorrectedSinglecopy.4_runs.tree**  

### MC1R Gene tree  
1.0_datasets -> 1.1c_read_trimming -> 1.3c_mapping -> 1.5c_genotyping -> 2.2c_MC1R_GeneTree  
**Final result: MC1R.POJA4_haplotypes_corrected.fasta.rooted.treefile**  

### ASTRAL species tree  
1.0_datasets -> 1.1c_read_trimming -> 1.3c_mapping -> 1.6.2c_phylogenetic_blocks -> 2.5c_phylogeny_ASTRAL  
**Final result (autosomal): species.Astralblocks.unlinked.partitioned.tre**  
**Final result (chromosome Z): species.Astralblocks.unlinked.partitioned.chrZ.nochrW.tre**  

### Stabeast3 species tree  
1.0_datasets -> 1.1c_read_trimming -> 1.3c_mapping -> 1.5c_genotyping -> 1.6.1c_phylogenetic_blocks -> 2.6c_Starbeast3_phylogeny  
**Final result: species.16runs.100MILL.tree**  

### Phylonet phylogenetic network  
1.0_datasets -> 1.1c_read_trimming -> 1.3c_mapping -> 1.5c_genotyping -> 1.6.1c_phylogenetic_blocks -> 2.8c_phylonet  
**Final result: Phylonet_best_network.tree**  

## Introgression analyses:  

### ABBA BABA test & Fbranch  
1.0_datasets -> 1.1c_read_trimming -> 1.3d_mapping -> 1.5d_genotyping (Stercorarius.RAZObowtie.merged.filteredSNPs.auto.stringent.max17.vcf.gz) -> 3.1d_Dsuite_ABBABABA  

### Genome scan of *f*~~dM~~  
1.0_datasets -> 1.1c_read_trimming -> 1.3d_mapping -> 1.5d_genotyping -> 3.1d_Dsuite_ABBABABA  

### *D*~~3~~, genome scan of genetic distances, and average divergence between samples  
1.0_datasets -> 1.1c_read_trimming -> 1.3c_mapping -> 1.6.2c_phylogenetic_blocks -> 3.2.1c_D3_distancematrix  

### mtDNA vs nuclear distances  
1.0_datasets -> 1.1c_read_trimming -> 1.3d_mapping -> 1.6.2d_phylogenetic_blocks -> 3.2.2d_mtDNA_vs_nuclear_distances  

### Mitonuclear introgression analysis  
1.0_datasets -> 1.1c_read_trimming -> 1.3c_mapping -> 3.5_mitonuclear  

### W chromosome homolog introgression analysis  
1.0_datasets -> 1.1c_read_trimming -> 1.3c_mapping -> 3.6_W_Z_homologs  
