---
layout: default
title: Choosing a PopGenomics Approach
parent: Best Principles
has_children: true
#editor_options: 
#  chunk_output_type: console
---



# Choosing a PopGenomics Approach  
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

## Choosing a Population Genomics Approach  
Deciding on a sequencing approach for a new project can be challenging, and depends primarily on the questions you are interested in, but also constraints set by funding, existing genomic resources, sample size, genome characteristics, and expertise in the group. Below we list a number of great review articles and empirical comparison studies that can help when making the decision. This topic was also discussed during the MarineOmics seminars (link).  
Presently, our site covers four approaches (WGS=whole genome (re)sequencing):  

1) high coverage individual WGS (hcWGS)  
 i) description  
2) low coverage individual WGS (lcWGS)  
 i) description
3) pooled individuals (Poolseq)  
4) reduced-representation sequencing 
 i) can be split into restriction-based methods (eg RADseq) and targeted capture methods (eg exome capture)


For each approach, we provide guiding principles and tutorials for bioinformatic processing from raw data to genotypes. While we currently don't discuss library preparation methods in depth, we do include advice on how to organize your sequencing runs in order to improve quality control downstream.  

It should be noted that you don't need to choose only one approach for all samples! If you are starting a research program on a species without any genomic resources, you may want to start with de novo RADseq to understand basic population structure and genetic diversity, then once a reference genome is assembled move to high/low coverage WGS. In the [MarineOmics seminar](link), Misha Matz recommended choosing one individual to use for genome assembly, doing high coverage WGS on a few other individuals that represent most of the variation in your taxa, and then WGS many other individuals at low coverage. The high coverage individuals can be useful for imputing missing genotypes in the low coverage samples.  



References:  

* Matz 2018: great intro to ecological genomics methods in species with few to no prior genomic resources, accessible to those just getting started.  
* Fuentes-Pardo and Ruzzante 2017: a VERY thorough review of the different approaches to generating population genomics data, with non-model species in mind. Some of the content about relative cost and available software is a little dated, but otherwise a very valuable reference.   
* Benjelloun et al 2019: empirical comparison of high coverage WGS, low coverage WGS, random variants (eg. from RADseq), exome capture, and commercial SNP chip for mammals (genome size=2.6Gb). For neutral genomic diversity, found 5k-10k random variants were enough to accurately estimate, but comercial panales/exome capture had acertainment bias. To detect selection and accurate estimates of LD, at least 1M variants were required. For the studied species, 5x coverage was just as accurate as 12x coverage. Low coverage WGS (< 5x) should use genotype likelihood methods, however the rate of false positives for heterozygotes is higher than classic genotyping methods at higher coverage.  
* Dorant et al. 2019: Case study comparing Poolseq vs. reduced representation seq (GBS) vs Rapture (targeted capture of RAD loci) for inferring population structure with weak genetic differentiation. All three detected same weak pop structure pattern, but Poolseq gave Fst values 3-5x higher than Rapture/GBS.     
 


This table briefly summarizes the pros and cons of five sequencing approaches and their appropriateness for answering specific questions. If you are interested in answering multiple questions (eg, neutral population structure and adaptive variation) then you would often want a method that can do both (eg, WGS).  

|  | hcWGS | lcWGS | RADseq | Poolseq | Targeted capture (eg, exome capture) |
|---|---|---|---|---|---|
| Pros | provides high quality, high density genotypes; can be used to improve development of a reference genome | provides high density genotypes at a reduced cost | generally cheaper per sample, density of genotypes can be tuned to fit question, does not require a reference | many individuals can be mixed for a low library prep cost, only option for some larval studies | allows consistent capture of loci between batches, good for sequencing targeted regions across many samples  |
| Cons | most expensive per sample, requires a reference genome, data requires large computational resources | requires a reference, can produce false positive heterozygotes, sensitive to batch effects, may be inappropriate to use individual SNP calls for some analyses | only 1-5% of genome covered which limits studying adaptive variation, de novo assembly can result in paralogs if not tuned or filtered well, can be hard to capture the same loci between batches | need higher sequencing depth so it can be more expensive than RAD (but less than other options), no info on individual genotypes (duh), requires a reference, requires some replicates to account for batch effects | expensive to start a new project and requires reference of some kind to design baits (unless using method like eecSeq), only gives info on targeted region so can result in ascertainment bias for some analyses |
| Population structure, genetic diversity summary stats | great, but if only interested in population structure then this may be $$$$ overkill | good, but may be limited to popgen analyses based on genotype likelihoods | great, esp. for many samples | good for methods that only use allele frequencies, won't work for individual based methods (eg Admixture) | ascertainment bias is likely to skew results, esp. for measures of genetic diversity |
| Demography (migration rates, population size through time) | great, may be overkill unless you incorporate haplotype-based methods | good, esp. with some high coverage samples to help with imputation | good for methods based on site frequency spectrum (eg moments), but not ideal for methods using extended haplotype or phasing info | active area of method development, but still not common except in cases of multiple temporal samples (cite) | likely not appropriate due to ascertainment bias |
| Signatures of selection, specific genomic regions of interest (GWAS) | great | good | Hotly debated, but unless your genome is way too big, WGS is the way to go | good if covering majority of the genome | only works for targeted regions |
| Phylogenetic inference | computationally challenging, but good at all divergence levels | good for tree shapes but not branch lengths, good for extracting organelle sequences and developing cost-effective tools (eg primers) for species and hybrid ID studies | good at shallow to medium divergence levels, but need to play with filtering parameters | not possible | an improvement over RAD as there will be less missing data among samples |
| Genetic crosses or mapping panels | likely $$$$ overkill | great if you have a genome | great as it allows many individuals | maybe | poor |
| molecular evolution (where accurate low-frequency alleles are required) | best | poor | poor | poor | good if only interested in targeted regions |
| Evolution of protein coding genes | great but requires decent annotations on genome | ? | poor, as you will miss many genes | ? | will miss introns and regulatory regions, but otherwise good |



