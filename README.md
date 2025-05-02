# Heteroplasmy in Human Mitochondrial Genome Mutation Rate Differences

This project investigates the discrepancy observed between mitochondrial genome mutation rates calculated from pedigree-based studies and substitution rates from phylogeny-based analyses. The primary aim was to assess hypotheses explaining this difference, specifically focusing on the role of heteroplasmy—where cells contain multiple mitochondrial DNA (mtDNA) genotypes—in inflating mutation rate estimates.

*This work was conducted to fulfill MDSC 508 and the Bachelors of Health Sciences in Bioinformatics Honours thesis at the University of Calgary.*

## Objectives

- Evaluate the reasons behind the divergent mutation rates observed in pedigree and phylogenetic studies.

- Investigate the impact of heteroplasmy on mutation rate estimates by examining mother-child mtDNA pairs.

- Identify and analyze mutation rate gradients within the mitochondrial genome.

## Methodology

### Data Collection and Curation: 
Collected and harmonized mitochondrial DNA variant data from four publicly available studies involving a total of 3,598 individuals. Data included mother-child pairs, variant positions, major and minor alleles, and minor allele frequencies (MAF).

### Mutation Rate Calculations: 
Developed and implemented R scripts provided in this repository to calculate mutation rates from mother-child paired mitochondrial DNA. Mutation rates were separately computed for blood and cheek tissue samples and combined datasets.

### Heteroplasmy Filtering: 
Utilized custom scripts to apply filters based on MAF to determine the impact of heteroplasmy on mutation rate estimates. Rates were recalculated at varying heteroplasmy thresholds.

### Mutation Gradient Analysis: 
Created and analyzed histograms (included in this repository) to identify spatial patterns and hotspots in mutation frequencies across the mitochondrial genome.

##Repository Contents

### R Scripts: 
Fully documented scripts for data processing, mutation rate calculation, heteroplasmy filtering, and histogram plotting.

### Figures: 
Generated histograms illustrating mutation hotspots and changes in mutation rate estimates due to heteroplasmy filtering.

### Datasets: 
Curated datasets used for mutation rate calculation and heteroplasmy analysis, clearly labelled and formatted.

## Key Findings

- Initial mutation rate estimates significantly exceeded literature values, showing a 14.9-to-22.4-fold increase when including heteroplasmic variants.

- Implementing a heteroplasmy filter reduced these estimates to align closely with higher-end literature rates (2.6 mutations/site/million years).

- Observed distinct mutation hotspots corresponding to known hypervariable regions and identified an unreported peak in cheek tissue around nucleotide position 11,000.

## Conclusion

The study highlights the importance of accounting for heteroplasmy in mitochondrial evolutionary studies and suggests that previously proposed evolutionary models need adjustments to accommodate heteroplasmic variance. Additionally, the findings indicate that tissue-specific mutation processes significantly influence mutation rate estimates.

## Future Directions

Further research should involve modelling mitochondrial inheritance at the cellular level, explicitly accounting for heteroplasmy, to improve accuracy in evolutionary rate estimations and deepen the understanding of mitochondrial genome dynamics.
