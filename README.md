# Asterophryinae_phylogenetic_data

Data for Asterophryinae frog phylogenetic data: nuclear and mitochondrial sequence data, site and locality, lifestyle, and frog communities

This repository contains the following data types:
- GenBank accession numbers for 3 nuclear and 2 mitochondrial loci for 236 frog samples
- Alignments and input files for BEAST2 and IQTREE phylogenetic analyses, with geological timings.
- Primer sequences 
- Metadata on anuran lifestyle, GPS location, site, locality, country
- Ancestral reconstructions of lifestyle
- Multispecies communities
- R Code to reproduce all figures and analyses

*This is the data respository for the manuscript:*

Ethan C. Hill, Mary J. Jarman, Claire J. Fraser, Diana F. Gao, Elizabeth R. Henry, Allison R. Fisher, Bulisa Iova, Allen Allison, and Marguerite A. Butler (in Review) A large nuclear and mitochondrial sequence dataset for phylogenetic analysis of the hyperdiverse Asterophryinae frogs of the New Guinea region with data on lifestyle, GPS coordinates, elevation, and multispecies communities with accompanying code. Submitted to Data in Brief 

*And related research article:*

Ethan C. Hill, Claire J. Fraser, Diana F. Gao, Mary J. Jarman, Elizabeth R. Henry, Bulisa Iova, Allen Allison, Marguerite A. Butler. 2022. Resolving the deep phylogeny: Implications for early adaptive radiation, cryptic, and present-day ecological diversity of Papuan microhylid frogs. Molecular Phylogenetics and Evolution. 177:107618. ISSN 1055-7903. https://doi.org/10.1016/j.ympev.2022.107618

### The repository contains: 
#### Tables
* Table 1: **Specimen metadata** including collection site, type locality, species citation, lifestyle based on microhabitat use and citation, gps, elevation, site name and number, geological terrane, GenBank accession numbers by locus.
* Table 2: **Primer sequences** redesigned for this study, based on Rivera et al. [15]. Product length reported in base pairs, and Tm in °C. Starting touchdown PCR temperature range in °C is reported in TD range. Outer indicates the outermost most forward/reverse primer for a locus.
* Table 3: **Data partitions and their best-fit evolutionary models** for sequence evolution using PartitionFinder2 allowing models to vary by locus and codon position (15 possible partitions).

#### Figures
* Figure 1: Time calibrated Bayesian inference phylogeny for 233 samples of Asterophryinae generated using BEAST2. Outgroups have been removed.
* Figure 2: Time calibrated maximum likelihood phylogeny generated for 233 samples of Asterophryinae using IQTREE. Outgroups have been removed.
* Figure 3: Nuclear-DNA-only time calibrated Bayesian inference phylogeny for all 236 samples of Asterophryinae generated with BEAST2.
* Figure 4: Mitochondrial-DNA-only time calibrated Bayesian inference phylogeny for all 236 samples of Asterophryinae pruned to 218 taxa generated with BEAST2.

#### BEAST 2 Files:
* alignment_asterophryinae_07192022.nex: aligned dataset of all 5 loci
* beast_infile_asterophryinae_07192022.xml: BEAST input file for all 5 loci
* beast_tree_asterophryinae.nex: BEAST consensus tree from all 5 loci
* alignment_asterophryinae_nuclear_07192022.nex: nuclear-only dataset alignment.
* beast_infile_asterophryinae_nuclear_07192022.xml: nuclear-only BEAST input file.
* beast_tree_asterophryinae_nuclear.nex: nuclear-only BEAST consensus tree.
* Alignment_asterophryinae_mitochondrial_07192022.nex: mitochondrial-only beast input file.
* beast_infile_asterophryinae_mitochondrial_07192022.xml: mitochondrial-only BEAST input file.
* beast_tree_asterophryinae_mitochondrial.nex: mitochondrial-only BEAST consensus tree.
* beast_218_tree_asterophryinae.nex: consensus BEAST phylogeny pruned to one species per site.

#### IQTREE Files:
* asterophryinae_07192022.phy: aligned dataset of all 5 loci for IQTREE analysis.
* asterophryinae_dates.txt: time calibrations for the phylogeny based on geological events.
* asterophryinae_partitions.nex: partition file for the IQTREE analysis.
* iqtree_timetree_asterophryinae.contree: consensus tree generated from the IQTREE analysis.
* iqtree_218_timetree_asterophryinae.contree: consensus IQTREE phylogeny pruned to one species per site.

#### R scripts which produce analyses and figures:
* asterophryinae_phylogeny_analysis_figures.R: script used for all analyses, to plot phylogenies, conduct ancestral state reconstruction, and calculate distances between collection sites. 
* clean_functions.R: accessory functions used in “asterophryinae_phylogeny_analysis_figures.R”
* gencolorABC.csv: color code specifications for genera used in "asterophryinae_phylogeny_analysis_figures.R"
