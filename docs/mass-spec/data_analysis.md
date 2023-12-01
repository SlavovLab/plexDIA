---
layout: default
title: Data Analysis
nav_order: 3
permalink: mass-spec/plexDIA_analysis
description: "plexDIA data searching and analysis: Tutorial on searching and analyzing multiplexed DIA data from plexDIA"
nav_exclude: false
---
{% include social-media-links.html %}

# Analyzing plexDIA data
{:.no_toc}

&nbsp;

{: .fs-5 .fw-300}
This section organizes libraries, tutorials and links to software and pipelines for searching, optimizing and analyzing plexDIA data.

* Will be replaced with the ToC, excluding the section header
{:toc}

&nbsp;

## Searching plexDIA data


### Libraries for searching plexDIA data

These libraries are for searching spectra from mTRAQ-labeled peptides. All libraries and FASTAs below include canonical and isoform protein sequences from Swissprot.

### Predicted libraries output from DIA-NN:

1.  [Human, Yeast, E. coli](https://drive.google.com/file/d/1k6PaBpth40Tci2snG8sWFG645Nub9iQw/view?usp=drive_link) & corresponding [Human, Yeast, E. coli FASTA](https://drive.google.com/file/d/1bFWZ2lptAYuQByCcNfhBG163_CE-iVQu/view?usp=drive_link) (source: [Derks et al, 2022](https://www.nature.com/articles/s41587-022-01389-w))

2.  [Human](https://drive.google.com/file/d/1srNY0Nz8b-oRISFf3XFxUI-XncmDjOFZ/view?usp=drive_link) & corresponding [Human FASTA](https://drive.google.com/file/d/1gBFWDbTQJCrWkK5rMUDxZhDpfsglWxVl/view?usp=drive_link) (source: [Derks et al, 2022](https://www.nature.com/articles/s41587-022-01389-w))


### Empirical, curated libraries:

1.  [Human, 100 cells of PDAC, melanoma, and monocytes from Q-Exactive](https://drive.google.com/file/d/1XPrTLq1WxXg7lfI3No1S9frOvE51V1Sx/view?usp=drive_link) & corresponding [Human FASTA](https://drive.google.com/file/d/1gBFWDbTQJCrWkK5rMUDxZhDpfsglWxVl/view?usp=drive_link) (source: [Derks et al, 2022](https://www.nature.com/articles/s41587-022-01389-w))
2.  [Human, THP-1 macrophages from timsTOF SCP](https://drive.google.com/file/d/1ldCjhKOhRpPfrEc7GQNaHztP_nwlrj1g/view?usp=drive_link) & corresponding [Human FASTA](https://drive.google.com/file/d/1gBFWDbTQJCrWkK5rMUDxZhDpfsglWxVl/view?usp=drive_link) (source: manuscript in preparation Derks et al, 2023)



## Tutorial on searching plexDIA data
Searching [plexDIA data](https://scp.slavovlab.net/Derks_et_al_2022) with [DIA-NN](https://github.com/vdemichev/DiaNN/releases/tag/1.8.1) is described in this [tutorial](https://youtu.be/0Wmg9LjDtgE). The tutorial is for DIA-NN v1.8.1 and searching with later versions of DIA-NN has some changes.  
* [Download Slides](https://plexdia.slavovlab.net/mass-spec/Searching-plexDIA-data-with-DIA-NN.pdf)

<iframe width="560" height="315" src="https://www.youtube.com/embed/0Wmg9LjDtgE" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

&nbsp;


These [detailed methods](https://www.nature.com/articles/s41587-022-01389-w#Sec12) are from: Derks, J., Leduc, A., Wallmann, G. *et al.* Increasing the throughput of sensitive proteomics by plexDIA. *Nat Biotechnol* (2022). [10.1038/s41587-022-01389-w][plexDIA_Nature],  [Preprint][plexDIA_Article], [Nature Research Briefing](https://www.nature.com/articles/s41587-022-01411-1)



[plexDIA_Article]: https://doi.org/10.1101/2021.11.03.467007 "Multiplexed data-independent acquisition by plexDIA"
[plexDIA_Nature]: https://doi.org/10.1038/s41587-022-01389-w "Derks, J., Slavov, N. et al. Increasing the throughput of sensitive proteomics by plexDIA. Nat Biotechnol (2022)"
[plexDIA_Code]: https://github.com/SlavovLab/plexDIA "plexDIA data analysis pipeline, GitHub repository from the Slavov Laboratory"



## Data pipelines for optimizing and processing plexDIA data
A [pipeline][plexDIA_Code] for analyzing plexDIA data and reproducing the analysis by [Derks et al][plexDIA_Nature] are available at the [plexDIA GitHub repository][plexDIA_Code].  


* [Pipeline for processing plexDIA data @ GitHub](https://github.com/SlavovLab/SPP)
* [Data-Driven Optimization of plexDIA by DO-MS](https://do-ms.slavovlab.net/),  [DO-MS @ GitHub](https://github.com/SlavovLab/DO-MS)


-------



&nbsp;  

&nbsp;

&nbsp;  

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;
