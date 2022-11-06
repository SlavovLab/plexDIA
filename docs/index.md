---
layout: default
title: Home
nav_order: 1
description: Sensitive, accurate and high-throughput proteomics by multiplexed data-independent acquisition (plexDIA). plexDIA for increasing proteomics throughput | Slavov Laboratory and single-cell proteomics center
permalink: /
---
{% include social-media-links.html %}

# multi<u>plex</u>ed <u>DIA</u>: *plexDIA*
<!-- {: .fs-6 .fw-300}  {: .fs-9 }   long_title: Multiplexed data-independent acquisition: plexDIA -->

## Increasing the throughput of sensitive proteomics by plexDIA

&nbsp;


[plexDIA Article][plexDIA_Nature]{: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[plexDIA code on GitHub][plexDIA_Code]{: .btn .fs-5 .mb-4 .mb-md-0 }
[Single-cell plexDIA](https://scp.slavovlab.net/plexDIA){: .btn .fs-5 .mb-4 .mb-md-0 }

------------

Current mass-spectrometry methods enable high-throughput proteomics of large sample amounts, but proteomics of low sample amounts remains limited in depth and throughput. To increase the throughput of sensitive proteomics, we developed an experimental and computational framework, [plexDIA][plexDIA_Nature], for simultaneously multiplexing the analysis of both peptides and samples. Multiplexed analysis with plexDIA increases throughput multiplicatively with the number of labels without reducing proteome coverage or quantitative accuracy. By using 3-plex nonisobaric mass tags, plexDIA enables quantifying 3-fold more protein ratios among nanogram-level samples. Using 1 hour active gradients and first-generation Q Exactive, plexDIA quantified about 8,000 proteins in each sample of labeled 3-plex sets. plexDIA also increases data completeness, reducing missing data over 2-fold across samples. When applied to single human cells, plexDIA quantified about 1,000 proteins per cell and achieved 98 % data completeness within a plexDIA set while using about 5 min of active chromatography per cell. These results establish a general framework for increasing the throughput of sensitive and quantitative protein analysis.


* Derks, J., Leduc, A., Wallmann, G. *et al.* Increasing the throughput of sensitive proteomics by plexDIA. *Nat Biotechnol* (2022). [10.1038/s41587-022-01389-w][plexDIA_Nature],  [Preprint][plexDIA_Article]
* Derks, J., Slavov N., [Strategies for increasing the depth and throughput of protein analysis by plexDIA](https://www.biorxiv.org/content/10.1101/2022.11.05.515287v1),  *bioRxiv* (2022). [10.1101/2022.11.05.515287](https://doi.org/10.1101/2022.11.05.515287), [GitHub](https://github.com/SlavovLab/plexDIA_perspective)

------------


[![plexDIA: Multiplexed data-independent acquisition for increasing proteomics throughput](https://scp.slavovlab.net/Figs/plexDIA_4.png){: width="100%" .center-image}](https://scp.slavovlab.net/plexDIA)


------------



## Perspectives on high-throughput multiplexed proteomics
* [Framework for multiplicative scaling of single-cell proteomics](https://www.nature.com/articles/s41587-022-01411-1), *Nature Biotechnology*
* [Increasing proteomics throughput](https://www.nature.com/articles/s41587-021-00881-z), *Nature Biotechnology*
* [Driving Single Cell Proteomics Forward with Innovation](https://pubmed.ncbi.nlm.nih.gov/34597050/), *Journal of Proteome Research*
* [Scaling up single-cell proteomics](https://doi.org/10.1016/j.mcpro.2021.100179), *Molecular and Cellular Proteomics*



## About the project

plexDIA is a project developed in the [Slavov Laboratory](http://slavovlab.net) at [Northeastern University](https://www.northeastern.edu/) in collaboration with Demichev and Rasler Laboratories at Charité, Universitätsmedizin. It was authored by [Jason Derks](https://slavovlab.net/people.htm), [Andrew Leduc](http://andrewdleduc.com/), [Harrison Specht](http://harrisonspecht.com), [R. Gray Huffman](https://slavovlab.net/people.htm), [Markus Ralser](https://www.crick.ac.uk/research/labs/markus-ralser), [Vadim Demichev](https://github.com/vdemichev) and [Nikolai Slavov](https://coe.northeastern.edu/people/slavov-nikolai/).   


Contact the authors by email: [nslavov\{at\}northeastern.edu](mailto:nslavov@northeastern.edu).

This project was supported by funding from the [NIH Director's Award](https://projectreporter.nih.gov/project_info_description.cfm?aid=9167004&icde=31336575) and by an [Allen Distinguished Investigator Award](https://alleninstitute.org/what-we-do/frontiers-group/distinguished-investigators/projects/tracking-proteome-dynamics-single-cells) from the Paul G. Allen Frontiers Group.


[plexDIA_Article]: https://doi.org/10.1101/2021.11.03.467007 "Multiplexed data-independent acquisition by plexDIA"
[plexDIA_Nature]: https://doi.org/10.1038/s41587-022-01389-w "Derks, J., Slavov, N. et al. Increasing the throughput of sensitive proteomics by plexDIA. Nat Biotechnol (2022)"
[plexDIA_Code]: https://github.com/SlavovLab/plexDIA "plexDIA data analysis pipeline repository"
