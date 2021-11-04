---
layout: default
title: Home
nav_order: 1
description: "plexDIA: Multiplexed data-independent acquisition for increasing proteomics throughput"
permalink: /
---
{% include social-media-links.html %}

# multi<u>plex</u>ed <u>DIA</u>: *plexDIA*
<!-- {: .fs-6 .fw-300}  {: .fs-9 } -->

## Increasing the throughput of sensitive proteomics by plexDIA

&nbsp;


[plexDIA Preprint][preprint_link]{: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[plexDIA code on GitHub][github_link]{: .btn .fs-5 .mb-4 .mb-md-0 }

------------

Current mass-spectrometry methods enable high-throughput proteomics of large sample amounts, but proteomics of low sample amounts remains limited in depth and throughput. We aimed to increase throughput for analyzing limited samples while achieving high proteome coverage and quantitative accuracy. We developed a general experimental and computational framework, plexDIA, for simultaneously multiplexing the analysis of both peptides and samples. Multiplexed analysis with plexDIA increases throughput multiplicatively with the number of labels without reducing proteome coverage or quantitative accuracy. By using using 3-plex nonisobaric mass tags, plexDIA enables quantifying 3-fold more protein ratios among nanogram-level samples. Using 1 hour active gradients and first-generation Q Exactive, plexDIA quantified about 8,000 proteins in each sample of labeled 3-plex sets. Furthermore, plexDIA increases the consistency of protein quantification, resulting in over 2-fold reduction of missing data across samples. We applied plexDIA to quantify proteome dynamics during the cell division cycle in cells isolated based on their DNA content. The high sensitivity and accuracy of plexDIA detected many classical cell cycle proteins and discovered new ones. These results establish a general framework for increasing the throughput of highly sensitive and quantitative protein analysis.  

------------


![plexDIA: Multiplexed data-independent acquisition for increasing proteomics throughput]({{site.baseurl}}/mass-spec/Figures/plexDIA.png){: width="100%" .center-image}

------------



## Perspectives on high-throughput multiplexed proteomics
* [Increasing proteomics throughput](https://www.nature.com/articles/s41587-021-00881-z), *Nature Biotechnology*
* [Driving Single Cell Proteomics Forward with Innovation](https://pubmed.ncbi.nlm.nih.gov/34597050/), *Journal of Proteome Research*



## About the project

plexDIA is a project developed in the [Slavov Laboratory](http://slavovlab.net) at [Northeastern University](https://www.northeastern.edu/) in collaboration with Demichev and Rasler Laboratories at Charité, Universitätsmedizin. It was authored by [Jason Derks](https://slavovlab.net/people.htm), [Andrew Leduc](http://andrewdleduc.com/), [Harrison Specht](http://harrisonspecht.com), [R. Gray Huffman](https://slavovlab.net/people.htm), [Markus Ralser](https://www.crick.ac.uk/research/labs/markus-ralser), [Vadim Demichev](https://github.com/vdemichev) and [Nikolai Slavov](https://coe.northeastern.edu/people/slavov-nikolai/).   


Contact the authors by email: [nslavov\{at\}northeastern.edu](mailto:nslavov@northeastern.edu).

This project was supported by funding from the [NIH Director's Award](https://projectreporter.nih.gov/project_info_description.cfm?aid=9167004&icde=31336575) and by an [Allen Distinguished Investigator Award](https://alleninstitute.org/what-we-do/frontiers-group/distinguished-investigators/projects/tracking-proteome-dynamics-single-cells) from the Paul G. Allen Frontiers Group.


[plexDIA_Article]: preprint_link "Multiplexed data-independent acquisition by plexDIA"
[plexDIA_Code]: github_link "plexDIA data analysis pipeline repository"
