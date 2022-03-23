# **plexDIA: Multiplexed data-independent acquisition**


<!--![GitHub release](https://img.shields.io/github/release/SlavovLab/DO-MS.svg)-->
![GitHub](https://img.shields.io/github/license/SlavovLab/DO-MS.svg)

* [Bulk plexDIA Website](https://plexDIA.slavovlab.net) | [Download data](https://plexDIA.slavovlab.net/mass-spec/data)
* [Single-cell plexDIA Website](https://scp.slavovlab.net/plexDIA) | [Download data](https://scp.slavovlab.net/Derks_et_al_2022)
* [Preprint](https://www.biorxiv.org/content/10.1101/2021.11.03.467007v2)


&nbsp;

<img src="https://scp.slavovlab.net/Figs/plexDIA_4.png" width="70%">



### Requirements

This application has been tested on R >= 3.5.0, OSX 10.14 / Windows 7/8/10. R can be downloaded from the main [R Project page](https://www.r-project.org/) or downloaded with the [RStudio Application](https://www.rstudio.com/products/rstudio/download/).



------------

## Reproducing the analysis

The main function is [plexDIA_Functions.R](https://github.com/SlavovLab/plexDIA/blob/main/code/plexDIA_Functions.R), and it uses the other functions in the [code directory](https://github.com/SlavovLab/plexDIA/blob/main/code/p).  


A. We have three main DATA directories which contain raw data files, DIA-NN reports, and data analysis results for each topic covered in the plexDIA paper:

	1) “DATA_Benchmarking”. This folder contains data from H. sapiens (U-937 and Jurkat), S. cerevisiae, and E. coli mixed-species experiments which was used to generate Figures 2-4 and Supplementary Figures 2-7.

	2) “DATA_Cell_Division_Cycle”. This folder contains data from U-937 (monocytes) which were FACS-sorted by cell-cycle phase, then combined into a plexDIA set to be run with V1 and V2 DIA methods. This data was used to generate Figure 5.

	3) “DATA_Single_cell”. This folder contains data from PDAC, Melanoma, and U-937 single cell plexDIA experiments, as well as 100-cell bulk runs which were used for benchmarking quantitative accuracy and generating a spectral library for searching single-cell data. This data was used to generate Figure 6 and Supplementary Figure 8 and 9.

	Note: each of these folders contains a "DIANN_outputs" subfolder, with all of the outputs that DIA-NN reports. The most relevant outputs for re-analysis will the Report.tsv. For more information about DIA-NN and its outputs please visit: https://github.com/vdemichev/DiaNN. For simplified DIA-NN reports which were used for generating figures and analysis, please download from MassIVE https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ae918c7ce5a94a4abd2c6b54a3806c9e (MSV000089093).


B. Installer for DIA-NN version 1.8.1 beta 16, which was used for this analysis.


C. Meta data for all raw mass-spec data can be found in “RawFile_Info.txt".


D. FASTA files used for searches can be found in the "FASTAs" folder. These FASTA files are from Swiss-Prot and contain canonical and isoform proteins, downloaded 2022.


E. Spectral libraries used for DIA-NN searches can be found in the "Spectral_libraries" folder.


F. DIA-NN search and library pipelines can be found in the "DIANN_pipelines" folder. These can be used to re-search raw data with the same settings that were used for generating the data.


G. Code used for generating figures and performing data analysis can be found here: https://github.com/SlavovLab/plexDIA









## About the project

<!--
DO-MS is a project...


The manuscript for this tool is published at the Journal of Proteome Research: [https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00039](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00039)
-->
The manuscript is freely available on bioRxiv: [Derks et al., 2022 (version 2)](https://www.biorxiv.org/content/10.1101/2021.11.03.467007v2)

For more information, contact [Slavov Laboratory](https://slavovlab.net) or directly [Prof. Nikolai Slavov](https://coe.northeastern.edu/people/slavov-nikolai/)

### License

The plexDIA code is distributed by an [MIT license](https://github.com/SlavovLab/DO-MS/blob/master/LICENSE).

### Contributing

Please feel free to contribute to this project by opening an issue or pull request.

<!--
### Data
All data used for the manuscript is available on [UCSD's MassIVE Repository](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ed5a1ab37dc34985bbedbf3d9a945535)
-->

<!--
### Figures/Analysis
Scripts for the figures in the DART-ID manuscript are available in a separate GitHub repository, [https://github.com/SlavovLab/DART-ID_2018](https://github.com/SlavovLab/DART-ID_2018)
-->

-------------

## Help!

For any bugs, questions, or feature requests,
please use the [GitHub issue system](https://github.com/SlavovLab/plexDIA/issues) to contact the developers.
