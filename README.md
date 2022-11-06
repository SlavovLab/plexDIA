# **plexDIA: Multiplexed data-independent acquisition**


<!--![GitHub release](https://img.shields.io/github/release/SlavovLab/DO-MS.svg)-->
![GitHub](https://img.shields.io/github/license/SlavovLab/DO-MS.svg)

* [Bulk plexDIA Website](https://plexDIA.slavovlab.net) &nbsp; | &nbsp; [Download data](https://plexDIA.slavovlab.net/mass-spec/data)
* [Single-cell plexDIA Website](https://scp.slavovlab.net/plexDIA) &nbsp; | &nbsp; [Download data](https://scp.slavovlab.net/Derks_et_al_2022)
* [Preprint](https://www.biorxiv.org/content/10.1101/2021.11.03.467007v2) &nbsp; | &nbsp; [*Nature Biotechnology* Article](https://doi.org/10.1038/s41587-022-01389-w)


&nbsp;

<img src="https://scp.slavovlab.net/Figs/plexDIA_4.png" width="70%">


**Perspective:** Derks, J., Slavov N., Strategies for increasing the depth and throughput of protein analysis by plexDIA,  [10.1101/2022.11.05.515287](https://doi.org/10.1101/2022.11.05.515287), [GitHub](https://github.com/SlavovLab/plexDIA_perspective)

&nbsp;

### Requirements

This application has been tested on R >= 3.5.0, OSX 10.14 / Windows 7/8/10. R can be downloaded from the main [R Project page](https://www.r-project.org/) or downloaded with the [RStudio Application](https://www.rstudio.com/products/rstudio/download/).



------------

## Reproducing the data analysis


1. Download all the data reports from the "search" section of MassIVE [MSV000089093](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ae918c7ce5a94a4abd2c6b54a3806c9e).

2. Download the two meta files from the "metadata" section of MassIVE (specifically, "[Meta_SingleCell.tsv](https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000089093/updates/2022-03-23_jderks_10a27189/metadata/Meta_SingleCell.tsv&forceDownload=true)" and "[Meta_Bulk_benchmarking.tsv](https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000089093/updates/2022-03-23_jderks_10a27189/metadata/Meta_Bulk_benchmarking.tsv&forceDownload=true)") [MSV000089093](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ae918c7ce5a94a4abd2c6b54a3806c9e).

3. Create a new R project, and copy the .Rmd, .R, .py, and .txt files from the Code section (https://github.com/SlavovLab/plexDIA) to it. Note: The main functions are in [plexDIA_Functions.R](https://github.com/SlavovLab/plexDIA/blob/main/code/plexDIA_Functions.R), but it uses the other functions in the [code directory](https://github.com/SlavovLab/plexDIA/blob/main/code), so please download all.

4. Please update the file paths in the .Rmd.

5. Run!



## About the project

<!--
DO-MS is a project...


The manuscript for this tool is published at the Journal of Proteome Research: [https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00039](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00039)
-->
The manuscript is freely available on bioRxiv: [Derks et al., 2022 (version 2)](https://www.biorxiv.org/content/10.1101/2021.11.03.467007v2). The peer reviewed version is available at *Nature Biotechnology*: [Derks et al., 2022](https://doi.org/10.1038/s41587-022-01389-w)

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
