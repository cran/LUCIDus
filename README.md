
# LUCIDus: Integreted clustering with multi-view data Version 3.0.1

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/LUCIDus?color=green)](https://cran.r-project.org/package=LUCIDus)
![](https://cranlogs.r-pkg.org/badges/grand-total/LUCIDus?color=blue)
[![](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)
<!-- badges: end -->



The **LUCIDus** package implements the statistical method LUCID proposed in the research paper [A Latent Unknown Clustering Integrating
Multi-Omics Data (LUCID) with Phenotypic Traits](https://doi.org/10.1093/bioinformatics/btz667)
(*Bioinformatics*, 2020). LUCID conducts integrated clustering by using multi-view data, including exposures, and omics data with/without outcome. **LUCIDus** features variable selection, incorporating missingness in omics data, visualization of the LUCID model via Sankey diagram, bootstrap inference, and functions for tuning model parameters.

LUCID version 3.0.1, a major update and enhancement from the original release, implements different integration strategies for multi-omics data with multiple layers including LUCID early integration, LUCID in parallel, and LUCID in serial. It also incorporates methods to deal with missingness in multi-omics data. The following DAG illustrates the three different LUCID models for three integration strategies.

![plot](./figure/fig1.png)

If you are interested in the integration of omic data to estimate mediator or latent structures, please check out [Conti
Lab](https://contilab.usc.edu/about/) to learn more.



## Installation

You can install the development version of LUCIDus 3.0.1 from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ContiLab-usc/LUCIDus-3.0",ref="main",auth_token = "xxx")
```
Note that this repo is now private, so only authorized users can download this package. Please go to [tokens](https://github.com/settings/tokens) to obtain your personal authorized token and input it into auth_token = "xxx" to download this package.

## Workflow
The following figure illustrate the workflow of LUCIDus 3.0.1.
![plot](./figure/fig2.png)


## Usage

Please refer to the
[tutorial](https://USCbiostats.github.io/LUCIDus/articles/LUCIDus.html).

## Citation

    #> 
    #> To cite LUCID methods, please use:
    #> 
    #>   Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi,
    #>   Graham Casey, Duncan C Thomas, David V Conti (2019). A latent unknown
    #>   clustering integrating multi-omics data (LUCID) with phenotypic
    #>   traits. Bioinformatics, btz667. URL
    #>   https://doi.org/10.1093/bioinformatics/btz667
    #> 
    #> To cite LUCIDus R package, please use:
    #> 
    #>   Yinqi Zhao (2022). LUCIDus: an R package to implement the LUCID
    #>   model. CRAN. R package version 2.2 URL
    #>   https://yinqi93.github.io/LUCIDus/index.html
    #> 
    #> To see these entries in BibTeX format, use 'print(<citation>,
    #> bibtex=TRUE)', 'toBibtex(.)', or set
    #> 'options(citation.bibtex.max=999)'.


