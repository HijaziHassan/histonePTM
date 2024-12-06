HistonePTM
<img src="man/figures/logo.png" align="right" width="240" height="277"/>
================
24 November 2024

- [Overview](#overview)
- [Installation](#installation)
- [Contributing](#contributing)
- [Getting help](#getting-help)
- [Workflows](#workflows)

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

## Overview

The goal of `histonePTM` is to make histone PTM analysis less tedious by
offering a whole workflow analysis or functions that help build a
workflow based on whatever software you are using.

Not only this, other functions allow retreiving data from the internet,
manipulate `mgf` files, visualize results in addition to some quality
control assessments.

Some functions rely heavily on other functions from well-established
packages.

## Installation

You can install the development version of histonePTM from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("HijaziHassan/histonePTM")
```

## Contributing

Any contribution is very welcomed. The first version is more adapted to
`Proline` software output. But it tried to generalize each function to
be generic and very flexible to be applicable for other software
outputs.

## Getting help

If you encouter any bug, a problem, a weired behavior, or have a feature
request, please open an
[issue](https://github.com/HijaziHassan/histonePTM/issues).

If you would like to discuss questions related to histone analysis using
mass spectrometry, please open a discussion here
[discussion](https://github.com/HijaziHassan/histonePTM/discussions)

## Workflows

### Analysis of DDA results using `analyzeHistone()`

If you are using [`Proline`](https://www.profiproteomics.fr/proline/)
software to validate identifications resulted from search engines such
as `Mascot`, the function `analyzeHistone()` will:

- **Isolate** histone peptides.
- **Normalize** intensities within peptide families to the total area or
  intensity
- **Abbreviate** histone peptides
- **Rename** PTMs
- **Calculate** coefficients of variations
- **Remove** and **store** duplications
- **Mark** with (**\***) and **store** identifications where the
  software assigns the same peak apex (i.e. same intensities) to
  isobaric positional isomers (i.e. K18acK23un and K18unK23ac) which
  nearly co-elute.
- **Filter** unwanted and/or identifications that are -more often than
  not- false positives (e.g. H3 K37mod).
- **Filter** some identifications if they are not quantified in a
  user-defined number of samples.

This results in 3 Excel file:

- **File 1**: Containing the raw data with sheets. Sheet 1 contains raw
  data of isolated histone peptides without any transformation of data.
  The rest of the sheets are filtered data from the first sheet original
  data. E.g. Only N-terminally labeled peptides.

- **File(s) 2**: Peptide-centric. An excel file per histone protein with
  each sheet containing identifications from same peptide.

- **File 3**: PTM-centric. An excel file summarizing PTMs with each
  sheet containing identifications with specific PTM.

All this with flexibility to:

- choose only to analyze (and output results) of user-defined histone
  protein (e.g. only `H3`).
- filter identification with cut-off threshold of missing values.
- output **File(s) 2** with either removing all unlabeled `me1`,
  `K37mod` (for H3K27R40 peptide) or both.
- group **File(s) 2** into one file or save each protein results in a
  separate file.

**Pre-requisites**

- `Proline` excel output file containing the sheets:

- `Best PSM from protein sets` which includes identification and their
  intensities in each sample. This assumes that IDs with multiple charge
  states are already summed using post-processing functionality inside
  `Proline`.

- `Search settings and infos` which includes information about `RAW`
  files’ names and their corresponding search result files’ names.

- An excel file containing at least two columns:

  - `SampleName`: custom samples names
  - `file`: names of `RAW` files. other optional columns: `Condition`,
    `BioReplicate`, and/or `TechReplicate` depending on the experimental
    design.

Some functions used to build-up this workflow among others are shown
below:

``` r
library(histonePTM)

# analyzeHistone(analysisfile, # file name
#                 metafile, #metafile name
#                 hist_prot= c('All','H3', 'H4', 'H2A', 'H2B'), #choose one these options
#                 labeling = c("PA", "TMA", "PIC_PA", "none") # allow reversing labeling when renaming PTMs
#                 NA_threshold, #numeric #optional
#                 norm_method = c('peptide_family', 'peptide_total'),
#                 extra_filter = c("none", 'no_me1', "K37un", "no_me1_K37un"), #optional
#                 output_result= c('single', 'multiple'), #optional
               
```

### Learn with Fun

Go histonic and get random quotes about histones from the literature!

``` r

litRicher()
#> [1] "Biological databases contribute perspective to evaluate the reasonableness of PTMs."
```

… or do your own literature from Pubmed!

``` r
#install.packages('rentrez')
library(rentrez)
litReview(start = 2023, 
          end = 2024 , 
          db = "pubmed",
          term = "Histone lactylation",
          save_file = FALSE)
#> [[1]]
#> # A tibble: 40 × 3
#>     year id       title                                                         
#>    <int> <chr>    <chr>                                                         
#>  1  2023 38707638 20(S)-ginsenoside Rh2 ameliorates ATRA resistance in APL by m…
#>  2  2023 38327665 Lactate regulates major zygotic genome activation by H3K18 la…
#>  3  2023 38162077 [Role of Histone Modifications in Acute Kidney Injury Progres…
#>  4  2023 38155775 Proteomic analysis identifies PFKP lactylation in SW480 colon…
#>  5  2023 38149461 Natural product fargesin interferes with H3 histone lactylati…
#>  6  2023 38101413 BACH1 changes microglial metabolism and affects astrogenesis …
#>  7  2023 38084701 Histone lactylation-derived LINC01127 promotes the self-renew…
#>  8  2023 38074159 Functioning and mechanisms of PTMs in renal diseases.         
#>  9  2023 38051708 Lactate-induced histone lactylation by p300 promotes osteobla…
#> 10  2023 38041125 Glucose transporter 3 (GLUT3) promotes lactylation modificati…
#> # ℹ 30 more rows
#> 
#> [[2]]
#> 
#> 2023 2024 
#>   20   20
#A list with two tables
    # 1 - Articles titles with their year of publication
    # 2 - number of publication per year
```

scrape MS-identified PTMs from Uniprot

``` r

ptm_Uniprot(Uniprot_accession = 'P62805', #H4_HUMAN
            save_file = FALSE)
#> # A tibble: 7 × 20
#>   type           begin end   xrefs_name xrefs_id    xrefs_url     evidences_code
#>   <chr>          <chr> <chr> <chr>      <chr>       <chr>         <chr>         
#> 1 PROTEOMICS_PTM 69    78    Proteomes  UP000005640 https://www.… ECO:0007829   
#> 2 PROTEOMICS_PTM 69    78    Proteomes  UP000005640 https://www.… ECO:0007829   
#> 3 PROTEOMICS_PTM 81    92    Proteomes  UP000005640 https://www.… ECO:0007829   
#> 4 PROTEOMICS_PTM 25    37    Proteomes  UP000005640 https://www.… ECO:0007829   
#> 5 PROTEOMICS_PTM 81    93    Proteomes  UP000005640 https://www.… ECO:0007829   
#> 6 PROTEOMICS_PTM 47    56    Proteomes  UP000005640 https://www.… ECO:0007829   
#> 7 PROTEOMICS_PTM 46    60    Proteomes  UP000005640 https://www.… ECO:0007829   
#> # ℹ 13 more variables: evidences_source_name <chr>, evidences_source_id <chr>,
#> #   evidences_source_url <chr>, PEP <chr>, peptide <chr>, unique <lgl>,
#> #   name <chr>, position <int>, sources <chr>, dbReferences_name <chr>,
#> #   dbReferences_id <chr>, dbReferences_url <chr>,
#> #   `dbReferences_properties_Localization probability` <chr>
```
