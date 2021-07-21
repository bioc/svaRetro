
<!-- README.md is generated from README.Rmd. Please do not edit this file directly. -->

# svaNUMT: R package for retrotransposed transcript detection from structural variants

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- badges: end -->

`svaRetro` contains functions for detecting retrotransposed transcripts
(RTs) from structural variant calls. It takes structural variant calls
in GRanges of breakend notation and identifies RTs by exon-exon
junctions and insertion sites. The candidate RTs are reported by events
and annotated with information of the inserted transcripts.

This package uses a breakend-centric event notation adopted from the
[`StructuralVariantAnnotation`](https://www.bioconductor.org/packages/release/bioc/html/StructuralVariantAnnotation.html)
package. More information about `VCF` objects and breakend-centric
GRanges object can be found by consulting the vignettes in the
corresponding packages with `browseVignettes("VariantAnnotation")` and
`browseVignettes("StructuralVariantAnnotation")`.

# Installation

`svaRetro` is currently under review on Bioconductor.

<!-- after acceptance 
[svaNUMT](https://bioconductor.org/packages/svaRetro) is currently available for download in the 'Devel' version of Bioconductor:


```r
# install.packages("BiocManager")
BiocManager::install("svaRetro")
```
-->

The development version can be installed from GitHub:

``` r
BiocManager::install("PapenfussLab/svaRetro")
```

# Workflow

Below is a workflow example for detecting RTs from a human SV callset.
This example is taken from the vignette of `svaRetro`.

``` r
library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(svaRetro)

RT_vcf <- readVcf(system.file("extdata", "diploidSV.vcf", package = "svaRetro"))
```

``` r
RT_gr <- StructuralVariantAnnotation::breakpointRanges(RT_vcf, nominalPosition=TRUE)
head(RT_gr)
#> GRanges object with 6 ranges and 12 metadata columns:
#>                                seqnames    ranges strand | paramRangeID
#>                                   <Rle> <IRanges>  <Rle> |     <factor>
#>   MantaINS:0:775:775:0:1:0_bp1        1     66365      + |           NA
#>      MantaINS:35:0:0:0:0:0_bp1        1   1004204      + |           NA
#>      MantaDEL:92:0:0:0:0:0_bp1        1   1161716      + |           NA
#>     MantaDEL:127:0:0:0:0:0_bp1        1   1162672      + |           NA
#>     MantaDEL:130:0:0:0:0:0_bp1        1   1183434      + |           NA
#>     MantaDEL:107:0:0:0:1:0_bp1        1   1302326      + |           NA
#>                                                                                                                                                                                                                                                                                                     REF
#>                                                                                                                                                                                                                                                                                             <character>
#>   MantaINS:0:775:775:0:1:0_bp1                                                                                                                                                                                                                                                           AATATAATATATAA
#>      MantaINS:35:0:0:0:0:0_bp1                                                                                                                                                                                                                                                                        G
#>      MantaDEL:92:0:0:0:0:0_bp1                                                                                                                                                                                                        CCTGTACGGTCAGGAGGAAACATGGCACCTCCCCTCTGGGGGCTCTTTCCAGAAACCCTCAACCC
#>     MantaDEL:127:0:0:0:0:0_bp1                                                                                          GGCGGGAAGGCGAGCTCGTGGCCAGGCCCTGCGGGAAGGCGAGCTCGTGGCCAGGCCCGGCGGGAAGGCGAGCTCGTGGCCAGGCCCGGCGGGAAGGCGAGCTCGTGGCCAGGCCCGGCGGGAAGGCGAGCTCGTGGCCAGGCCCTGCGGGAAGGCGAGCTCGTGGCCAGGCCCT
#>     MantaDEL:130:0:0:0:0:0_bp1 CAGGCTGGATCTCCAACTCTGACCTACAGGCAGGAAAGTGGGCAGCCCTGGGAGGCTGGACTGAGGGAGGCTGGACTTCCCACTCAGGCCTACACGCAGGAAAATGGGCAGCCCTGGGAGGCTGGACCGAGGGAGGCTGGGCCTCCCACTCCACCCTACAGGCCAGGACACGGGCAGCCCTGGGAGGCTAGACCGAGGGAGGCTGGGCCTCCCATCTACCCTACAGGCCGGGACACAGGCAGCCCTGGGAGGCTGTACCGAGGG
#>     MantaDEL:107:0:0:0:1:0_bp1                                                          GAATGAGTGGATTGGTGAGTGAATTGGTGAGTTGAATTGGTGTGTGTAGTGGATGAGTGTGGATGAATGTGAATTGGCGAGTATGGATGTGTGAATTGGTGAGTGTGAATGTGTGGATTGGTGAGTGAATTGGTGAGTTGAATTGGTGTGTGTAGTGTGGATGAGTGTGAATTGGCGAGTGTGGATGAGTGTGAATTGGTGAGTGTG
#>                                                                                                                                                                            ALT
#>                                                                                                                                                                    <character>
#>   MantaINS:0:775:775:0:1:0_bp1 ATATATATATTATTATATAATATATATTATATAATATATTTTATTATATAATATAATATATATTATATAATATAATATATTTTATTATATAAATATATATTATATTATATAATATAATATATATTAATATAAATATATATTAT
#>      MantaINS:35:0:0:0:0:0_bp1                                                                 GGCCACGCGGGCTGTGCAGATGCAGGTGCGGCGGGGCGGGGCCACGCGGGCTGTGAAGGTGCAGGTGCGGCGGGGCAGA
#>      MantaDEL:92:0:0:0:0:0_bp1                                                                                                                                              CT
#>     MantaDEL:127:0:0:0:0:0_bp1                                                                                                                                               G
#>     MantaDEL:130:0:0:0:0:0_bp1                                                                                                                                               C
#>     MantaDEL:107:0:0:0:1:0_bp1                                                                                                                                      GCAGTGTGAA
#>                                     QUAL      FILTER                 sourceId
#>                                <numeric> <character>              <character>
#>   MantaINS:0:775:775:0:1:0_bp1       999    MaxDepth MantaINS:0:775:775:0:1:0
#>      MantaINS:35:0:0:0:0:0_bp1       999        PASS    MantaINS:35:0:0:0:0:0
#>      MantaDEL:92:0:0:0:0:0_bp1       999        PASS    MantaDEL:92:0:0:0:0:0
#>     MantaDEL:127:0:0:0:0:0_bp1       440        PASS   MantaDEL:127:0:0:0:0:0
#>     MantaDEL:130:0:0:0:0:0_bp1       643        PASS   MantaDEL:130:0:0:0:0:0
#>     MantaDEL:107:0:0:0:1:0_bp1       999        PASS   MantaDEL:107:0:0:0:1:0
#>                                                     partner      svtype
#>                                                 <character> <character>
#>   MantaINS:0:775:775:0:1:0_bp1 MantaINS:0:775:775:0:1:0_bp2         INS
#>      MantaINS:35:0:0:0:0:0_bp1    MantaINS:35:0:0:0:0:0_bp2         INS
#>      MantaDEL:92:0:0:0:0:0_bp1    MantaDEL:92:0:0:0:0:0_bp2         DEL
#>     MantaDEL:127:0:0:0:0:0_bp1   MantaDEL:127:0:0:0:0:0_bp2         DEL
#>     MantaDEL:130:0:0:0:0:0_bp1   MantaDEL:130:0:0:0:0:0_bp2         DEL
#>     MantaDEL:107:0:0:0:1:0_bp1   MantaDEL:107:0:0:0:1:0_bp2         DEL
#>                                    svLen
#>                                <numeric>
#>   MantaINS:0:775:775:0:1:0_bp1       129
#>      MantaINS:35:0:0:0:0:0_bp1        78
#>      MantaDEL:92:0:0:0:0:0_bp1       -63
#>     MantaDEL:127:0:0:0:0:0_bp1      -174
#>     MantaDEL:130:0:0:0:0:0_bp1      -263
#>     MantaDEL:107:0:0:0:1:0_bp1      -197
#>                                                                                                                                                                        insSeq
#>                                                                                                                                                                   <character>
#>   MantaINS:0:775:775:0:1:0_bp1 TATATATATTATTATATAATATATATTATATAATATATTTTATTATATAATATAATATATATTATATAATATAATATATTTTATTATATAAATATATATTATATTATATAATATAATATATATTAATATAAATATATATTAT
#>      MantaINS:35:0:0:0:0:0_bp1                                                                 GCCACGCGGGCTGTGCAGATGCAGGTGCGGCGGGGCGGGGCCACGCGGGCTGTGAAGGTGCAGGTGCGGCGGGGCAGA
#>      MantaDEL:92:0:0:0:0:0_bp1                                                                                                                                              T
#>     MantaDEL:127:0:0:0:0:0_bp1                                                                                                                                               
#>     MantaDEL:130:0:0:0:0:0_bp1                                                                                                                                               
#>     MantaDEL:107:0:0:0:1:0_bp1                                                                                                                                      CAGTGTGAA
#>                                   insLen    HOMLEN
#>                                <numeric> <numeric>
#>   MantaINS:0:775:775:0:1:0_bp1       142         0
#>      MantaINS:35:0:0:0:0:0_bp1        78        10
#>      MantaDEL:92:0:0:0:0:0_bp1         1         0
#>     MantaDEL:127:0:0:0:0:0_bp1         0         9
#>     MantaDEL:130:0:0:0:0:0_bp1         0         7
#>     MantaDEL:107:0:0:0:1:0_bp1         9         0
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome
```

Note that `StructuralVariantAnnotation` requires the `GRanges` object to
be composed entirely of valid breakpoints. Please consult the vignette
of the `StructuralVariantAnnotation` package for ensuring breakpoint
consistency.

### Identifying Retrotransposed Transcripts

The package provides `rtDetect` to identify RTs using the provided SV
calls. This is achieved by detecting intronic deletions, which are
breakpoints at exon-intron (and intron-exon) boundaries of a transcript.
Fusions consisting of an exon boundary and a second genomic location are
reported as potential insertion sites. Due to the complexity of RT
events, insertion sites can be discovered on both left and right sides,
only one side, or none at all.

``` r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#> Loading required package: GenomicFeatures
#> Loading required package: AnnotationDbi
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:AnnotationDbi':
#> 
#>     select
#> The following object is masked from 'package:VariantAnnotation':
#> 
#>     select
#> The following objects are masked from 'package:Biostrings':
#> 
#>     collapse, intersect, setdiff, setequal, union
#> The following object is masked from 'package:XVector':
#> 
#>     slice
#> The following object is masked from 'package:matrixStats':
#> 
#>     count
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following objects are masked from 'package:GenomicRanges':
#> 
#>     intersect, setdiff, union
#> The following object is masked from 'package:GenomeInfoDb':
#> 
#>     intersect
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
hg19.genes <- TxDb.Hsapiens.UCSC.hg19.knownGene
RT_vcf <- readVcf(system.file("extdata", "diploidSV.vcf", package = "svaRetro"))
RT_gr <- StructuralVariantAnnotation::breakpointRanges(RT_vcf, nominalPosition=TRUE)
RT <- rtDetect(RT_gr, hg19.genes, maxgap=50, minscore=0.3)
```

The output is a list of `GRanges` object consisting of two sets of
`GRanges` calls, `insSite` and `junctions`, containing candidate
insertion sites and exon-exon junctions respectively. Candidate
insertion sites are annotated by the source transcripts and whether
exon-exon junctions are detected for the source transcripts. RT junction
breakends are annotated by the UCSC exon IDs, corresponding transcripts,
and NCBI gene symbols.

``` r
RT$SKA3
#> $junctions
#> GRanges object with 14 ranges and 16 metadata columns:
#>                                  seqnames    ranges strand | paramRangeID
#>                                     <Rle> <IRanges>  <Rle> |     <factor>
#>    MantaDEL:245251:6:6:0:0:0_bp2       13  21729832      - |           NA
#>    MantaDEL:245251:5:8:0:0:0_bp2       13  21732061      - |           NA
#>    MantaDEL:245251:5:9:0:0:0_bp2       13  21734038      - |           NA
#>   MantaDEL:245251:7:10:0:0:0_bp2       13  21735929      - |           NA
#>   MantaDEL:245251:4:11:0:0:0_bp2       13  21742127      - |           NA
#>                              ...      ...       ...    ... .          ...
#>    MantaDEL:245251:5:9:0:0:0_bp1       13  21732261      + |           NA
#>   MantaDEL:245251:7:10:0:0:0_bp1       13  21734126      + |           NA
#>   MantaDEL:245251:4:11:0:0:0_bp1       13  21736014      + |           NA
#>    MantaDEL:245251:3:4:0:0:0_bp1       13  21742538      + |           NA
#>    MantaDEL:245251:2:3:0:0:0_bp1       13  21746642      + |           NA
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              REF
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      <character>
#>    MantaDEL:245251:6:6:0:0:0_bp2 TCTGCAACAGATACAAATAACAAATATCAATTTAATAAAATTAAAAGCCATTAAGACAAATGACACAATACTGTGGCTATATATTTTACACTTATAAAATAATTGAGGATAGATTCCCACTGATATCATTAAACTGGATAATTCGGGAATCTGAGATTCAGGGATCACAAGTTCTATATCAAAAGATAGAGACAGGCTATTAACTTAAGCTGGCAAATGTCAATTAAAAACAAAATTTTTACCAATATTCAAATGTTAATTTTTTTTTTTTTTTTTTAGATGGAGTCTCGCTCTGTTGCCAGGCTGGAGTGCAGTGGCATGATCTTGGCTCACTGCAACCTCCGCCTCCCAGATTCAAGCGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGACTACAGGTGTGCACCACCACACCTGGCTAATTTTTGTATTTTTAGTAGAGACGAGGTTTCACCATGTTGGTCAGGATGGTCAAATGTTAATTTTTAAATGTCCTCCTCAAATAACACATGAACTTTCTTTACAAAGGTAACATACTCAC
#>    MantaDEL:245251:5:8:0:0:0_bp2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               G
#>    MantaDEL:245251:5:9:0:0:0_bp2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               A
#>   MantaDEL:245251:7:10:0:0:0_bp2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               T
#>   MantaDEL:245251:4:11:0:0:0_bp2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               A
#>                              ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ...
#>    MantaDEL:245251:5:9:0:0:0_bp1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               A
#>   MantaDEL:245251:7:10:0:0:0_bp1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               T
#>   MantaDEL:245251:4:11:0:0:0_bp1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               A
#>    MantaDEL:245251:3:4:0:0:0_bp1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               A
#>    MantaDEL:245251:2:3:0:0:0_bp1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               T
#>                                          ALT      QUAL      FILTER
#>                                  <character> <numeric> <character>
#>    MantaDEL:245251:6:6:0:0:0_bp2           T       999        PASS
#>    MantaDEL:245251:5:8:0:0:0_bp2       <DEL>       999        PASS
#>    MantaDEL:245251:5:9:0:0:0_bp2       <DEL>       525        PASS
#>   MantaDEL:245251:7:10:0:0:0_bp2       <DEL>       539        PASS
#>   MantaDEL:245251:4:11:0:0:0_bp2       <DEL>       999        PASS
#>                              ...         ...       ...         ...
#>    MantaDEL:245251:5:9:0:0:0_bp1       <DEL>       525        PASS
#>   MantaDEL:245251:7:10:0:0:0_bp1       <DEL>       539        PASS
#>   MantaDEL:245251:4:11:0:0:0_bp1       <DEL>       999        PASS
#>    MantaDEL:245251:3:4:0:0:0_bp1       <DEL>       999        PASS
#>    MantaDEL:245251:2:3:0:0:0_bp1       <DEL>       999        PASS
#>                                                    sourceId
#>                                                 <character>
#>    MantaDEL:245251:6:6:0:0:0_bp2  MantaDEL:245251:6:6:0:0:0
#>    MantaDEL:245251:5:8:0:0:0_bp2  MantaDEL:245251:5:8:0:0:0
#>    MantaDEL:245251:5:9:0:0:0_bp2  MantaDEL:245251:5:9:0:0:0
#>   MantaDEL:245251:7:10:0:0:0_bp2 MantaDEL:245251:7:10:0:0:0
#>   MantaDEL:245251:4:11:0:0:0_bp2 MantaDEL:245251:4:11:0:0:0
#>                              ...                        ...
#>    MantaDEL:245251:5:9:0:0:0_bp1  MantaDEL:245251:5:9:0:0:0
#>   MantaDEL:245251:7:10:0:0:0_bp1 MantaDEL:245251:7:10:0:0:0
#>   MantaDEL:245251:4:11:0:0:0_bp1 MantaDEL:245251:4:11:0:0:0
#>    MantaDEL:245251:3:4:0:0:0_bp1  MantaDEL:245251:3:4:0:0:0
#>    MantaDEL:245251:2:3:0:0:0_bp1  MantaDEL:245251:2:3:0:0:0
#>                                                         partner      svtype
#>                                                     <character> <character>
#>    MantaDEL:245251:6:6:0:0:0_bp2  MantaDEL:245251:6:6:0:0:0_bp1         DEL
#>    MantaDEL:245251:5:8:0:0:0_bp2  MantaDEL:245251:5:8:0:0:0_bp1         DEL
#>    MantaDEL:245251:5:9:0:0:0_bp2  MantaDEL:245251:5:9:0:0:0_bp1         DEL
#>   MantaDEL:245251:7:10:0:0:0_bp2 MantaDEL:245251:7:10:0:0:0_bp1         DEL
#>   MantaDEL:245251:4:11:0:0:0_bp2 MantaDEL:245251:4:11:0:0:0_bp1         DEL
#>                              ...                            ...         ...
#>    MantaDEL:245251:5:9:0:0:0_bp1  MantaDEL:245251:5:9:0:0:0_bp2         DEL
#>   MantaDEL:245251:7:10:0:0:0_bp1 MantaDEL:245251:7:10:0:0:0_bp2         DEL
#>   MantaDEL:245251:4:11:0:0:0_bp1 MantaDEL:245251:4:11:0:0:0_bp2         DEL
#>    MantaDEL:245251:3:4:0:0:0_bp1  MantaDEL:245251:3:4:0:0:0_bp2         DEL
#>    MantaDEL:245251:2:3:0:0:0_bp1  MantaDEL:245251:2:3:0:0:0_bp2         DEL
#>                                      svLen      insSeq    insLen    HOMLEN
#>                                  <numeric> <character> <numeric> <numeric>
#>    MantaDEL:245251:6:6:0:0:0_bp2      -542                     0         1
#>    MantaDEL:245251:5:8:0:0:0_bp2     -2110        <NA>         0         2
#>    MantaDEL:245251:5:9:0:0:0_bp2     -1776        <NA>         0         4
#>   MantaDEL:245251:7:10:0:0:0_bp2     -1802        <NA>         0         1
#>   MantaDEL:245251:4:11:0:0:0_bp2     -6112        <NA>         0         2
#>                              ...       ...         ...       ...       ...
#>    MantaDEL:245251:5:9:0:0:0_bp1     -1776        <NA>         0         4
#>   MantaDEL:245251:7:10:0:0:0_bp1     -1802        <NA>         0         1
#>   MantaDEL:245251:4:11:0:0:0_bp1     -6112        <NA>         0         2
#>    MantaDEL:245251:3:4:0:0:0_bp1     -3939        <NA>         0         2
#>    MantaDEL:245251:2:3:0:0:0_bp1     -3870        <NA>         0         2
#>                                       exon                              txs
#>                                  <integer>                           <list>
#>    MantaDEL:245251:6:6:0:0:0_bp2    176912            uc001unt.3,uc001unv.3
#>    MantaDEL:245251:5:8:0:0:0_bp2    176913            uc001unt.3,uc001unv.3
#>    MantaDEL:245251:5:9:0:0:0_bp2    176914 uc001unt.3,uc001unu.3,uc001unv.3
#>   MantaDEL:245251:7:10:0:0:0_bp2    176915 uc001unt.3,uc001unu.3,uc001unv.3
#>   MantaDEL:245251:4:11:0:0:0_bp2    176916 uc001unt.3,uc001unu.3,uc001unv.3
#>                              ...       ...                              ...
#>    MantaDEL:245251:5:9:0:0:0_bp1    176913 uc001unt.3,uc001unu.3,uc001unv.3
#>   MantaDEL:245251:7:10:0:0:0_bp1    176914 uc001unt.3,uc001unu.3,uc001unv.3
#>   MantaDEL:245251:4:11:0:0:0_bp1    176915 uc001unt.3,uc001unu.3,uc001unv.3
#>    MantaDEL:245251:3:4:0:0:0_bp1    176916 uc001unt.3,uc001unu.3,uc001unv.3
#>    MantaDEL:245251:2:3:0:0:0_bp1    176917 uc001unt.3,uc001unu.3,uc001unv.3
#>                                   exons gene_symbol
#>                                  <list>      <list>
#>    MantaDEL:245251:6:6:0:0:0_bp2 176912        SKA3
#>    MantaDEL:245251:5:8:0:0:0_bp2 176913        SKA3
#>    MantaDEL:245251:5:9:0:0:0_bp2 176914        SKA3
#>   MantaDEL:245251:7:10:0:0:0_bp2 176915        SKA3
#>   MantaDEL:245251:4:11:0:0:0_bp2 176916        SKA3
#>                              ...    ...         ...
#>    MantaDEL:245251:5:9:0:0:0_bp1 176913        SKA3
#>   MantaDEL:245251:7:10:0:0:0_bp1 176914        SKA3
#>   MantaDEL:245251:4:11:0:0:0_bp1 176915        SKA3
#>    MantaDEL:245251:3:4:0:0:0_bp1 176916        SKA3
#>    MantaDEL:245251:2:3:0:0:0_bp1 176917        SKA3
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome
#> 
#> $insSite
#> GRanges object with 4 ranges and 17 metadata columns:
#>                                 seqnames    ranges strand | paramRangeID
#>                                    <Rle> <IRanges>  <Rle> |     <factor>
#>     MantaBND:245251:0:3:0:0:0:0       13  21746762      + |           NA
#>   MantaDEL:245251:5:6:0:0:0_bp2       13  21731995      - |           NA
#>     MantaBND:245251:0:3:0:0:0:1       11 108585702      - |           NA
#>   MantaDEL:245251:5:6:0:0:0_bp1       13  21729260      + |           NA
#>                                         REF             ALT      QUAL
#>                                 <character>     <character> <numeric>
#>     MantaBND:245251:0:3:0:0:0:0           T T[11:108585702[        49
#>   MantaDEL:245251:5:6:0:0:0_bp2           T           <DEL>       283
#>     MantaBND:245251:0:3:0:0:0:1           T  ]13:21746762]T        49
#>   MantaDEL:245251:5:6:0:0:0_bp1           T           <DEL>       283
#>                                      FILTER                    sourceId
#>                                 <character>                 <character>
#>     MantaBND:245251:0:3:0:0:0:0        PASS MantaBND:245251:0:3:0:0:0:0
#>   MantaDEL:245251:5:6:0:0:0_bp2        PASS   MantaDEL:245251:5:6:0:0:0
#>     MantaBND:245251:0:3:0:0:0:1        PASS MantaBND:245251:0:3:0:0:0:1
#>   MantaDEL:245251:5:6:0:0:0_bp1        PASS   MantaDEL:245251:5:6:0:0:0
#>                                                       partner      svtype
#>                                                   <character> <character>
#>     MantaBND:245251:0:3:0:0:0:0   MantaBND:245251:0:3:0:0:0:1         BND
#>   MantaDEL:245251:5:6:0:0:0_bp2 MantaDEL:245251:5:6:0:0:0_bp1         DEL
#>     MantaBND:245251:0:3:0:0:0:1   MantaBND:245251:0:3:0:0:0:0         BND
#>   MantaDEL:245251:5:6:0:0:0_bp1 MantaDEL:245251:5:6:0:0:0_bp2         DEL
#>                                     svLen      insSeq    insLen    HOMLEN
#>                                 <numeric> <character> <numeric> <numeric>
#>     MantaBND:245251:0:3:0:0:0:0        NA                     0         0
#>   MantaDEL:245251:5:6:0:0:0_bp2     -2734        <NA>         0         0
#>     MantaBND:245251:0:3:0:0:0:1        NA                     0         0
#>   MantaDEL:245251:5:6:0:0:0_bp1     -2734        <NA>         0         0
#>                                  exons                              txs
#>                                 <list>                           <list>
#>     MantaBND:245251:0:3:0:0:0:0 176918            uc001unt.3,uc001unu.3
#>   MantaDEL:245251:5:6:0:0:0_bp2 176911 uc001unt.3,uc001unu.3,uc001unv.3
#>     MantaBND:245251:0:3:0:0:0:1     NA                               NA
#>   MantaDEL:245251:5:6:0:0:0_bp1     NA                               NA
#>                                        rtFound rtFoundSum gene_symbol
#>                                         <list>  <logical>      <list>
#>     MantaBND:245251:0:3:0:0:0:0      TRUE,TRUE       TRUE        SKA3
#>   MantaDEL:245251:5:6:0:0:0_bp2 TRUE,TRUE,TRUE       TRUE        SKA3
#>     MantaBND:245251:0:3:0:0:0:1             NA       <NA>          NA
#>   MantaDEL:245251:5:6:0:0:0_bp1             NA       <NA>          NA
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome
```

## Visualising breakpoint pairs via circos plots

One way of visualising RT is by circos plots. Here we use the package
[`circlize`](https://doi.org/10.1093/bioinformatics/btu393) to
demonstrate the visualisation of insertion site and exon-exon junctions.

To generate a simple circos plot of RT event with SKA3 transcript:

``` r
library(circlize)
rt_chr_prefix <- c(RT$SKA3$junctions, RT$SKA3$insSite)
seqlevelsStyle(rt_chr_prefix) <- "UCSC"
pairs <- breakpointgr2pairs(rt_chr_prefix)
pairs
```

To see supporting breakpoints clearly, we generate the circos plot
according to the loci of event.

``` r
circos.initializeWithIdeogram(
    data.frame(V1=c("chr13", "chr11"),
               V2=c(21720000,108585000),
               V3=c(21755000,108586000),
               V4=c("q12.11","q24.3"),
               V5=c("gneg","gpos50")))
circos.genomicLink(as.data.frame(S4Vectors::first(pairs)), as.data.frame(S4Vectors::second(pairs)))
```

![](README-unnamed-chunk-8-1.png)<!-- -->

``` r
circos.clear()
```

<!-- # Citation

You can cite `svaNUMT` [here]()

```
@ARTICLE{svaNUMT,
  title    = "",
  author   = "",
  journal  = "",
  volume   = ,
  number   = ,
  pages    = ,
  month    = ,
  year     = ,
  url      = ,
  doi      = ,
  pmc      = 
}
```
-->
