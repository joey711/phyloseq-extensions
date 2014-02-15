
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>


# edgeR with phyloseq

Define the extension. Do not need to have anything loaded yet to do this.


```r
################################################################################ 
#' Convert phyloseq OTU count data into DGEList for edgeR package
#' 
#' Further details.
#' 
#' @param physeq (Required).  A \code{\link{phyloseq-class}} or
#'  an \code{\link{otu_table-class}} object. 
#'  The latter is only appropriate if \code{group} argument is also a 
#'  vector or factor with length equal to \code{nsamples(physeq)}.
#'  
#' @param group (Required). A character vector or factor giving the experimental
#'  group/condition for each sample/library. Alternatively, you may provide
#'  the name of a sample variable. This name should be among the output of
#'  \code{sample_variables(physeq)}, in which case
#'  \code{get_variable(physeq, group)} would return either a character vector or factor.
#'  This is passed on to \code{\link[edgeR]{DGEList}},
#'  and you may find further details or examples in its documentation.
#'  
#' @param method (Optional). The label of the edgeR-implemented normalization to use.
#'  See \code{\link[edgeR]{calcNormFactors}} for supported options and details. 
#'  The default option is \code{'RLE'}, which is a scaling factor method 
#'  proposed by Anders and Huber (2010).
#'  At time of writing, the \link[edgeR]{edgeR} package supported 
#'  the following options to the \code{method} argument:
#'  
#'  \code{c('TMM', 'RLE', 'upperquartile', 'none')}.
#'
#' @param ... Additional arguments passed on to \code{\link[edgeR]{DGEList}}
#' 
#' @examples
#' 
phyloseq_to_edgeR = function(physeq, group, method = "RLE", ...) {
    require("edgeR")
    require("phyloseq")
    # Enforce orientation.
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    x = as(otu_table(physeq), "matrix")
    # Add one to protect against overflow, log(0) issues.
    x = x + 1
    # Check `group` argument
    if (identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1) {
        # Assume that group was a sample variable name (must be categorical)
        group = get_variable(physeq, group)
    }
    # Define gene annotations (`genes`) as tax_table
    taxonomy = tax_table(physeq, errorIfNULL = FALSE)
    if (!is.null(taxonomy)) {
        taxonomy = data.frame(as(taxonomy, "matrix"))
    }
    # Now turn into a DGEList
    y = DGEList(counts = x, group = group, genes = taxonomy, remove.zeros = TRUE, 
        ...)
    # Calculate the normalization factors
    z = calcNormFactors(y, method = method)
    # Check for division by zero inside `calcNormFactors`
    if (!all(is.finite(z$samples$norm.factors))) {
        stop("Something wrong with edgeR::calcNormFactors on this data,\n         non-finite $norm.factors, consider changing `method` argument")
    }
    # Estimate dispersions
    return(estimateTagwiseDisp(estimateCommonDisp(z)))
}
################################################################################ 
```


# Citations

If you find this extension or tutorial useful in your work, please cite the following:

### Differential Abundance for Microbiome Data
McMurdie and Holmes (2014) Waste Not, Want Not: Why Rarefying Microbiome Data is Inadmissible.
PLoS Computational Biology *in press*

### phyloseq
McMurdie and Holmes (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.
[PLoS ONE. 8(4):e61217](http://dx.plos.org/10.1371/journal.pone.0061217)

### edgeR
Robinson MD, McCarthy DJ, Smyth GK (2009) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.
[Bioinformatics (Oxford, England) 26: 139–140](http://bioinformatics.oxfordjournals.org/content/26/1/139.short)


# Load example data and try it out


```r
date()
```

```
## [1] "Sat Feb 15 10:54:08 2014"
```

```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.7.17'
```

```r
library("edgeR")
packageVersion("edgeR")
```

```
## [1] '3.4.2'
```



I use here the same example dataset as in the DESeq2 vignette, the publicly available data from a study on colorectal cancer:

[Genomic analysis identifies association of Fusobacterium with colorectal carcinoma](http://genome.cshlp.org/content/22/2/292.long).
Kostic, A. D., Gevers, D., Pedamallu, C. S., Michaud, M., Duke, F., Earl, A. M., et al. (2012). *Genome research*, 22(2), 292-298. 

This work was published ahead of print in [Genome Research](http://genome.cshlp.org/) alongside a highly-related article from a separate group of researchers (long-live reproducible observations!): [Fusobacterium nucleatum infection is prevalent in human colorectal carcinoma](http://genome.cshlp.org/content/22/2/299.long). In case you are interested. For the purposes of example, however, we will stick to the data from the former study, with data available at the [microbio.me/qiime](http://www.microbio.me/qiime/) server.

Study ID:  `1457`

Project Name:  `Kostic_colorectal_cancer_fusobacterium`

Study Abstract:  The tumor microenvironment of colorectal carcinoma is a complex community of genomically altered cancer cells, nonneoplastic cells, and a diverse collection of microorganisms. Each of these components may contribute to carcino genesis; however, the role of the microbiota is the least well understood. We have characterized the composition of the microbiota in colorectal carcinoma using whole genome sequences from nine tumor/normal pairs. Fusobacterium sequences were enriched in carcinomas, confirmed by quantitative PCR and 16S rDNA sequence analysis of 95 carcinoma/normal DNA pairs, while the Bacteroidetes and Firmicutes phyla were depleted in tumors. Fusobacteria were also visualized within colorectal tumors using FISH. These findings reveal alterations in the colorectal cancer microbiota; however, the precise role of Fusobacteria in colorectal carcinoma pathogenesis requires further investigation.



```r
filepath = system.file("extdata", "study_1457_split_library_seqs_and_mapping.zip", 
    package = "phyloseq")
kostic = microbio_me_qiime(filepath)
```

```
## Found biom-format file, now parsing it... 
## Done parsing biom... 
## Importing Sample Metdadata from mapping file...
## Merging the imported objects... 
## Successfully merged, phyloseq-class created. 
##  Returning...
```

```r
# Remove the samples for which the DIAGNOSIS was not included
kosticB = subset_samples(kostic, DIAGNOSIS != "None")
dge = phyloseq_to_edgeR(kosticB, group = "DIAGNOSIS")
# Perform binary test
et = exactTest(dge)
# Extract values from test results
tt = topTags(et, n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
res = tt@.Data[[1]]
alpha = 0.01
sigtab = res[(res$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kosticB)[rownames(sigtab), 
    ], "matrix"))
head(sigtab)
```

```
##         Kingdom         Phylum                 Class             Order
## 156697 Bacteria   Fusobacteria  Fusobacteria (class)   Fusobacteriales
## 211205 Bacteria     Firmicutes            Clostridia     Clostridiales
## 114860 Bacteria     Firmicutes               Bacilli   Lactobacillales
## 49616  Bacteria Proteobacteria   Gammaproteobacteria    Pasteurellales
## 38061  Bacteria   Fusobacteria  Fusobacteria (class)   Fusobacteriales
## 484430 Bacteria Proteobacteria Epsilonproteobacteria Campylobacterales
##                    Family         Genus                 Species  logFC
## 156697   Fusobacteriaceae  Leptotrichia Leptotrichia trevisanii  4.116
## 211205    Lachnospiraceae          <NA>                    <NA> -3.551
## 114860   Streptococcaceae Streptococcus                    <NA>  3.778
## 49616     Pasteurellaceae   Haemophilus  Haemophilus influenzae  3.129
## 38061    Fusobacteriaceae  Leptotrichia                    <NA>  3.423
## 484430 Campylobacteraceae Campylobacter  Campylobacter gracilis  2.759
##        logCPM    PValue       FDR  Kingdom         Phylum
## 156697 11.776 1.491e-25 3.735e-22 Bacteria   Fusobacteria
## 211205 10.528 3.187e-22 3.153e-19 Bacteria     Firmicutes
## 114860 10.975 3.777e-22 3.153e-19 Bacteria     Firmicutes
## 49616  10.586 3.569e-19 2.235e-16 Bacteria Proteobacteria
## 38061  11.062 5.579e-19 2.795e-16 Bacteria   Fusobacteria
## 484430  9.998 8.940e-18 3.732e-15 Bacteria Proteobacteria
##                        Class             Order             Family
## 156697  Fusobacteria (class)   Fusobacteriales   Fusobacteriaceae
## 211205            Clostridia     Clostridiales    Lachnospiraceae
## 114860               Bacilli   Lactobacillales   Streptococcaceae
## 49616    Gammaproteobacteria    Pasteurellales    Pasteurellaceae
## 38061   Fusobacteria (class)   Fusobacteriales   Fusobacteriaceae
## 484430 Epsilonproteobacteria Campylobacterales Campylobacteraceae
##                Genus                 Species
## 156697  Leptotrichia Leptotrichia trevisanii
## 211205          <NA>                    <NA>
## 114860 Streptococcus                    <NA>
## 49616    Haemophilus  Haemophilus influenzae
## 38061   Leptotrichia                    <NA>
## 484430 Campylobacter  Campylobacter gracilis
```


As expected from the original study abstract and title, a *Fusobacterium* OTU was most-significantly differentially abundant between the cancerous and healthy samples.

Here is a bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.


```r
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size = 6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 




---

## Other extensions for the phyloseq package:

#### [DESeq](DESeq.html)

#### [DESeq2](DESeq2.html)

#### [edgeR](edgeR.html)

#### [extensions-index](extensions-index.html)

#### [index](index.html)

