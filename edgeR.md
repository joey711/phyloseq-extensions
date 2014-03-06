
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
## [1] "Wed Mar  5 16:27:53 2014"
```

```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.7.21'
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

This work was published ahead of print in [Genome Research](http://genome.cshlp.org/) alongside a highly-related article from a separate group of researchers (hooray for reproducible observations!): [Fusobacterium nucleatum infection is prevalent in human colorectal carcinoma](http://genome.cshlp.org/content/22/2/299.long). In case you are interested. For the purposes of example, however, we will stick to the data from the former study, with data available at the [microbio.me/qiime](http://www.microbio.me/qiime/) server.

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
```


### Independent filtering

Independent filtering is done automatically in DESeq2, but not in DESeq. Here we will filter OTUs for which the variance across all samples is very low, and we'll do this before ever passing the data to DESeq. 


```r
kosticB
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2505 taxa and 185 samples ]
## sample_data() Sample Data:       [ 185 samples by 71 sample variables ]
## tax_table()   Taxonomy Table:    [ 2505 taxa by 7 taxonomic ranks ]
```

```r
kosticp = transformSampleCounts(kosticB, function(x) {
    x/sum(x)
})
hist(log10(apply(otu_table(kosticp), 1, var)), xlab = "log10(variance)", breaks = 50, 
    main = "A large fraction of OTUs have very low variance")
```

![plot of chunk var-plot-ind-filter](figure/var-plot-ind-filter.png) 

```r
varianceThreshold = 1e-05
keepOTUs = names(which(apply(otu_table(kosticp), 1, var) > varianceThreshold))
kosticB = prune_taxa(keepOTUs, kosticB)
kosticB
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 195 taxa and 185 samples ]
## sample_data() Sample Data:       [ 185 samples by 71 sample variables ]
## tax_table()   Taxonomy Table:    [ 195 taxa by 7 taxonomic ranks ]
```


Here we've used an arbitrary but not-unreasonable variance threshold of 10<sup>-5</sup>. It is important to keep in mind that this filtering is independent of our downstream test. The sample classifications were not used.

Now let's use our newly-defined function to convert the phyloseq data object `kosticB` into an edgeR "DGE" data object, called `dge`.


```r
dge = phyloseq_to_edgeR(kosticB, group = "DIAGNOSIS")
# Perform binary test
et = exactTest(dge)
# Extract values from test results
tt = topTags(et, n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
res = tt@.Data[[1]]
alpha = 0.001
sigtab = res[(res$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kosticB)[rownames(sigtab), 
    ], "matrix"))
dim(sigtab)
```

```
## [1] 62 18
```

```r
head(sigtab)
```

```
##         Kingdom         Phylum                Class            Order
## 156697 Bacteria   Fusobacteria Fusobacteria (class)  Fusobacteriales
## 114860 Bacteria     Firmicutes              Bacilli  Lactobacillales
## 211205 Bacteria     Firmicutes           Clostridia    Clostridiales
## 2394   Bacteria  Bacteroidetes        Flavobacteria Flavobacteriales
## 38061  Bacteria   Fusobacteria Fusobacteria (class)  Fusobacteriales
## 49616  Bacteria Proteobacteria  Gammaproteobacteria   Pasteurellales
##                   Family          Genus                 Species  logFC
## 156697  Fusobacteriaceae   Leptotrichia Leptotrichia trevisanii  4.370
## 114860  Streptococcaceae  Streptococcus                    <NA>  3.830
## 211205   Lachnospiraceae           <NA>                    <NA> -3.599
## 2394   Flavobacteriaceae Capnocytophaga Capnocytophaga ochracea  3.148
## 38061   Fusobacteriaceae   Leptotrichia                    <NA>  3.554
## 49616    Pasteurellaceae    Haemophilus  Haemophilus influenzae  3.191
##        logCPM    PValue       FDR  Kingdom         Phylum
## 156697  12.82 2.931e-27 5.716e-25 Bacteria   Fusobacteria
## 114860  11.95 2.598e-22 2.533e-20 Bacteria     Firmicutes
## 211205  11.52 1.172e-21 7.621e-20 Bacteria     Firmicutes
## 2394    10.86 1.193e-20 5.814e-19 Bacteria  Bacteroidetes
## 38061   12.04 4.827e-20 1.882e-18 Bacteria   Fusobacteria
## 49616   11.55 9.061e-20 2.945e-18 Bacteria Proteobacteria
##                       Class            Order            Family
## 156697 Fusobacteria (class)  Fusobacteriales  Fusobacteriaceae
## 114860              Bacilli  Lactobacillales  Streptococcaceae
## 211205           Clostridia    Clostridiales   Lachnospiraceae
## 2394          Flavobacteria Flavobacteriales Flavobacteriaceae
## 38061  Fusobacteria (class)  Fusobacteriales  Fusobacteriaceae
## 49616   Gammaproteobacteria   Pasteurellales   Pasteurellaceae
##                 Genus                 Species
## 156697   Leptotrichia Leptotrichia trevisanii
## 114860  Streptococcus                    <NA>
## 211205           <NA>                    <NA>
## 2394   Capnocytophaga Capnocytophaga ochracea
## 38061    Leptotrichia                    <NA>
## 49616     Haemophilus  Haemophilus influenzae
```


Here is a bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.


```r
library("ggplot2")
packageVersion("ggplot2")
```

```
## [1] '0.9.3.1'
```

```r
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


As expected from the original study abstract and title, *Fusobacterium* OTUs were most-significantly differentially abundant between the cancerous and healthy samples. If you look closely, two different genera of the Fusobacteria phylum were among the most significantly different, *Leptotrichia* (the winner) as well as *Fusobacterium*.

---

## Paired tests

As mentioned above, the design of this experiment is 95 carcinoma/normal pairs, where each pair comes from the same patient. Although the previous tests are valid, they are conservative in that they do not use this extra information regarding the sample-pairs, and in that sense have forfeited extra power. There is support in edgeR for paired tests, and this is officially described in [one of the edgeR user guides](http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf). It is also demonstrated here in the following.


```r
Diagnosis = get_variable(kosticB, "DIAGNOSIS")
Patient = get_variable(kosticB, "ANONYMIZED_NAME")
# Notice that we have some patients without one of the pairs.
tapply(Patient, Diagnosis, length)
```

```
## Healthy   Tumor 
##      95      90
```

```r
length(levels(Patient))
```

```
## [1] 97
```

```r
any(tapply(Diagnosis, Patient, length) > 2)
```

```
## [1] FALSE
```

```r
sum(tapply(Diagnosis, Patient, length) < 2)
```

```
## [1] 9
```

```r
# Keep only patients with both healthy and cancer samples
keepPatients = names(which(tapply(Diagnosis, Patient, function(x) {
    length(unique(x))
}) == 2))
kosticBpair = subset_samples(kosticB, ANONYMIZED_NAME %in% keepPatients)
Diagnosis = get_variable(kosticBpair, "DIAGNOSIS")
Patient = get_variable(kosticBpair, "ANONYMIZED_NAME")
# With that accounting out of the way, define the design matrix
design = model.matrix(~Patient + Diagnosis)
```


Must estimate the dispersions, as always. This is different than in the function shown above.


```r
# Add one to protect against overflow, log(0) issues.
x = as(otu_table(kosticBpair), "matrix") + 1L
taxonomy = data.frame(as(tax_table(kosticBpair), "matrix"))
# Now turn into a DGEList
x = DGEList(counts = x, group = Diagnosis, genes = taxonomy, remove.zeros = TRUE)
# Calculate the normalization factors and estimate dispersion
x = calcNormFactors(x, method = "RLE")
x = estimateGLMCommonDisp(x, design)
x = estimateGLMTrendedDisp(x, design)
x = estimateGLMTagwiseDisp(x, design)
```


As in the edgeR User's Guide, we proceed to fit a linear model and test for the treatment effect. Note that we can omit the coefficient argument to `glmLRT` because the "treatment effect" (in this case the tissue diagnosis) is the last coeffcient in the model.


```r
fit <- glmFit(x, design)
lrt <- glmLRT(fit)
topTags(lrt)
```

```
## Coefficient:  DiagnosisTumor 
##         Kingdom         Phylum                  Class            Order
## 180285 Bacteria     Firmicutes             Clostridia    Clostridiales
## 72853  Bacteria     Firmicutes             Clostridia    Clostridiales
## 174809 Bacteria  Bacteroidetes            Bacteroidia    Bacteroidales
## 322235 Bacteria  Bacteroidetes            Bacteroidia    Bacteroidales
## 194648 Bacteria     Firmicutes             Clostridia    Clostridiales
## 181056 Bacteria     Firmicutes             Clostridia    Clostridiales
## 469709 Bacteria  Bacteroidetes            Bacteroidia    Bacteroidales
## 186029 Bacteria Actinobacteria Actinobacteria (class) Coriobacteriales
## 157566 Bacteria  Bacteroidetes            Bacteroidia    Bacteroidales
## 248140 Bacteria  Bacteroidetes            Bacteroidia    Bacteroidales
##                   Family            Genus                 Species  logFC
## 180285   Ruminococcaceae Faecalibacterium                    <NA> -1.631
## 72853    Ruminococcaceae Faecalibacterium                    <NA> -1.526
## 174809    Bacteroidaceae      Bacteroides                    <NA> -1.362
## 322235    Bacteroidaceae      Bacteroides   Bacteroides uniformis -1.178
## 194648              <NA>             <NA>                    <NA> -1.149
## 181056   Ruminococcaceae Faecalibacterium                    <NA> -1.401
## 469709    Bacteroidaceae      Bacteroides       Bacteroides dorei -1.309
## 186029 Coriobacteriaceae      Collinsella Collinsella aerofaciens -1.298
## 157566     Rikenellaceae        Alistipes                    <NA> -1.176
## 248140    Bacteroidaceae      Bacteroides      Bacteroides caccae -1.190
##        logCPM    LR    PValue       FDR
## 180285  15.71 92.98 5.278e-22 1.029e-19
## 72853   13.22 73.89 8.255e-18 8.049e-16
## 174809  15.03 54.24 1.772e-13 1.152e-11
## 322235  13.30 50.89 9.786e-13 4.771e-11
## 194648  11.86 49.92 1.602e-12 6.249e-11
## 181056  14.09 48.04 4.182e-12 1.359e-10
## 469709  15.92 43.90 3.450e-11 9.142e-10
## 186029  11.75 43.74 3.750e-11 9.142e-10
## 157566  13.71 43.06 5.311e-11 1.151e-09
## 248140  14.16 41.75 1.035e-10 2.018e-09
```

This test detects OTUs that are differentially abundant in the tumor colon mucosa relative to the healthy colon mucosa (control), adjusting for baseline differences between the patients. This test can be viewed as a generalization of a paired t-test.

Re-make the plot.


```r
respair = topTags(lrt, n = nrow(x), adjust.method = "BH", sort.by = "PValue")
respair = respair@.Data[[1]]
alpha = 0.001
sigtab = respair[(respair$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kosticB)[rownames(sigtab), 
    ], "matrix"))
dim(sigtab)
```

```
## [1] 37 19
```

```r
head(sigtab)
```

```
##         Kingdom        Phylum       Class         Order          Family
## 180285 Bacteria    Firmicutes  Clostridia Clostridiales Ruminococcaceae
## 72853  Bacteria    Firmicutes  Clostridia Clostridiales Ruminococcaceae
## 174809 Bacteria Bacteroidetes Bacteroidia Bacteroidales  Bacteroidaceae
## 322235 Bacteria Bacteroidetes Bacteroidia Bacteroidales  Bacteroidaceae
## 194648 Bacteria    Firmicutes  Clostridia Clostridiales            <NA>
## 181056 Bacteria    Firmicutes  Clostridia Clostridiales Ruminococcaceae
##                   Genus               Species  logFC logCPM    LR
## 180285 Faecalibacterium                  <NA> -1.631  15.71 92.98
## 72853  Faecalibacterium                  <NA> -1.526  13.22 73.89
## 174809      Bacteroides                  <NA> -1.362  15.03 54.24
## 322235      Bacteroides Bacteroides uniformis -1.178  13.30 50.89
## 194648             <NA>                  <NA> -1.149  11.86 49.92
## 181056 Faecalibacterium                  <NA> -1.401  14.09 48.04
##           PValue       FDR  Kingdom        Phylum       Class
## 180285 5.278e-22 1.029e-19 Bacteria    Firmicutes  Clostridia
## 72853  8.255e-18 8.049e-16 Bacteria    Firmicutes  Clostridia
## 174809 1.772e-13 1.152e-11 Bacteria Bacteroidetes Bacteroidia
## 322235 9.786e-13 4.771e-11 Bacteria Bacteroidetes Bacteroidia
## 194648 1.602e-12 6.249e-11 Bacteria    Firmicutes  Clostridia
## 181056 4.182e-12 1.359e-10 Bacteria    Firmicutes  Clostridia
##                Order          Family            Genus
## 180285 Clostridiales Ruminococcaceae Faecalibacterium
## 72853  Clostridiales Ruminococcaceae Faecalibacterium
## 174809 Bacteroidales  Bacteroidaceae      Bacteroides
## 322235 Bacteroidales  Bacteroidaceae      Bacteroides
## 194648 Clostridiales            <NA>             <NA>
## 181056 Clostridiales Ruminococcaceae Faecalibacterium
##                      Species
## 180285                  <NA>
## 72853                   <NA>
## 174809                  <NA>
## 322235 Bacteroides uniformis
## 194648                  <NA>
## 181056                  <NA>
```

```r
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
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
    ggtitle("Log Fold Change of Significant OTUs in a Paired Test")
```

![plot of chunk paired-ggplot](figure/paired-ggplot.png) 


### Other covariates available

As a side note, there are many other interesting patient-sample covariates available in the `sample_data`, which you can access with `sample_data()` and `get_variable()`.


```r
sample_variables(kostic)
```

```
##  [1] "X.SampleID"                    "BarcodeSequence"              
##  [3] "LinkerPrimerSequence"          "NECROSIS_PERCENT"             
##  [5] "TARGET_SUBFRAGMENT"            "ASSIGNED_FROM_GEO"            
##  [7] "EXPERIMENT_CENTER"             "TITLE"                        
##  [9] "RUN_PREFIX"                    "AGE"                          
## [11] "NORMAL_EQUIVALENT_PERCENT"     "FIBROBLAST_AND_VESSEL_PERCENT"
## [13] "DEPTH"                         "TREATMENT"                    
## [15] "AGE_AT_DIAGNOSIS"              "COMMON_NAME"                  
## [17] "HOST_COMMON_NAME"              "BODY_SITE"                    
## [19] "ELEVATION"                     "REPORTS_RECEIVED"             
## [21] "CEA"                           "PCR_PRIMERS"                  
## [23] "COLLECTION_DATE"               "ALTITUDE"                     
## [25] "ENV_BIOME"                     "SEX"                          
## [27] "PLATFORM"                      "RACE"                         
## [29] "BSP_DIAGNOSIS"                 "STUDY_CENTER"                 
## [31] "COUNTRY"                       "CHEMOTHERAPY"                 
## [33] "YEAR_OF_DEATH"                 "ETHNICITY"                    
## [35] "ANONYMIZED_NAME"               "TAXON_ID"                     
## [37] "SAMPLE_CENTER"                 "SAMP_SIZE"                    
## [39] "YEAR_OF_BIRTH"                 "ORIGINAL_DIAGNOSIS"           
## [41] "AGE_UNIT"                      "STUDY_ID"                     
## [43] "EXPERIMENT_DESIGN_DESCRIPTION" "Description_duplicate"        
## [45] "DIAGNOSIS"                     "BODY_HABITAT"                 
## [47] "SEQUENCING_METH"               "RUN_DATE"                     
## [49] "HISTOLOGIC_GRADE"              "LONGITUDE"                    
## [51] "ENV_MATTER"                    "TARGET_GENE"                  
## [53] "ENV_FEATURE"                   "KEY_SEQ"                      
## [55] "BODY_PRODUCT"                  "TUMOR_PERCENT"                
## [57] "LIBRARY_CONSTRUCTION_PROTOCOL" "REGION"                       
## [59] "RUN_CENTER"                    "TUMOR_TYPE"                   
## [61] "BSP_NOTES"                     "RADIATION_THERAPY"            
## [63] "INFLAMMATION_PERCENT"          "HOST_SUBJECT_ID"              
## [65] "PC3"                           "LATITUDE"                     
## [67] "OSH_DIAGNOSIS"                 "STAGE"                        
## [69] "PRIMARY_DISEASE"               "HOST_TAXID"                   
## [71] "Description"
```




---

## Other extensions for the phyloseq package:

#### [DESeq](DESeq.html)

#### [DESeq2](DESeq2.html)

#### [edgeR](edgeR.html)

#### [extensions-index](extensions-index.html)

#### [index](index.html)

