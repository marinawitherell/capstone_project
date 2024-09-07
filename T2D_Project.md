T2D_Project
================
2024-09-07

## Loading the packages

``` r
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

## Upload the data

``` r
mydata<- readRDS("raw_data/T2D_1000.RDS")
mydata
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 12062 taxa and 1000 samples ]
    ## sample_data() Sample Data:       [ 1000 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]

## Normalization (saveRDS)

``` r
normalized_data <- transform_sample_counts(mydata, function(x) x / sum(x))
```

## Ordination Plots

``` r
## Taxa Ordination
GP.ord <- ordinate(normalized_data, "PCoA")
plot_ordination(normalized_data, GP.ord, type="taxa", color="Phylum", title="Taxa in Normalized Data")
```

![](T2D_Project_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
## Full Ordination
GP.ord <- ordinate(normalized_data, "PCoA")
plot_ordination(normalized_data, GP.ord, title="sample", color="sample_body_site")
```

![](T2D_Project_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
## New Ordination
new_ord <- ordinate(
  physeq = normalized_data,
  method = "CAP",
  distance = "bray",
  formula = ~ sample_body_site
)

plot_ordination(normalized_data, new_ord, color="sample_body_site", title="Samples in Normalized Data")
```

![](T2D_Project_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

## Permanova

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-4

``` r
#calculate distance
myphyloseq_bray <- phyloseq::distance(normalized_data, method="bray")

sampledf <- data.frame(sample_data(normalized_data))

# adonis(myphyloseq_bray~sample_body_site, data=sampledf)

beta <- betadisper(myphyloseq_bray, sampledf$sample_body_site)
permutest(beta)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
    ## Groups      1  0.7971 0.79709 50.727    999  0.001 ***
    ## Residuals 998 15.6820 0.01571                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Hypothesis Testing

How different is the diversity between skin and stool samples, or what
is their beta diversity? If it is found significant, what are the most
prevalent taxa present in each type of sample?

``` r
## The above output has a p-value of 0.001. Since this is less than 0.05, we can reject the null hypothesis that the two sample body sites (skin and stool) have the same taxa.

library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
## Finding the top 10 taxa of each site
feces_phy <- subset_samples(normalized_data, sample_body_site=="feces")
nasal_phy <- subset_samples(normalized_data, sample_body_site == "nasal cavity")

taxa_skin_df = as.data.frame(tax_table(nasal_phy))
taxa_stool_df = as.data.frame(tax_table(feces_phy))

## find the row sums of otu
skin_rows = rowSums(otu_table(nasal_phy))
stool_rows = rowSums(otu_table(feces_phy))

## add row sums as a column to the taxa table
skin_rows <- as.data.frame(skin_rows)
stool_rows <- as.data.frame(stool_rows)
taxa_skin_df = merge(taxa_skin_df, skin_rows, by=0)
taxa_stool_df = merge(taxa_stool_df, stool_rows, by=0)

## format
row.names(taxa_skin_df) = taxa_skin_df[,1]
taxa_skin_df = taxa_skin_df[,-1]
row.names(taxa_stool_df) = taxa_stool_df[,1]
taxa_stool_df = taxa_stool_df[,-1]

## organize taxa by the abundance 
taxa_skin_df = taxa_skin_df[order(taxa_skin_df$skin_rows, decreasing = TRUE), ]
taxa_stool_df = taxa_stool_df[order(taxa_stool_df$stool_rows, decreasing = TRUE), ]

## get otu tables
otu_skin = as.data.frame(otu_table(nasal_phy))
otu_stool = as.data.frame(otu_table(feces_phy))

## function to filter rows
relevent_rows = vector()
filter_func <- function(df){
    for(i in 1:nrow(df)){
      row <- df[i,]
      row_name <- row.names(df)[i]
      count = 0
      for (x in row){
        if (x > 0){
          count <- count+1
        }
      }
      if (count > 3){
        relevent_rows <- c(relevent_rows, row_name)
      }
    }
  return(relevent_rows)
}


## function to filter columns
relevent_cols = vector()
filter_func2 <- function(df){
    for(i in 1:ncol(df)){
      col <- df[,i]
      col_name <- colnames(df)[i]
      if (colSums(df)[i] > 0){
        relevent_cols <- c(relevent_cols, col_name)
      }
    }
  return(relevent_cols)
}

## filter rows of each sample body site
rel_rows_skin <- filter_func(otu_skin)
rel_rows_stool <- filter_func(otu_stool)

## get filtered taxa table
filtered_skin_taxa =  taxa_skin_df[rownames(taxa_skin_df) %in% rel_rows_skin, ]
filtered_stool_taxa =  taxa_stool_df[rownames(taxa_stool_df) %in% rel_rows_stool, ]

## top 10 taxa
skin_10 = filtered_skin_taxa[1:10, ]
stool_10 = filtered_stool_taxa[1:10, ]

## filter top 10 into phyloseq object
my_phy1 <- subset_taxa(normalized_data, 
                        rownames(tax_table(normalized_data)) %in% rownames(skin_10) | rownames(tax_table(normalized_data)) %in% rownames(stool_10))

## get otu table
phy_otu1 = as.data.frame(otu_table(my_phy1))

## find columns with sums > 0
rel_col1 = filter_func2(phy_otu1)

## filter the phyloseq object
my_phy1 <- subset_samples(my_phy1, rownames(sample_data(my_phy1)) %in% rel_col1)


## top 60 taxa
skin_60 = filtered_skin_taxa[1:60, ]
stool_60 = filtered_stool_taxa[1:60, ]

## filter top 60 into phyloseq object
my_phy2 <- subset_taxa(normalized_data, 
                        rownames(tax_table(normalized_data)) %in% rownames(skin_60) | rownames(tax_table(normalized_data)) %in% rownames(stool_60))

## get otu table
phy_otu2 = as.data.frame(otu_table(my_phy2))

## find columns with sums > 0
rel_col2 = filter_func2(phy_otu2)

## filter the phyloseq object
my_phy2 <- subset_samples(my_phy2, rownames(sample_data(my_phy2)) %in% rel_col2)
```

``` r
## plotting the top 10 otus
GP.ord <- ordinate(my_phy1, "PCoA")

p1 = plot_ordination(my_phy1, GP.ord, type="taxa", color="Order", title="A") +
guides(shape = guide_legend(override.aes = list(size = 2)),
               color = guide_legend(override.aes = list(size = 2))) +
        theme(legend.title = element_text(size = 7), 
              legend.text  = element_text(size = 7),
              legend.key.size = unit(0.5, "lines"),
              legend.box.spacing = margin(0.5),
              axis.title=element_text(size=10),
              plot.title = element_text(face="bold")) +
  geom_point(size = 4)

p2 = plot_ordination(my_phy1, GP.ord, type="Sample", color="sample_body_site", title = "B")+
guides(shape = guide_legend(override.aes = list(size = 2)),
               color = guide_legend(override.aes = list(size = 2))) +
        theme(legend.title = element_text(size = 7), 
              legend.text  = element_text(size = 7),
              legend.key.size = unit(0.5, "lines"),
              legend.box.spacing = margin(0.5),
              axis.title=element_text(size=10),
              plot.title = element_text(face="bold"))

## top 60 otus
GP.ord <- ordinate(my_phy2, "PCoA")
p3 = plot_ordination(my_phy2, GP.ord, type="taxa", color="Order", title="C")+
  guides(shape = guide_legend(override.aes = list(size = 2)),
                 color = guide_legend(override.aes = list(size = 2))) +
          theme(legend.title = element_text(size = 7), 
                legend.text  = element_text(size = 7),
                legend.key.size = unit(0.5, "lines"),
                legend.box.spacing = margin(0.5),
                axis.title=element_text(size=10),
                plot.title = element_text(face="bold")) +
  geom_point(size = 3)

p4 = plot_ordination(my_phy2, GP.ord, type="Sample", color="sample_body_site", title = "D")+
guides(shape = guide_legend(override.aes = list(size = 2)),
               color = guide_legend(override.aes = list(size = 2))) +
        theme(legend.title = element_text(size = 7), 
              legend.text  = element_text(size = 7),
              legend.key.size = unit(0.5, "lines"),
              legend.box.spacing = margin(0.5),
              axis.title=element_text(size=10),
              plot.title = element_text(face="bold"))

grid.arrange(arrangeGrob(p1, top = "Taxa", left="Top 20 OTUs"), arrangeGrob(p2, top = "Samples"), arrangeGrob(p3, left="Top 60 OTUs"), p4,
                         widths=c(6,5), heights = c(2,2))
```

![](T2D_Project_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
library(vegan)

#calculate distance
myphyloseq_bray <- phyloseq::distance(my_phy1, method="bray")

sampledf <- data.frame(sample_data(my_phy1))

#adonis(myphyloseq_bray~sample_body_site, data=sampledf)

beta <- betadisper(myphyloseq_bray, sampledf$sample_body_site)
permutest(beta)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df Sum Sq  Mean Sq     F N.Perm Pr(>F)  
    ## Groups      1  0.143 0.143020 5.635    999  0.024 *
    ## Residuals 997 25.305 0.025381                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#calculate distance
myphyloseq_bray <- phyloseq::distance(my_phy2, method="bray")

sampledf <- data.frame(sample_data(my_phy2))

#adonis(myphyloseq_bray~sample_body_site, data=sampledf)

beta <- betadisper(myphyloseq_bray, sampledf$sample_body_site)
permutest(beta)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df Sum Sq  Mean Sq      F N.Perm Pr(>F)   
    ## Groups      1  0.272 0.272039 15.031    999  0.002 **
    ## Residuals 998 18.062 0.018098                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
