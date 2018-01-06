Phyloseqify the Fratjuice Samples
================
Richard Sprague
1/6/2018

Load the prerequisites. This code block could be eliminated if these were included in an environment file.

``` r
#library(phyloseq)
#library(actino)
#library(tidyverse)

# ugly hack to avoid modifying existing project structure
SAMPLEPATH <- "../samples/ubiome/ubiome-processed-outputs"
```

Read the mapping data into a more familiar "Actino-style" mapfile. \[See the crappy documentation for the [Actino package](https://github.com/richardsprague/actino)\].

``` r
results <- readxl::read_excel(file.path(SAMPLEPATH, "/results.xlsx"))
kits <- read.delim(file.path(SAMPLEPATH,"../kits.tsv"))
samples <- read.delim(file.path(SAMPLEPATH,"../../samples.tsv"))


# not quite working: want to add other 'samples' metadata (eg. "Geo" field) to the mapfile
# results[match(results$SeqID,kits$seq_id),]
# paste(samples$fraternity,samples$location)
# match(kits$sample_id,samples$sample_id)]

mapping_df <- data.frame(SSR=results$SeqID, id=results$Sample_ID, Site = kits$tube)
```

Save it as an Excel file to be later readable by the actino functions. You only have to run this chunk once. It's commented now because I already created the map file.

``` r
#library(xlsx)
#xlsx::write.xlsx(mapping_df, file = file.path(SAMPLEPATH,"/json/frat_map_file.xlsx"), row.names = FALSE)
```

Now create the Phyloseq object

``` r
sample_files<-actino::just_json_files_in(file.path(SAMPLEPATH,"/json"))

frat_map_file <- paste0(SAMPLEPATH,"/json/frat_map_file.xlsx")
fratjuice.phylum <- actino::phyloseq_from_JSON_at_rank(sample_files, frat_map_file, rank = "phylum", count.normalized = TRUE)
fratjuice.genus <- actino::phyloseq_from_JSON_at_rank(sample_files, frat_map_file, rank = "genus", count.normalized = TRUE)
phyloseq::sample_data(fratjuice.genus)$id <- mapping_df$id
fratjuice.species <- actino::phyloseq_from_JSON_at_rank(sample_files, frat_map_file, rank = "species", count.normalized = TRUE)
phyloseq::sample_data(fratjuice.genus)$id <- mapping_df$id

phyloseq::sample_data(fratjuice.phylum)
```

    ##           SSR  Site
    ## 300909 300909   Gut
    ## 300921 300921 Spare
    ## 300945 300945   Gut
    ## 300948 300948 Spare
    ## 300960 300960   Gut
    ## 300975 300975 Spare
    ## 300987 300987   Gut
    ## 301026 301026 Spare
    ## 301029 301029   Gut
    ## 301032 301032 Spare

plot the Phylum level

``` r
phyloseq::sample_data(fratjuice.phylum)$id <- mapping_df$id
phyloseq::plot_bar(fratjuice.phylum, fill = "Phylum", x = "id") + 
  ggplot2::scale_y_continuous(labels=function(x)x/10000) + ggplot2::ylab("Abundance (%)")
```

![](frat-phyloseq_files/figure-markdown_github/unnamed-chunk-5-1.png)

There are 537 unique genus-level taxa found in these samples, so showing all of them would take up too much space. Let's just plot the most common ones:

``` r
frat.topgenus <- phyloseq::prune_taxa(phyloseq::taxa_sums(fratjuice.genus)>50000, fratjuice.genus)
phyloseq::plot_bar(frat.topgenus, fill = "Genus", x = "id") + 
  ggplot2::scale_y_continuous(labels=function(x)x/10000) + ggplot2::ylab("Abundance (%)")
```

![](frat-phyloseq_files/figure-markdown_github/unnamed-chunk-6-1.png)
