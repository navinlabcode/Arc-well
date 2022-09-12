# Calculation of Breadth of coverage

This [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline calculates the breadth of coverage from a downsampled bam file.

# Installation

Make sure to install and properly run snakemake and clone this repository

# Usage

All `.bam` files must be inside a folder called `data`. after that run

```
snakemake --snakefile /PATH/TO/SNAKEFILE/ds_breadth.smk --cores 20 
```

# Dependencies

- snakemake
- samtools
- sambamba
- bedtools

# Pipeline steps

This pipelines uses `sambamba` to mark and remove duplicates, randomly downsample the bam file to `x` number of reads (default 800k), which can be changed inside the pipeline params and trims the reads to `x` base pairs (default 50) which can also be changed within the snakefile parameters

# Output

A file `.covhist` from bedtools genomeCoverageBed that can be read within R for analysis with the following code

```
library(tidyverse)

gini.index <- function(x, n)
{
  this.order <- order(x)
  x <- x[this.order]
  n <- n[this.order]
  # Fraction of population with given wealth
  f <- n / sum(n)
  # Cumulative total wealth
  s <- cumsum(x*n)
  # Delayed cumulative total wealth
  ds <- c(0, s[-length(s)])
  # Total wealth
  m <- s[length(s)]
  # Formula from wikipedia
  # https://en.wikipedia.org/w/index.php?title=Gini_coefficient&oldid=968643673#Discrete_probability_distribution
  1 - sum(f * (ds + s)) / m
}

calc_coverage <- function(path) {
  
  inpaths <- Sys.glob(paste0(path, "*.covhist.txt"))
  
  coverage.stats <- tibble(bed_path=inpaths) %>%
    mutate(cellname = str_extract(basename(bed_path), "^[^.]*")) %>%
    group_by(cellname) %>%
    summarize(.groups="keep",
              read_tsv(bed_path, 
                       col_names=c("refname", "depth", "count",
                                   "refsize", "frac"),
                       col_types=cols(col_character(), col_double(),
                                      col_double(), col_double(),
                                      col_double())),
    ) %>%
    filter(refname=="genome") %>%
    summarize(breadth = 1 - frac[depth==0],
              gini_index = gini.index(depth, count),
              .groups="keep")
  
}

# breadth and gini coefficient ----

# EXAMPLE_DATA
EXAMPLE_DATA_cov <- calc_coverage(path = "/PATH/TO/EXAMPLE_DATA/covfile/") %>% 
  mutate(sample = "EXAMPLE_DATA",
         tech = "wafer-ffpe")
```


