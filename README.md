# Sequencing-depth-bias-on-microdiversity-estimates

Estimation of sequencing depth influence on microdiversity estimates


## Required tools: Enveomics, BBTools, inStrain, Bowtie2, bedtools and samtools programs

1) Enveomics scripts collection can be downloaded from http://enve-omics.ce.gatech.edu/enveomics/download
```
git clone git://github.com/lmrodriguezr/enveomics.git enveomics 
```
2) BBTools scripts collection can be downloaded from https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/

3) inStrain can be downloaded following the Installation instructions from https://instrain.readthedocs.io/en/latest/installation.html

```
git clone https://github.com/MrOlm/instrain.git

cd instrain

pip install .
```
4) Bowtie2 can be downloaded from https://github.com/BenLangmead/bowtie2

5) Beedtools can be downloaded from https://bedtools.readthedocs.io/en/latest/content/installation.html

6) Samtools can be downloaded from https://github.com/samtools/samtools


## Illumina metagenomic trimming

```
bbduk.sh in1=Short.MG1.R1.fastq in2=Short.MG1.R2.fastq out1=Short.MG1.trimmed.R1.fastq out2=Short.MG1.trimmed.R2.fastq ref=adapters.fa ktrim=r k=28 mink=12 hdist=1 tbo=t tpe=t qtrim=rl trimq=20 minlength=100
```

## Convert metagenomes from .fastq to .fasta format

Short-reads:
```
cat Short.MG1.R1.fastq | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}' > Short.MG1.R1.fasta
cat Short.MG1.R2.fastq | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}' > Short.MG1.R2.fasta
```

Long-reads:
```
cat Long.MG1.fastq | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}' > Long.MG1.fasta
```

## Short- and long-read random subsamplings of 1,5,10,20,30,40,50,60,70,80,90% of reads

Short-reads:
```
FastA.subsample.pl -f 1,5,10,20,30,40,50,60,70,80,90 Short.MG1.R1.fasta
```
  a) Selection of same reads for both paired-end files ==> This is an example using the 10% of reads:

```
grep ">" Short.MG1.R1.fasta.10.0000-1.fa > IDS.txt
sed -i 's/-1/-2/g' IDS.txt
FastA.filter.pl IDS.txt Short.MG1.R1.fasta.10.0000-1.fa > Short.MG1.R2.fasta.10.0000-1.fa
```

Long-reads:
```
FastA.subsample.pl -f 1,5,10,20,30,40,50,60,70,80,90 Long.MG1.fasta
```

## Fragment long-reads in 200 bps length 

```
shred.sh in=Long.MG1.fasta.10.0000-1.fa out=Long.Fragmented.MG1.fasta.10.0000-1.fa length=200 minlength=0
```

## Mapping with Bowtie2 (example with long-reads)

```
cat MAGs* > CAT.ALL.MAGs.fasta

bowtie2-build CAT.ALL.MAGs.fasta CAT.ALL.MAGs.db

bowtie2 --reorder --no-unal -f -p 36 -x CAT.ALL.MAGs.db -U Long.Fragmented.MG1.fasta.10.0000-1.fa > Long.Fragmented.MG1.sam
```

## Filtering reads to 95% identity and TAD80 estimation

```
sam.filter.rb -m Long.Fragmented.MG1.sam -o filter95.Long.Fragmented.MG1.sam

samtools view -b filter95.Long.Fragmented.MG1.sam | samtools sort -l 9 -@ 36 -o filter95.Long.Fragmented.MG1.bam

bedtools genomecov -ibam filter95.Long.Fragmented.MG1.bam -bga > filter95.Long.Fragmented.MG1.bam.bg

BedGraph.tad.rb -i filter95.Long.Fragmented.MG1.bam.bg -r 0.8 -s > TAD80.filter95.Long.Fragmented.MG1.txt
```

## Nucleotide diversity quantification using inStrain

```
inStrain profile filter95.Long.Fragmented.MG1.sam MAG1.fasta -o IS_MAG1_Long.Fragmented.MG1 --pairing_filter non_discordant --skip_mm_profiling
```

## Average Nucleotide Identity of reads (ANIr)

```
anir.rb -g MAG1.fasta -m filter95.Long.Fragmented.MG1.sam --m-format sam -a fix --tab ANIr_MAG1_Long.Fragmented.MG1.tsv
```

## Prepare inputs for calculate rarefaction curves of nucleotide diversity and sequencing depth

The inputs for the R script are two tables with the nucleotide diversity or ANIr values and sequencing depth for each species across metagenomes (see example tables)

```
library(R6)

Rarefact <- R6Class("Rarefact", list(
  # Definition
  index  = NULL,
  curves = list(),
  modelf = "log-square",
  dup    = function() {
    x <- Rarefact$new()
    x$index  <- self$index
    x$curves <- self$curves
    x$modelf <- self$modelf
    invisible(x)
  },
  create = function(spp, sites) {
    k <- seq_len(length(spp) * length(sites))
    index <- matrix(k, nrow = length(spp), dimnames = list(spp, sites))
    self$index  = index
    self$curves = list()
    invisible(self)
  },
  
  # Ancilliary Slots
  model.fun = list(),
  
  # Ancillary Getters
  spp    = function() rownames(self$index),
  sites  = function() colnames(self$index),
  curve  = function(sp, site) self$curves[[ self$index[sp, site] ]],
  ranges = function(spp, sites, what = c("logratio", "absolute")) {
    what <- match.arg(what)
    if (missing(spp))   spp   <- self$spp()
    if (missing(sites)) sites <- self$sites()
    
    x <- c()
    y <- c()
    for (i in spp) {
      for (j in sites) {
        xy <- self$curve(i, j)
        if (nrow(xy) == 0) next
        if (what == "logratio") xy[, 2] <- xy[, 2] / tail(xy[, 2], 1)
        x <- c(x, xy[, 1])
        y <- c(y, xy[, 2])
      }
    }
    sel <- x > 0 & y > 0
    x <- x[sel]
    y <- y[sel]
    if(what == "logratio") y <- log2(y)
    list(xlim = range(x), ylim = range(y))
  },
  
  # Ancillary Setters
  add = function(sp, site, seq, div) {
    k <- self$index[sp, site]
    self$curves[[ k ]] <- cbind(seq, div)
    invisible(self)
  },
  
  # Analyses
  filter = function(min.seq) {
    x <- self$dup()
    for (i in seq_along(self$curves)) {
      sel <- self$curves[[i]][, 1] >= min.seq
      x$curves[[i]] <- self$curves[[i]][sel, , drop = FALSE]
    }
    invisible(x)
  },
  model = function(sp, site, warn = TRUE) {
    k <- self$index[sp, site]
    if (is.null(self$model.fun[k][[1]])) {
      if (nrow(self$curves[[k]]) < 3) {
        if (warn) warning("Model cannot be estimated with fewer than 3 points")
        return(NULL)
      } else {
        x <- self$curves[[k]][, 1]
        y <- self$curves[[k]][, 2]
        if (self$modelf == "linear") {
          self$model.fun[[k]] <- lm(y ~ x)
        } else if (self$modelf == "log") {
          self$model.fun[[k]] <- lm(y ~ log(x))
        } else if (self$modelf == "square") {
          self$model.fun[[k]] <- lm(y ~ poly(x, 2))
        } else if (self$modelf == "log-square") {
          self$model.fun[[k]] <- lm(y ~ poly(log(x), 2))
        } else {
          stop("Unknown value for modelf: ", self$modelf)
        }
      }
    }
    invisible(self$model.fun[[k]])
  },
  predict = function(sp, site, seq, ...) {
    m <- self$model(sp, site)
    y <- predict(m, data.frame(x = seq), ...)
    response <- rownames(attr(terms(m), "factors"))[1]
    
    if (response == "log(y)") {
      exp(y)
    } else if (response == "exp(y)") {
      log(y)
    } else {
      y
    }
  },
  
  # Graphics
  plot = function(
    spp, sites, what = c("logratio", "absolute"), plot.opts = list(), ...
  ) {
    what <- match.arg(what)
    if (missing(spp))   spp   <- self$spp()
    if (missing(sites)) sites <- self$sites()
    
    # Plotting canvas
    ranges <- self$ranges(what = what)
    if (is.null(plot.opts[["xlim"]]))
      plot.opts[["xlim"]] <- ranges[["xlim"]]
    if (is.null(plot.opts[["log"]]))
      plot.opts[["log"]] <- ifelse(what == "logratio", "x", "")
    if (is.null(plot.opts[["ylim"]]))
      plot.opts[["ylim"]] <- ranges[["ylim"]]
    if (is.null(plot.opts[["xlab"]]))
      plot.opts[["xlab"]] <- "Sequencing Depth (X)"
    if (is.null(plot.opts[["ylab"]]))
      plot.opts[["ylab"]] <- ifelse(
        what == "logratio", "Log2(Diversity Ratio)", "Diversity"
      )
    if (is.null(plot.opts[["main"]]))
      plot.opts[["main"]] <- paste0(
        ifelse(what == "logratio", "Relative ", ""), "Diversity Rarefaction"
      )
    plot.opts[["type"]] <- "n"
    plot.opts[["x"]] <- 0
    do.call("plot", plot.opts)
    if(what == "logratio") abline(h = 0, lty = 3, col = "gray")
    
    # Data curves
    for (i in spp) {
      for (j in sites) {
        d <- self$curve(i, j)
        if (nrow(d) <= 3) next
        r <- cor(d)[1, 2]
        col <- rgb((1 - r) / 2, 0.5, 0.5, 0.75)
        y <- d[, 2]
        if (what == "logratio") y <- log2(y / tail(y, n = 1))
        lines(d[, 1], y, col = col, ...)
      }
    }
  },
  plot.error = function(spp, sites, modelf, plot.opts = list(), ...) {
    if (missing(spp))   spp   <- self$spp()
    if (missing(sites)) sites <- self$sites()
    if (!missing(modelf)) {
      self$model.fun <- list()
      self$modelf <- modelf
      if (is.null(plot.opts[["main"]])) plot.opts[["main"]] <- self$modelf
    }
    
    plot.opts[["x"]] <- 0
    plot.opts[["t"]] <- "n"
    if (is.null(plot.opts[["xlim"]])) plot.opts[["xlim"]] <- c(0, 1)
    if (is.null(plot.opts[["ylim"]])) plot.opts[["ylim"]] <- c(-0.2, 1.2)
    if (is.null(plot.opts[["xaxs"]])) plot.opts[["xaxs"]] <- "i"
    if (is.null(plot.opts[["yaxs"]])) plot.opts[["yaxs"]] <- "i"
    if (is.null(plot.opts[["xlab"]])) plot.opts[["xlab"]] <- "observed (minmax)"
    if (is.null(plot.opts[["ylab"]]))
      plot.opts[["ylab"]] <- "predicted (minmax of observed)"
    do.call(plot, plot.opts)
    
    for (i in spp) {
      for (j in sites) {
        d <- self$curve(i, j)
        if (nrow(d) < 3) next
        range <- diff(range(d[, 2]))
        lines(
          (d[, 2] - min(d[, 2])) / range,
          (self$predict(i, j, d[, 1]) - min(d[, 2])) / range,
          ...
        )
      }
    }
    abline(0, 1, col = 2, lwd = 3)
  }
))

read.rarefact <- function(file.seq, file.div) {
  seq <- read.table(file.seq, header = TRUE)
  div <- read.table(file.div, header = TRUE)
  stopifnot(dim(seq) == dim(div), seq[, 1:2] == div[, 1:2])
  stopifnot(colnames(seq) == colnames(div))
  
  x <- Rarefact$new()$create(unique(div[, 2]), colnames(div)[-2:-1])
  for (i in x$spp()) {
    for (j in x$sites()) {
      x <- x$add(i, j, seq[seq[, 2] == i, j], div[div[, 2] == i, j])
    }
  }
  return(x)
}

# Example use:
# Load data using `read.rarefact`
rare <- read.rarefact("Seqdeph.Table.csv", "inStrain.Table.csv")

# Plot all the data as is
rare$plot(lwd = 2, plot.opts = list(las = 1))

# Filter by depth ≥ 10X and plot
rare$filter(10)$plot()
```

## Interpolation of nucleotide diversity to a standardized coverage (e.g. 200X)

```
# Plot error of different models
rare.f <- rare$filter(10)
layout(matrix(3:3, nrow = 3))
rare.f$plot.error(modelf = "linear")
rare.f$plot.error(modelf = "log")
rare.f$plot.error(modelf = "square")
rare.f$plot.error(modelf = "log-square")

# Reset `rare.f` to make sure the `modelf` value is as expected
rare.f <- rare$filter(10)
rare.f$modelf <- "log-square"

# Predict the diversity of the first species in the first site at sequencing depth of 200X

spp <- rare.f$spp()
sites <- rare.f$sites()

predicted_values <- matrix(NA, nrow = length(spp), ncol = length(sites), dimnames = list(spp, sites))

for (i in seq_along(spp)) {
  for (j in seq_along(sites)) {
    # Verificar si hay suficientes puntos para ajustar el modelo
    if (nrow(rare.f$curves[[rare.f$index[spp[i], sites[j]]]]) >= 3) {
      # Intentar obtener el valor predicho para la especie i en la muestra j
      tryCatch({
        predicted_values[i, j] <- rare.f$predict(spp[i], sites[j], 200)
      }, error = function(e) {
        # Capturar el error y continuar con la siguiente iteración
        next
      })
    }
  }
}

print(predicted_values)
```

## Example of rarefaction curve graph:

![figure](Example.Rarefaction.pdf)
