#!/usr/bin/env Rscript

# currently ignore hets and missing data
create.mhaps <-
  function(GT,
           mhap.loc,
           label.column = "NAME",
           chrom.column = "CHROM",
           pos.column = "POS") {
    samples <-
      colnames(GT)[!(colnames(GT) %in% c(chrom.column, pos.column))]
    GT <- merge(GT, mhap.loc, by = c(chrom.column, pos.column))
    
    GT[GT == "0/0"] <- "W"
    GT[GT == "1/1"] <- "M"
    
    
    mhaps <-
      # for each sample column
      apply(GT[, samples], 2, function(sample)
        
        # for each microhaplotype
        tapply(sample, GT[, label.column], function(genotype)
          
          # merge together the genotypes
          paste0(genotype, collapse = "")))
    
    
    invalid.index <- grepl("[^WM]", mhaps, perl = TRUE)
    mhaps[invalid.index] <- NA
    
    mhaps
  }


prepare.groupings <-
  function(mhaps,
           metadata,
           sample.column,
           group.column) {
    metadata.order <- match(colnames(mhaps), metadata[, sample.column])
    
    metadata[metadata.order, group.column]
  }


calculate.mhap.prop <- function(mhaps, group) {
  mhap.prop.list <-
    # for each microhaplotypes row
    apply(mhaps, 1, function(mhap)
      
      # for each region
      tapply(mhap, group, function(microhaplotype)
        
        # calculate the proportion of microhaplotypes
        as.data.frame(table(microhaplotype) / sum(
          table(microhaplotype)
        ))))
  
  
  # for every microhaplotype
  lapply(mhap.prop.list, function(mhap)
    
    # combine each region data into a table
    do.call(rbind, lapply(names(mhap), function(region)
      
      # with region names added as a new column
      cbind(region, mhap[[region]]))))
}


widen.mhap.prop <- function(mhap.prop) {
  # for every microhaplotype
  lapply(mhap.prop, function(mhap) {
    # convert to wide format
    mhap.prop.wide <-
      reshape(mhap,
              direction = "wide",
              idvar = "region",
              timevar = "microhaplotype")
    
    # remove extra column names
    colnames(mhap.prop.wide) <-
      sub("Freq.", "", colnames(mhap.prop.wide), fixed = TRUE)
    
    mhap.prop.wide
  })
}


combine.mhap.prop <- function(mhap.prop.wide) {
  # calculate how many zeroes needed to pad a table
  mat.ncol <- max(vapply(mhap.prop.wide, length, numeric(1)))
  
  mhap.prop.mat <-
    # combine every microhaplotype
    do.call(rbind, lapply(names(mhap.prop.wide), function(mhap) {
      # select the microhaplotype table
      mhap.wide <- mhap.prop.wide[[mhap]]
      
      # if the table is not wide enough
      last.ncol <- length(mhap.wide)
      if (last.ncol != mat.ncol)
        # add padding to the table
        mhap.wide[(1 + last.ncol):mat.ncol] <- 0
      
      # make column names uniform
      colnames(mhap.wide) <- 1:mat.ncol
      
      # add microhaplotype name to the table
      cbind(mhap, mhap.wide)
    }))
  
  # add proper column names
  colnames(mhap.prop.mat) <-
    c("Amplicon_name", "Region", paste0("Genotype.", 1:(mat.ncol - 1)))
  
  # remove row names
  rownames(mhap.prop.mat) <- NULL
  
  mhap.prop.mat
}


create.mhap.metadata <-
  function(mhap.loc,
           label.column = "NAME",
           chrom.column = "CHROM",
           pos.column = "POS") {
    # create microhaplotype metadata
    mhap.name <- mhap.loc[, label.column]
    mhap.chrom <-
      unique(mhap.loc[, c(label.column, chrom.column)])[, chrom.column]
    
    mhap.start <- tapply(mhap.loc[, pos.column], mhap.name, min)
    mhap.end <- tapply(mhap.loc[, pos.column], mhap.name, max)
    
    # reorder because tapply orders by factor rules
    mhap.start <-
      mhap.start[match(unique(mhap.name), names(mhap.start))]
    mhap.end <- mhap.end[match(unique(mhap.name), names(mhap.end))]
    
    data.frame(
      Amplicon_name = unique(mhap.name),
      Chr = mhap.chrom,
      Start = mhap.start,
      End = mhap.end
    )
  }


calculate.mhap.frequency <-
  function(GT,
           mhap.loc,
           metadata,
           sample.column,
           group.column,
           label.column = "NAME") {
    mhaps <- create.mhaps(GT, mhap.loc, label.column = label.column)
    
    group <-
      prepare.groupings(mhaps, metadata, sample.column, group.column)
    
    mhap.prop <- calculate.mhap.prop(mhaps, group)
    
    mhap.prop.wide <- widen.mhap.prop(mhap.prop)
    
    mhap.prop.mat <- combine.mhap.prop(mhap.prop.wide)
    # fill missing data with zeroes
    mhap.prop.mat[is.na(mhap.prop.mat)] <- 0
    
    mhap.metadata <- create.mhap.metadata(mhap.loc)
    
    # combine microhaplotype metadata and frequencies
    mhap.freq <-
      merge(mhap.metadata,
            mhap.prop.mat,
            by = "Amplicon_name",
            sort = FALSE)
    
    mhap.freq[do.call(order, mhap.freq[, c("Chr", "Start", "End", "Region")]), ]
  }


metadata <-
  read.delim("Pv4_samples.txt",
             check.names = FALSE,
             encoding = "UTF-8")

sample.column <- "Sample"
group.column <- "Population"


# genotype of major alleles
GT <- read.delim("panel.major.GT.txt.gz", check.names = FALSE)
chrom.pos <- strsplit(rownames(GT), ":", fixed = TRUE)
GT["CHROM"] <- vapply(chrom.pos, "[[", character(1), 1)
GT["POS"] <- as.numeric(vapply(chrom.pos, "[[", character(1), 2))


# 3_10 (panel)
mhap.loc <- read.delim("3_10_panel.txt")

mhap.freq <-
  calculate.mhap.frequency(GT,
                           mhap.loc,
                           metadata,
                           sample.column,
                           group.column)
mhap.freq[, "Amplicon_name"] <-
  paste0("3_10_", mhap.freq[, "Amplicon_name"])

mhap.freq.file <- "pjInput_major_3_10_mhap.txt"
write.table(
  mhap.freq,
  file = mhap.freq.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
