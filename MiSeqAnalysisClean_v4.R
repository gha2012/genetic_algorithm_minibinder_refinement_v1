finish <- FALSE
while (finish == F) {
########### Install all necessary packages #############
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 
   pkgs <- c("readxl", "tidyverse", "ShortRead", "ggplot2", "ggbreak",
             "plyr", "grid", "svDialogs", "shiny", "openxlsx")
   
  for(i in 1:length(pkgs)){
    if(!require(pkgs[i], character.only = T)){
      install.packages(pkgs[i])
      require(pkgs[i], character.only = T)
    }else{
      require(pkgs[i], character.only = T)
    }
  }


########### Select all necessary files ##############
   
   ## Select working directory
   dlg_message("Please set your working directory")
   wd <- dlg_dir(
     title = "Select your working directory",
     multiple = FALSE,
     filters = dlg_filters["All", ],
     gui = .GUI
   )$res
   setwd(wd)
   
  
  ## Name experiment
  ExpName <- dlgInput("Enter the experiment name", "e.g. Sort xy")$res
   
   
  ## Open sequencing file(s)
  dlg_message("Please choose the sequencing file you want to analyze.")
  file <- dlg_open(
    default = getwd(),
    title = "Select your sequencing file",
    multiple = FALSE,
    filters = dlg_filters["All", ],
    gui = .GUI
  )$res
  filename <- dlgInput("Enter the sample name", "Sample name")$res
  
  # data frame with all NGS files
  files <- NULL
  filenames <- NULL
  files[1] <- file
  filenames[1] <- filename
  
  i <- 1
  x <- T
  while (x == T) {
    res <- dlg_message("Do you want to analyze an additional sequencing file?
                       It has to be of the same type (e.g. both genomic DNA)!", "yesno")$res
    if (res == "yes") {
      i = i + 1
      file <- dlg_open(
        default = getwd(),
        title = "Select your sequencing file",
        multiple = FALSE,
        filters = dlg_filters["All", ],
        gui = .GUI
      )$res
      filename <- dlgInput("Enter the sample name", "Sample name")$res
      files[i] <- file
      filenames[i] <- filename
    } else {
      x = F
    }
  }
    
  
  ## primer sequence upstream of the binders
  primer <- dlgInput("Please enter the binding sequence of the forward primer used for the first round of amplification",
                     "ATCG")$res
  
  gap <- dlgInput("Please enter the number of nucleotides between the primer and the binder sequence",
                     "e.g. 4")$res
  
  # the primer sequence will be cut off upstream of the binder
  # gap is the length of shared sequence between the primer and the binder that has to subtracted from the analyzed sequence
  
  
  
  ## Open reference file
  dlg_message("Please select the reference file you want to match your sequencing results to.")
  refF <- read_excel(
    dlg_open(
      default = getwd(),
      title = "Select the reference file",
      multiple = FALSE,
      filters = dlg_filters["All", ],
      gui = .GUI
    )$res
  )
  
  # select column of reference file which holds the target names
  dlg_message("Please use the console/ pop-up window to select the column of your reference file that specifies the binders' target.")
  target <- select.list(names(refF), title = "Select the column that specifies the targets", graphics = TRUE)
  
  # select column of reference file which holds the binder sequences that shall be matched with NGS sequences
  dlg_message("Please use the console/ pop-up window to select the column of your reference file that specifies the binders' names.")
  binders <- select.list(names(refF), title = "Select the column that holds the binders' names", graphics = TRUE)
  
  # select column of reference file which holds the binder sequences that shall be matched with NGS sequences
  dlg_message("Please use the console/ pop-up window to select the column of your reference file that specifies the binders' DNA sequences.")
  dna_seqs <- select.list(names(refF), title = "Select the column that holds the binders' DNA sequences", graphics = TRUE)
  
  
  
  ## Define output folder to save result files
  dlg_message("Please select the output folder where you wish to save the results.")
  output <- dlg_dir(
    title = "Select output folder",
    multiple = FALSE,
    filters = dlg_filters["All", ],
    gui = .GUI
  )$res
  
  

########### Actual analysis ###############
  dlg_message("The data analysis might take a while, please confirm to start.")
  
  ## Identify the sequence length that uniquely identifies each binder of this reference file
  matches <- NULL
  matches <- data.frame(refF[, c(1, 2)])
  matches$count <- NA
  k = 0
  z = 0
  sub = 20
  while (all(sapply(matches$count, function(x)
    is.numeric(x) & all(x) %in% 0)) == FALSE) {
    sub = sub + 1
    ref <- substr(refF[[dna_seqs]], 1, sub)
    matches$count <- 0
    k = 0
    z = 0
    for (i in 1:length(ref)) {
      for (j in 1:length(ref)) {
        if (i != j) {
          if (ref[i] == ref[j]) {
            matches[i, 3] <- matches[i, 3] + 1
          }
          else{
            k = k + 1
          }
        }
        else{
          z = z + 1
        }
      }
    }
    print(sub)
  }
  
  
  ## NGS data quality filter
  numFiles <- length(files)
  allSeqs <- data.frame() # stores all sequences that passed the quality filter Q>30
  allStats <- data.frame() # stores statistics of NGS analyses of all files

  for (k in 1:numFiles) {
    fastQ <- readFastq(files[k])
    cycles <- width(fastQ)
    
    seqs <- sub(paste0(".*", primer), "", sread(fastQ))   # remove sequence of primer upstream of binders
    seqs <- substr(seqs, (as.numeric(gap) + 1), cycles)   # remove remaining sequence upstream of binders
    # seq holds all NGS sequences reduced to potential binder sequences
    
    qual <- quality(fastQ)
    rel <- alphabetScore(qual) / cycles
    seqQual <- bind_cols(seqs, rel)                 # all sequences with quality score
    colnames(seqQual) <- c("Sequences", "QScores")
    seqClean <- seqQual[seqQual$QScores >= 30, ]    # only sequences with quality score >30
    
    allStats[1,k] <- nrow(seqQual)  # number of reads
    allStats[2,k] <- nrow(seqClean)  # number of reads with quality score >30
    allStats[3,k] <- round(nrow(seqClean) / nrow(seqQual) * 100, 2)  # percentage of sequences with quality score >30
    colnames(allStats)[k] <- filenames[k]
    
    allSeqs <- merge(data.frame(allSeqs, row.names=NULL), data.frame(seqClean$Sequences, row.names=NULL), 
          by = 0, all = TRUE)[-1]
    colnames(allSeqs)[k] <- filenames[k]
  }
  rownames(allStats)[1] <- "Total number of reads"
  rownames(allStats)[2] <- "Number of reads above quality treshold"
  rownames(allStats)[3] <- "Percentage of reads that passed the quality treshold"
  
  matches <- NULL
  matches <- data.frame(refF[c(target, binders)])
  ## Match sequencing results with library sequences
  for (k in 1:numFiles) {
    p <- k + 2
    matches[, p] <- 0
    
    if (allStats[2, k] > 0) {
      # General statistics
      ref <- substr(refF[[dna_seqs]], 1, sub)  # take reference sequences and reduce length to what's necessary to uniquely identify a binder
      MiSeq <- substr(na.exclude(allSeqs[, k]), 1, sub)  # take NGS sequences and reduce length to what's necessary to uniquely identify a binder
      
      check <- match(MiSeq, ref)
      allStats[4, k] <- length(na.exclude(unique(check))) # number binders identified (= number of unique sequences)
      
      binderTot <- as.numeric(nrow(matches)) # total number of binders from reference file
      allStats[5, k] <- format(round(as.numeric(allStats[4, k]) / binderTot * 100, 2), nsmall = 2) # Percentage of identified binders from total binders in reference file
      allStats[6, k] <- format(round(length(na.exclude(check)) / nrow(allSeqs[k]) * 100, 2), nsmall = 2) # Percentage of matched clean (Q>30) sequences
      
      # Count for each sequence
      for (j in 1:length(MiSeq)) {
        for (i in 1:length(ref)) {
          if (ref[i] == MiSeq[j]) {
            matches[i, p] <- matches[i, p] + 1
            i = length(ref)                         # stop inner loop if sequence is identified
          }
        }
      }
      colnames(matches)[p] <- filenames[k]
      
      minx <- matches[matches[, p] > 0, ] # select only binders with counts > 0
      min <- round_any(min(minx[, p]), 1, f = ceiling) # lowest count (> 0)
      minZ <- nrow(subset(matches, matches[, p] %in% c(min(minx[, p])))) # number of binders with lowest this count
      max <- round_any(max(matches[, p]), 1, f = ceiling) # highest count
      
      allStats[7, k] <- min # lowest count (> 0)
      allStats[8, k] <- minZ # number of binders with lowest this count
      allStats[9, k] <- max # highest count
      
      print(k)
    }
    else{
      allStats[4, k] <- 0 # number binder identified (= number of unique sequences)
      allStats[5, k] <- 0
      allStats[6, k] <- 0
      
      p <- k + 2
      colnames(matches)[p] <- filenames[k]
      
      min <- 0
      minZ <- 0
      max <- 0
      
      allStats[7, k] <- min # lowest count (> 0)
      allStats[8, k] <- minZ # number of binders with lowest this count
      allStats[9, k] <- max # highest count
      
      print(k)
    }
  }
  
  rownames(allStats)[4] <- "Number of identified binder"
  rownames(allStats)[5] <- "Percentage of identified binder of total reference binder"
  rownames(allStats)[6] <- "Percentage of matched clean sequences"
  rownames(allStats)[7] <- "Lowest count of reads"
  rownames(allStats)[8] <- "Number of binders with the lowest count"
  rownames(allStats)[9] <- "Highest count of reads"
  
  
  
########### Create output #############
  
  ## save results as .csv
  if (file.exists(file.path(output, ExpName))){
    setwd(file.path(output, ExpName))
  } else {
    dir.create(file.path(output, ExpName))
    setwd(file.path(output, ExpName))
  }
  write.table(matches, paste(ExpName, sep = ""), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  # table format of read counts to use for https://jsb-lab.bio/xyplot/
  write.csv(allStats, paste(ExpName, ".csv", sep = ""))
  # saves important analysis statistics
  
  ## plot and save results as boxplots with individual data points
  for (k in 1:numFiles) {
    p <- k + 2
    max <- round_any(max(matches[,p]), 100, f = ceiling)
    
    ggplot() +
      geom_boxplot(aes(x = matches[,1], y = matches[,p]), position = "dodge", outlier.shape = NA) +
      geom_jitter(aes(x = matches[,1], y = matches[,p])) +  # Jitter for all points
      labs(
        x = "Binder group",
        y = "Count per binder sequence",
        title = "Sequence representation",
        caption = paste(
          "\n",
          allStats[2,k], " reads were analyzed of which\n",
          allStats[6, k], "% could be matched to a binder,\n",
          "representing ", allStats[4,k], " out of ", binderTot, " total binders (", allStats[5,k], "%).",
          sep = ""
        )
      ) +
      theme_classic() +
      theme(plot.caption = element_text(size = 10L, hjust = 0)) +
      scale_y_continuous(
        breaks = c(seq(0, max, (max/5))),
        limits = c(0, max),
        expand = c(0, 5)
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
      )
    
    ggsave(paste("Binder-representation_", filenames[k], "_indvPoints_QClean.png", sep = ""),
           width = 7, height = 6, dpi = 300,
           plot = last_plot())
    
    ## plot and save results as boxplots with boxplots only
    ggplot() +
      geom_boxplot(aes(x = matches[,1], y = matches[,p]), position = "dodge", outlier.shape = T) +
      labs(
        x = "Binder group",
        y = "Count per binder sequence",
        title = "Sequence representation",
        caption = paste(
          "\n",
          allStats[2,k], " reads were analyzed of which\n",
          allStats[6, k], "% could be matched to a binder,\n",
          "representing ", allStats[4,k], " out of ", binderTot, " total binders (", allStats[5,k], "%).",
          sep = ""
        )
      ) +
      theme_classic() +
      theme(plot.caption = element_text(size = 10L, hjust = 0)) +
      scale_y_continuous(
        breaks = c(seq(0, max, (max/5))),
        limits = c(0, max),
        expand = c(0, 5)
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
      )
    
    ggsave(paste("Binder-representation_", filenames[k], "_BoxPlot_QClean.png", sep = ""),
           width = 7, height = 6, dpi = 300,
           plot = last_plot())
  }
  
  dlg_message("The analysis is finished. Please find the results in the defined output folder.")
  
  ##
  finish <- TRUE
}
