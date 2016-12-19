ReadMrd <- function(file) {
  # Loads mrd files into a mirDeep object
  #
  # Args:
  #   file: The path/file/file object to read
  # Returns:
  #   mirDeep object
  
  # Read the file into a character vector
  f.handle <- file(file, "r")
  lines <- readLines(f.handle)
  close(f.handle)
  
  # split the records into a list
  record.pos <- which(substr(lines, 1, 1) == ">")
  record.pos <- c(record.pos, length(lines))
  record.list <- list()
  last <- NA
  for (pos in record.pos) {
    if (is.na(last)) {
      last <- pos
      next
    }
    
    name <- substring(lines[last],2)
    record.list[[name]] <- lines[(last+1):(pos-1)]
    last <- pos
  }
  
  #parse records
  record.list <- lapply(record.list, ParseMrdRecord)
  record.list
}
ParseMrdRecord <- function(recordLines) {
  # parses mrd output record into an mrd object
  
  out.list <- list()
  
  # split on whitespace
  record <- strsplit(recordLines, "\\s+")
  
  # parse header
  
  if(record[[1]][1] == "total") { 
    # mirbase style records
    # here is where mirDeep gets dicey.  If there is an exact match of multiple
    # miRNAs in a single primary sequence we will get a line for all of them. 
    # This means there will be multiple and unpredictable lines of reads. 
    i <- 1
    while(length(record[[i]]) > 2) {
      name <- record[[i]][1]
      value <- record[[i]][4]
      if(is.na(value)) value <- 0
      out.list[[name]] <- value
      i <- i + 1
    }
    begin <- i
  } else {
    # mirdeep discovery style records
    out.list$score_total <- as.integer(record[[1]][3])
    out.list$score_star <- as.integer(record[[2]][5])
    out.list$score_mature <- as.integer(record[[3]][5])
    out.list$score_mfe <- as.integer(record[[4]][4])
    out.list$score_randfold <- as.integer(record[[5]][4])
    out.list$score_seed <- as.integer(record[[6]][5])
    out.list$miRNA_family <- record[[7]][5]
    out.list$total <- as.integer(record[[8]][4])
    out.list$mature <- as.integer(record[[9]][4])
    out.list$loop <- as.integer(record[[10]][4])
    out.list$star <- as.integer(record[[11]][4])
    begin <- 12
  }
  # parse details. First find out the structure of the detail header
  if (record[[begin]][1] == "exp") {
    out.list$exp <- record[[begin]][2]
    begin <- begin + 1
  }
  out.list$pri_seq <- record[[begin]][2]
  out.list$pri_struct <- record[[begin+1]][2]
  begin <- begin +2
  
  c(out.list, ParseMrdDetail(record[-(1:begin)]))
  
}
ParseMrdDetail <- function(detailrecord) {
  # filter out rows without any information first
  detailrecord <- detailrecord[lapply(detailrecord, length) > 0]
  # Pull out the depth values
  depth <- lapply(detailrecord, function (x) {
    (strsplit(x[1], "\\_x"))[[1]][2]
  })
  depth <- as.integer(unlist(depth))
  
  # Make a matrix of the alignments
  alignments <- lapply(detailrecord, "[", 2)
  if (is.null(unlist(lapply(alignments[!is.na(depth)], strsplit, "")))) {
    return(list(depth=NULL, alignments=NULL))
  } 
  alignments <- t(matrix(unlist(lapply(alignments[!is.na(depth)], strsplit, "")), nrow=nchar(alignments[[1]])))
  # remove rows that don't contain a depth as they are not useful
  mirna.matches <- unlist(lapply(detailrecord[is.na(depth)], "[", 1))
  depth <- depth[!is.na(depth)]
  total.depth <- colSums((alignments != ".") * depth)
  mismatch.depth <- colSums(apply(alignments, c(1,2), grepl, pattern='[A-Z]') * depth)
  mismatched.reads <- sum(mismatch.depth)
  list(depth=depth, 
       alignments=alignments, 
       total.depth = total.depth, 
       mismatch.depth=mismatch.depth, 
       mismatched.reads=mismatched.reads, 
       mirna.matches=mirna.matches)
}
PlotMrd <- function(mrd, absolute=F) {
  # Plots an mrd record
  #
  # Args:
  #   mrd: The mrd object
  #   id: The id of the record you are interested in
  #   absolute: Plot as the absolute depth or fractional?
  # Returns:
  #   a ggplot object
  
  library(ggplot2)
  temp <- with(mrd, data.frame(pri_seq=unlist(strsplit(pri_seq, "")),
                                     pri_struct=unlist(strsplit(pri_struct, "")), 
                                     depth=total.depth))
  # if we are working in fractional depth, scale the depth
  temp$depth <- with(temp, 
                     if (absolute)
                       depth 
                     else 
                       depth/max(depth)
  )
  temp$loc <- 1:nrow(temp)
  
  # There are some magic numbers to adjust the positioning of the graph so that
  # we don't get the bases or folding overlapping with the graph.
  plot <- ggplot(temp, aes(loc, depth)) + geom_bar(stat="identity") + 
    geom_text(aes(label=pri_seq), y=1.07*max(temp$depth), size=3) +
    geom_text(aes(label=pri_struct), y=1.04*max(temp$depth), size=3) +
    theme(legend.position="none") + 
    coord_cartesian(ylim=c(-0.025*max(temp$depth), 1.1*max(temp$depth)))
  plot
}
PlotMrdCompare <- function(mrd1, mrd2, absolute=F) {
  # Plots an mrd record
  #
  # Args:
  #   mrd: The mrd object
  #   id: The id of the record you are interested in
  #   absolute: Plot as the absolute depth or fractional?
  # Returns:
  #   a ggplot object
  
  library(ggplot2)
  temp <- with(mrd1, data.frame(pri_seq=unlist(strsplit(pri_seq, "")),
                               pri_struct=unlist(strsplit(pri_struct, "")), 
                               depth1=total.depth))
  temp$depth2 <- mrd2$total.depth
  # if we are working in fractional depth, scale the depth
  if (!absolute) {
    temp$depth1 <- temp$depth1/max(temp$depth1)
    temp$depth2 <- temp$depth2/max(temp$depth2)
  }
  temp$loc <- 1:nrow(temp)
  
  # There are some magic numbers to adjust the positioning of the graph so that
  # we don't get the bases or folding overlapping with the graph.
  plot <- ggplot(temp, aes(loc, depth1)) + geom_line() + geom_line(aes(y=depth2,), color="red")
  plot
}

PlotMrdMM <- function(mrd, absolute=F) {
  # Plots an mrd record
  #
  # Args:
  #   mrd: The mrd object
  #   id: The id of the record you are interested in
  #   absolute: Plot as the absolute depth or fractional?
  # Returns:
  #   a ggplot object
  
  library(ggplot2)
  temp <- with(mrd, data.frame(depth=total.depth-mismatch.depth, MM=F))
  temp <- rbind(temp, with(mrd, data.frame(depth=mismatch.depth, MM=T)))
  temp$loc <- 1:length(mrd$total.depth)
  
  
  # if we are working in fractional depth, scale the depth
  temp$depth <- with(temp, 
                     if (absolute)
                       depth 
                     else 
                       depth/max(depth)
  )
  
  # There are some magic numbers to adjust the positioning of the graph so that
  # we don't get the bases or folding overlapping with the graph.
  plot <- ggplot(temp, aes(loc, depth, fill=MM)) + geom_bar(stat="identity")
    
  plot
}

CompareMRD <- function(mrd1, mrd2, id) {
  # Compare two different mirDeep analyses
  #
  # Args:
  #   mrd1: The mrd1 object
  #   mrd2: The mrd2 object
  #   id: The id of the record you are interested in
  # Returns:
  #   a ks test object
  ecdf1 <- cumsum(mrd1[[id]]$depth)
  ecdf1 <- ecdf1/max(ecdf1)
  ecdf2 <- cumsum(mrd2[[id]]$depth)
  ecdf2 <- ecdf2/max(ecdf2)
  
  ks.test(ecdf1,ecdf2)
}

CompareIDs <- function(data1, data2) {
  # Compare two different mirDeep analyses
  #
  # Args:
  #   mrd1: The mrd1 object
  #   mrd2: The mrd2 object
  #   id: The id of the record you are interested in
  # Returns:
  #   a ks test object
  ecdf1 <- cumsum(data1$total.depth)
  ecdf1 <- ecdf1/max(ecdf1)
  ecdf2 <- cumsum(data2$total.depth)
  ecdf2 <- ecdf2/max(ecdf2)
  ks.test(ecdf1,ecdf2)
}

loadMirBaseHTML <- function(html.in) {
  out.list <- list()
  library(XML)
  doc <- htmlParse(html.in)
  mirbase.table <- xpathSApply(doc, "//div[contains(@class, 'readSequence')]", xmlValue)
  mirbase.table <- unlist(strsplit(mirbase.table, "\\n"))
  mirbase.table <- strsplit(mirbase.table, "\\s+")
  # filter out the blank rows
  mirbase.table <- mirbase.table[lapply(mirbase.table, length) > 0]
  # grab the names
  out.list$mirna.matches <- mirbase.table[[1]][c(2,3)]
  out.list$pri_seq <- mirbase.table[[length(mirbase.table)-1]]
  out.list$pri_struct <- mirbase.table[[length(mirbase.table)]][1]
  
  mirbase.table <- mirbase.table[2:(length(mirbase.table)-2)]
  
  # depth vector
  depth <- lapply(mirbase.table, "[", 2)
  depth <- as.integer(unlist(depth))
  
  # Make a matrix of the alignments
  alignments <- lapply(mirbase.table, "[", 1)
  alignments <- t(matrix(unlist(lapply(alignments, strsplit, "")), nrow=nchar(alignments[[1]])))
  total.depth <- colSums((alignments != ".") * depth)
  
  c(out.list, list(depth=depth, 
       alignments=alignments, 
       total.depth = total.depth))
}

ParseMirbaseDetail <- function(detailrecord) {
  # filter out rows without any information first
  detailrecord <- detailrecord[lapply(detailrecord, length) > 0]
  # Pull out the depth values
  depth <- lapply(detailrecord, function (x) {
    (strsplit(x[1], "\\_x"))[[1]][2]
  })
  depth <- as.integer(unlist(depth))
  
  # Make a matrix of the alignments
  alignments <- lapply(detailrecord, "[", 2)
  if (is.null(unlist(lapply(alignments[!is.na(depth)], strsplit, "")))) {
    return(list(depth=NULL, alignments=NULL))
  } 
  alignments <- t(matrix(unlist(lapply(alignments[!is.na(depth)], strsplit, "")), nrow=nchar(alignments[[1]])))
  # remove rows that don't contain a depth as they are not useful
  mirna.matches <- unlist(lapply(detailrecord[is.na(depth)], "[", 1))
  depth <- depth[!is.na(depth)]
  total.depth <- colSums((alignments != ".") * depth)
  mismatch.depth <- colSums(apply(alignments, c(1,2), grepl, pattern='[A-Z]') * depth)
  mismatched.reads <- sum(mismatch.depth)
  list(depth=depth, 
       alignments=alignments, 
       total.depth = total.depth, 
       mismatch.depth=mismatch.depth, 
       mismatched.reads=mismatched.reads, 
       mirna.matches=mirna.matches)
}

ComputeNucleotideVariations <- function(mrd.rec) {
  library(foreach)
  library(iterators)
  library(plyr)
  primary.seq <- unlist(strsplit(mrd.rec$pri_seq, ""))
  if(is.null(mrd.rec$alignments)) {
    return(mrd.rec)
  }
  mod.table <- primary.seq == t(mrd.rec$alignments) # modified bases are False
  mod.table <- mod.table | t(mrd.rec$alignments == ".") # ignore bases not in the sequence
  mod.indexes <- which(!(mod.table | t(mrd.rec$alignments == ".")), arr.ind = T)
  mod.indexes <- data.frame(mod.indexes)
  out.table <- foreach(idx=iter(mod.indexes, by="row"), .combine=rbind) %do% {
    nucleotide <- mrd.rec$alignments[idx$col, idx$row]
    pos <- idx$row
    count <- mrd.rec$depth[idx$col]
    data.frame(pos=pos, nucleotide=nucleotide, count=count)
  }
  # each length polymorphism results in a different line so we can collapse these
  out.table <- ddply(out.table, .(pos, nucleotide), .fun=function(x) { data.frame(count=sum(x$count))})
  mrd.rec$snv <- out.table
  mrd.rec
}
ComputeLengthVariations <- function(mrd.rec) {
  library(foreach)
  library(iterators)
  library(plyr)
  primary.seq <- unlist(strsplit(mrd.rec$pri_seq, ""))
  if(is.null(mrd.rec$alignments)) {
    return(mrd.rec)
  }
  model <- unlist(strsplit(mrd.rec$exp, ""))
  
  # These models mark the positions of canonical 5 and 3 miRNA pieces.  M is for the other miRNAs
  model5 <- model == "5"
  model3 <- model == "3"
  modelM <- model == "M"

  # var.matrix has the positions which are not in the sequences
  var.matrix <- t(mrd.rec$alignments == ".")
  # xor here will give us a false if there are alignments outside of the
  # canonical sequence OR alignments which lack bases in the canonical.
  # these var matrixes 
  var5 <- xor(model5, var.matrix)
  var3 <- xor(model3, var.matrix)
  varM <- xor(modelM, var.matrix)
  
  # summary of lengths
  lengths <- data.frame(length=rowSums(!t(var.matrix)), count=mrd.rec$depth)
  lengths <- ddply(lengths, .(length), .fun=function(x) { data.frame(count=sum(x$count))})
  
  mrd.rec$lengths <- lengths
  
  # bin reads into 5 or 3
  if (!is.null(model5)) test5 <- colSums(!var5)
  if (!is.null(model3)) test3 <- colSums(!var3)
  if (!is.null(modelM)) testM <- colSums(!varM)
  bin5 <- test5 < test3
  bin3 <- test3 < test5
  binM <- (testM <= test5) & (testM <= test3)
  
  depth5 <- (!var5) %*% matrix(mrd.rec$depth * bin5, ncol=1)
  depth3 <- (!var3) %*% matrix(mrd.rec$depth * bin3, ncol=1)
  depthM <- (!varM) %*% matrix(mrd.rec$depth * binM, ncol=1)
  

  mrd.rec$length.poly.5 <- as.vector(depth5)
  mrd.rec$length.poly.3 <- as.vector(depth3)
  mrd.rec$length.poly.M <- as.vector(depthM)
  #if we have no 5p or 3p pieces the following will give a warning.  The results
  #are fine so we don't need to worry about these.
  suppressWarnings(mrd.rec$range.5p <- range(which(model5)))
  suppressWarnings(mrd.rec$range.3p <- range(which(model3)))
  suppressWarnings(mrd.rec$range.M <- range(which(modelM)))
  mrd.rec
}
PlotNucleotideVariants <- function(mrd) {
  plot.data <- mrd$snv
  if(is.finite(sum(mrd$range.3p)) & is.finite(sum(mrd$range.5p))) {
    cutoff <- (mrd$range.3p[1] + mrd$range.5p[2])/2
    plot.data$plot <- ifelse(plot.data$pos < cutoff, "5p", "3p")
    plot.data$pos.corrected <- ifelse(plot.data$pos < cutoff, plot.data$pos-mrd$range.5p[1], plot.data$pos-mrd$range.3p[1])
    plot.data$AF <- plot.data$count/mrd$total.depth[plot.data$pos]
  } else {
    if (is.finite(sum(mrd$range.3p))) {
      plot.data$plot <- "3p"
      range <- mrd$range.3p
    } else if (is.finite(sum(mrd$range.5p))) {
      plot.data$plot <- "5p"
      range <- mrd$range.5p
    } else if (is.finite(sum(mrd$range.M))) {
      plot.data$plot <- "M"
      range <- mrd$range.M
    } else {
      stop("no data in range to plot")
    }
    plot.data$pos.corrected <- plot.data$pos-range[1]
    plot.data$AF <- plot.data$count/mrd$total.depth[plot.data$pos]
  }
  ggplot(plot.data, aes(pos.corrected, count, fill=nucleotide)) + geom_bar(stat="identity") + facet_wrap(~plot, ncol=1)
}
PlotLengthDistribution <- function(mrd) {
  ggplot(mrd$lengths, aes(factor(length), count)) + geom_bar(stat="identity")
}
PlotLengthPolymorphisms <- function(mrd) {
  library(reshape2)
  if(is.finite(sum(mrd$range.3p)) & is.finite(sum(mrd$range.5p))) {
    cutoff <- (mrd$range.3p[1] + mrd$range.5p[2])/2
    plot.data <- data.frame(pos=1:length(mrd$length.poly.5), poly.5=mrd$length.poly.5, poly.3=mrd$length.poly.3)
    plot.data <- melt(plot.data, "pos")
    suppressWarnings(plot.data$cat <- factor(cut(plot.data$pos, c(-Inf, mrd$range.5p[1]-0.5, mean(mrd$range.5p), mrd$range.5p[2]-0.5, 
                                          cutoff,
                                          mrd$range.3p[1]-0.5, mean(mrd$range.3p), mrd$range.3p[2]-0.5, Inf), 
                         #c("55e", "55r", "53r", "53e", "35e", "35r", "33r", "33e")
                         c("5' extension", "5' retraction", "3' retraction", "3' extension",
                           "5' extension", "5' retraction", "3' retraction", "3' extension")
                         )))
    plot.data$pos <- ifelse(plot.data$pos < cutoff, plot.data$pos-mrd$range.5p[1], plot.data$pos-mrd$range.3p[1])
    plot.data <- plot.data[plot.data$value > 0, ]
    plot.data$pos <- ifelse(plot.data$cat %in% c("5' retraction", "3' extension"), plot.data$pos + 1, plot.data$pos)
    ggplot(plot.data, aes(pos, value, fill=cat)) + geom_bar(stat="identity") + facet_wrap(~variable, ncol=1, scales="free_y")
  } else {
    plot.data <- data.frame(pos=1:length(mrd$length.poly.5), poly.5=mrd$length.poly.5, poly.3=mrd$length.poly.3, poly.M=mrd$length.poly.M)
    if (is.finite(sum(mrd$range.3p))) {
      range <- mrd$range.3p
    } else if (is.finite(sum(mrd$range.5p))) {
      range <- mrd$range.5p
    } else if (is.finite(sum(mrd$range.M))) {
      range <- mrd$range.M
    } else {
      stop("no data in range to plot")
    }
    plot.data <- melt(plot.data, "pos")
    plot.data$cat <- factor(cut(plot.data$pos, c(-Inf, range[1]-0.5, mean(range), range[2]-0.5, Inf), 
                                                 c("5' extension", "5' retraction", "3' retraction", "3' extension")))
    plot.data$pos <- plot.data$pos-range[1]
    plot.data$pos <- ifelse(plot.data$cat %in% c("5' retraction", "3' extension"), plot.data$pos + 1, plot.data$pos)
    plot.data <- plot.data[plot.data$value > 0, ]
  }
  ggplot(plot.data, aes(pos, value, fill=cat)) + geom_bar(stat="identity") + facet_wrap(~variable, ncol=1, scales="free_y")
}
