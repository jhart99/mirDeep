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