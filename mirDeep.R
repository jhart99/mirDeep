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
  record.list <- lapply(record.list, function(x) {
    out.list <- list()
    record <- strsplit(x, "\\s+")

    out.list$total <- as.integer(record[[1]][4])
    out.list$left <- as.integer(record[[2]][4])
    out.list$right <- as.integer(record[[3]][4])
    out.list$other <- as.integer(record[[4]][4])
    if (record[[5]][1] == "exp") {
      out.list$exp <- record[[5]][2]
      out.list$pri_seq <- record[[6]][2]
      out.list$pri_struct <- record[[7]][2]
      begin <- 8
    } else {
      out.list$exp <- NA
      out.list$pri_seq <- record[[5]][2]
      out.list$pri_struct <- record[[6]][2]
      begin <- 7
    }
    out.list$depth <- vector("integer", length(out.list$pri_seq))
    for (line in record[-(1:begin)]) {
      if (length(line) < 2 ) next
      id <- strsplit(line[1], "\\_x")
      count <- id[[1]][2]
      # This next line parses the matches for stuff that isn't "." to give a
      # vector of true/false values.  This is multiplied by the observed count
      # to make an integer vector and then added to the existing totals
      counts <- (strsplit(line[2], "")[[1]] != ".") * as.integer(count)
      out.list$depth <- out.list$depth + counts
    }
    out.list
  })
  
  record.list
}

PlotMrd <- function(mrd, id) {
  # Plots an mrd record
  #
  # Args:
  #   mrd: The mrd object
  #   id: The id of the record you are interested in
  # Returns:
  #   a ggplot object
  
  library(ggplot2)
  temp <- with(mrd[[id]], data.frame(exp=unlist(strsplit(exp, "")), 
                                     pri_seq=unlist(strsplit(pri_seq, "")),
                                     pri_struct=unlist(strsplit(pri_struct, "")), 
                                     depth=depth))
  temp$loc <- 1:nrow(temp)
  
  plot <- ggplot(temp, aes(loc, depth)) + geom_line() + 
    geom_point(aes(color=exp)) + 
    geom_text(aes(label=pri_seq, color=exp), y=1.07*max(temp$depth), size=3) +
    geom_text(aes(label=pri_struct, color=exp), y=1.04*max(temp$depth), size=3) +
    theme(legend.position="none") + 
    coord_cartesian(ylim=c(-0.025*max(temp$depth), 1.1*max(temp$depth))) +
    ggtitle(id)
  plot
}
