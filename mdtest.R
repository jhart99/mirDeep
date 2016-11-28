test <- ReadMrd("miRBase.mrd")
test2 <- ReadMrd("mirDeep2.mrd")

test[[1]]  # to see the structure of the individual records
test[['hsa-let-7d']]  # you can also pull out indiviual miRNAs using this format

barplot(test[[1]]$total.depth)  # you can plot the depth using a barplot
plot(test[[1]]$total.depth) # or the normal plotting system

PlotMrd(test[[1]]) # pretty plots


# plot a miRNA showing mismatches (possibly RNA editing)
PlotMrdMM(test[['hsa-let-7d']])

# compare two miRNAs
CompareIDs(test[[1]], test[[2]])

# Load miRNA data from miRbase (just save the page and load it)
toot <- loadMirBaseHTML("miRNA Search Results.html")
PlotMrd(toot)
# or you can just copy and paste the URL
library(RCurl)
toot <- loadMirBaseHTML("http://www.mirbase.org/cgi-bin/get_read.pl?acc=MI0000269")
PlotMrd(toot)

