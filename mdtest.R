test <- ReadMrd("miRBase.mrd")

test[[1]]  # to see the structure of the individual records
test$`hsa-let-7d`  # you can also pull out indiviual miRNAs using this format

barplot(test[[1]]$depth)  # you can plot the depth using a barplot
plot(test[[1]]$depth) # or the normal plotting system

<<<<<<< HEAD
PlotMrd(test, 'hsa-let-7d') # pretty plots
=======
PlotMrd(test, foo) # pretty plots
>>>>>>> 03fd82379d4681fefbbefeda7ed113bb14b9cc69
