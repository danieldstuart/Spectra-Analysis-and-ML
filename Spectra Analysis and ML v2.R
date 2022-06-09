# List of packages 
packages = c("dplyr", "rfUtilities", "caret", "R.utils", "pheatmap", "MALDIquantForeign", "MLeval", "MASS", "MALDIquant", "caretEnsemble", "corrplot", "readMzXmlData")
library(caret)
## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

setwd("~/")
setwd("~/R/Spectra Analysis and ML/Place Data Here") #make sure your data is placed in the correct folder separated into folders by class

#foldernum <- readline(prompt="How many data folders/classes do you have? ")
#foldernum <- as.integer(foldernum)
foldernum = 6 #sets the number of folders to pull from, the number of classes
folders <- list.files(full.names = T, include.dirs = T) #pulls out all of the folder names
classes <- sub("..", "", folders) #create class column to append to data

for (x in 1:foldernum){
  setwd("~/R/Spectra Analysis and ML/Place Data Here")
  setwd(folders[x]) #changes working directory to next folder 
  directory <- list.files(pattern="*.mzxml", recursive = T) #pulls out spectra within selected folder
  s1 <- importMzXml(directory)
  outtest <- exists("s")
  if(outtest == FALSE){
  s <- s1 #combines all the pulled spectra together
  classifier <- as.character(c(rep(classes[x], length(s1)))) #creates classifier table
  }
  s <- c(s,s1)
  classifier1 <- as.character(c(rep(classes[x], length(s1)))) #creates classifer1 table as temporary store to append
  classifier <- as.character(c(classifier, classifier1)) #combines class names from each iteration
}

classifier <- as.factor(classifier) #converts classes to factors

## Select mass range of interest
minmass <- 500
maxmass <- 2400
#minmass <- readline(prompt="Enter minimum m/z value: ")
#minmass <- as.integer(minmass)
#maxmass <- readline(prompt="Enter maximum m/z value: ")
#maxmass <- as.integer(maxmass)

## Process and align all the spectra
spectra <- smoothIntensity(s, method="MovingAverage", halfWindowSize=2)
spectra <- removeBaseline(spectra, method="SNIP", iterations=100)
spectra <- alignSpectra(spectra, halfWindowSize=20, SNR=4, tolerance=0.002, warpingMethod="lowess")

plot(spectra[[1]], xlim=c(minmass, maxmass), ylim=c(0, 2000)) #plots the first spectra 

spectra <- calibrateIntensity(spectra, method = "TIC") #normalize intensity

## Detect peaks for all spectra and plot result for first spectra
peaks <- detectPeaks(spectra, method="MAD", halfWindowSize=20, SNR=5)

## Equalize similar peaks, remove inconsistent peaks, combine into dataset of peak(m/z) and intensity 
peaks <- binPeaks(peaks, method = "strict", tolerance=0.002)
peaks <- filterPeaks(peaks, minFrequency=0.1)
imatrix <- intensityMatrix(peaks, spectra)

## Keep only peaks between the min and maxmass
cols <- as.integer(colnames(imatrix))
col_to_keep <- cols > minmass & cols < maxmass
imatrix <- imatrix[,col_to_keep]

## Add classification column to dataset
dataset <- data.frame(imatrix, classifier, stringsAsFactors = TRUE)
datasetsave <- dataset

## Split data into training(dataset) and testing(validation)
validation_index <- createDataPartition(dataset$classifier, p=0.70, list=FALSE)
validation <- dataset[-validation_index,]
dataset <- dataset[validation_index,]

## Train control 
control <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions = "final", classProbs = TRUE, index = createFolds(dataset$classifier, 10), allowParallel =  TRUE)

## Set metric for model success, Kappa can be good for low % of samples in 1 class
metric <- "Accuracy"

## Train models of chosen type's
algorithms <- c("rf", "nnet", "kknn", "lda") #algorithms to be tested
model_list <- caretList(classifier~., data=dataset, metric=metric, trControl=control, methodList=algorithms)

# Sets working directory so results can be output
setwd("~/")
setwd("~/R/Spectra Analysis and ML/Results Output Here")

graphics.off()

## Output results of the training and compare success metrics
results <- resamples(model_list)
summary(results)
pdf("models comparisson.pdf", width = 7, height = 7, onefile = FALSE) #saves results as a pdf of given name with given size (in inches)
dotplot(results) #plot to be saved
dev.off() #completes the pdf save

for (methodname in algorithms){
  ## Predict testing data with chosen model
  p <- predict(model_list[[methodname]], validation)
  cm <- confusionMatrix(p, validation$classifier)
  print(cm)
  
  tryCatch({
    
  ## Output heatmap of prediction results
  mat <- as.table(cm)
  pdf(paste(methodname, "heatmap.pdf"), width = 7, height = 7, onefile = F) #saves results as a pdf of given name with given size (in inches)
  map <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, show_colnames = TRUE, show_rownames = TRUE, legend = TRUE, fontsize_number = 20)
  dev.off()
  
  ## Output correlation plot of prediction results
  co <- cor(mat)
  pdf(paste(methodname, "correlation plot.pdf"), width = 7, height = 7, onefile = FALSE) #saves results as a pdf of given name with given size (in inches)
  corrplot(co, method = "circle")
  dev.off()
  
  }, error=function(e){cat("EROR :", conditionMessage(e), "\n")})
}

# Outputs a pdf of important variables
pdf(paste("rf", "important variables.pdf"), width = 7, height = 7, onefile = FALSE)
impvar <- varImp(model_list[["rf"]], scale = TRUE)
plot(impvar, top = 20)
dev.off()

# Trains an LDA model using the whole data set
lda_model <- lda(classifier ~ ., data = datasetsave)
lda_predict <- predict(lda_model)

# Plots LDA results
ggplot(cbind(datasetsave, lda_predict$x),
       aes(y = LD1, x = LD2, colour = classifier)) + 
  stat_ellipse(aes(fill = classifier), geom = "polygon", alpha = .2, level = 0.99) +
  geom_point() +
  ggtitle("LDA") +
  theme(legend.position = "right")

#Predict LDA using validation data set to blind
val <- validation
val$classifier <- "Blind"
lda_val_predict <- predict(object = lda_model, newdata = val)

# Plots LDA results with blind
ggplot(cbind(datasetsave, lda_predict$x), 
       aes(y = LD1, x = LD2, colour = classifier)) +
    stat_ellipse(aes(fill = classifier), geom = "polygon", alpha = .3, level = 0.99) +
    geom_point() +
    geom_point(data = cbind(val, lda_val_predict$x), colour = "black") +
    ggtitle("LDA") +
    theme(legend.position = "right")

# Outputs a pdf of LDA results
pdf(paste("LDA", "results.pdf"), width = 7, height = 7, onefile = FALSE)
ggplot(cbind(datasetsave, lda_prediction$x),
       aes(y = LD1, x = LD2, colour = classifier)) +
  stat_ellipse(aes(fill = classifier), geom = "polygon", alpha = .2, level = 0.99) +
  geom_point() +
  ggtitle("LDA") +
  theme(legend.position = "right")
dev.off()

## Plotting with plotly
library(plotly)

# Setting up lda data for plotly
lda_data <- cbind(datasetsave, lda_predict$x)
LD1 <- lda_data$LD1
LD2 <- lda_data$LD2
LD3 <- lda_data$LD3
class <- lda_data$classifier
lda_dat <- data.frame(LD1, LD2, LD3, class, stringsasfactors = TRUE)

# Setting up plot properties
t1 <- list(
  family = "Times New Roman",
  size = 24,
  color = "black"
)

t2 <- list(
  family = "Times New Roman",
  size = 18,
  color = "black")

axx <- list(title = "LD2")
axy <- list(title = "LD1")
axz <- list(title = "LD3")

# Plotly figure
fig <- plot_ly(lda_dat, x = LD2, y = LD1, z = LD3, type="scatter3d", mode="markers", color=classifier, size=3)%>%
  layout(title = list(text = "Bacterial Species LDA", font = t1, y = 0.95), scene = list(xaxis=axx,yaxis=axy,zaxis=axz), legend = list(title = list(text = "Species"), font = t2))
# Plot
fig

setwd("~/")
#graphics.off()