#!/usr/bin/R
# this script is to plot epistasis of the profiles
library(R.methodsS3);
library(R.utils);

################# get the scaling ################

dir.create("epistasis", showWarnings = FALSE);
system("cd epistasis");

binStrToDec <- 
  function (binaryArray)
  {
    convertedDec = sum(2^(which(rev(binaryArray) == 1) - 1));
    return(convertedDec);
  }

findRm0 <- 
  function (binaryArray0, binaryArray1, binaryArray2)
  {
    length = length(binaryArray1);
    if (length(binaryArray2) != length || length(binaryArray0) != length) {
      stop("two strings are not equally long!");
    }
    tempBinary = array(, length);
    for (m in c(1:length)) {
      if (binaryArray1[m] != binaryArray0[m]) {
        tempBinary[m] = binaryArray1[m];
      } else if (binaryArray2[m] != binaryArray0[m]) {
        tempBinary[m] = binaryArray2[m];
      } else {
        tempBinary[m] = binaryArray0[m];
      }
    }
    return(sum(2^(which(rev(tempBinary) == 1) - 1)));
  }

calculateEpistasis <-
  function (stringLength, scalingMethod, min, max)
  {
    decVector = c(0:(2^stringLength - 1));
    distanceVector = c(0:stringLength);
    vectorLength = length(decVector);
    binVector = array(, dim = c(vectorLength, stringLength));
    scaledVector = array(, vectorLength);
    for (i in c(1:vectorLength)) {
      tempBin = sapply(strsplit(paste(intToBits(decVector[i])),""),'[[',2);
      tempBin = rev(tempBin[1:stringLength]);
      binVector[i,] = tempBin[];
    }
    
    binVectorAsInt = array(, dim = c(vectorLength, stringLength));
    for (i in c(1:stringLength)) {
      binVectorAsInt[,i] = strtoi(binVector[,i]);      
    }

    if (scalingMethod == "logLinearField") {
      scaledVector <- sapply(decVector, function(x) {min * (max / min) ^ (x/(vectorLength-1))});
      
      Rm0 = NULL;
      Rm1 = NULL;
      Rm2 = NULL;
      epista = NULL;
      mCal = 1;
      for (i in c(1:vectorLength)) {
        WT = binStrToDec(binVectorAsInt[i,]);
        mutationPool = array(, dim = c(stringLength, stringLength));
        for (j in c(1:stringLength)) {
          mutationPool[j,] = binVectorAsInt[i,];
          mutationPool[j,j] = 1 - mutationPool[j,j]; 
        }
        for (k in c(1:(stringLength - 1))) {
          for (l in c(1:(stringLength - k))) {
            Rm1[mCal] = binStrToDec(mutationPool[k,]);
            Rm2[mCal] = binStrToDec(mutationPool[l,]);
            Rm0[mCal] = findRm0(binVectorAsInt[i,], mutationPool[k,], mutationPool[l,]);
            epista[mCal] = (scaledVector[(Rm1[mCal]+1)] - scaledVector[(WT+1)]) * (scaledVector[(Rm2[mCal] + 1)] - scaledVector[(WT+1)]) - (scaledVector[(Rm0[mCal]+1)] - scaledVector[(WT+1)]);
            mCal = mCal + 1;
          }        
        }
      }

    } else if (scalingMethod == "linearField") {
      scaledVector <- sapply(decVector, function(x) {min + (max - min) * (x/(vectorLength-1))});
  
      
      epista = NULL;
      mCal = 1;
      for (i in c(1:vectorLength)) {
        WT = binStrToDec(binVectorAsInt[i,]);
        mutationPool = array(, dim = c(stringLength, stringLength));
        for (j in c(1:stringLength)) {
          mutationPool[j,] = binVectorAsInt[i,];
          mutationPool[j,j] = 1 - mutationPool[j,j]; 
        }
        for (k in c(1:(stringLength - 1))) {
          for (l in c(1:(stringLength - k))) {
            Rm1 = binStrToDec(mutationPool[k,]);
            Rm2 = binStrToDec(mutationPool[l,]);
            Rm0 = findRm0(binVectorAsInt[i,], mutationPool[k,], mutationPool[l,]);
            epista[mCal] = (scaledVector[(Rm1+1)] - scaledVector[(WT+1)]) * (scaledVector[(Rm2 + 1)] - scaledVector[(WT+1)]) - (scaledVector[(Rm0+1)] - scaledVector[(WT+1)]);
            mCal = mCal + 1;
          }        
        }
      }
      
    } else if (scalingMethod == "logLinearProfile") {
      scaledVector <- sapply(distanceVector, function(x) {min * (max / min) ^ (x / stringLength)});
      
      stop("logLinearProfile is not implemented yet");
    } else {
      stop("specify how to scale your fields and profiles!");
    }
    
    
    return(epista);
  }


epistasis = calculateEpistasis(5, "logLinearField", 1e-2, 1e1);

plot(density(epistasis));



