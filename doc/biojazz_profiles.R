#!/usr/bin/R
# this scirpt is used to plot the transformation of profiles in BioJazz model

library(R.utils);

################ generate the plots ################


plotScaling <-
    function (stringLength, scaling, min, max)
    {
        decVector = c(0:(2^stringLength - 1));
        distanceVector = c(0:stringLength);
        vectorLength = length(decVector);
        binVector = array(, vectorLength);
        scaledVector = array(, vectorLength);
        for (i in c(1:vectorLength)) {
            tempBin = sapply(strsplit(paste(intToBits(decVector[i])),""),'[[',2);
            tempBin = rev(tempBin[1:stringLength]);
            binVector[i] = paste(tempBin[], collapse="");

        }

        if (scaling == "logLinearField") {
            scaledVector <- sapply(decVector, function(x) {min * (max / min) ^ (x/(vectorLength-1))});


        } else if (scaling == "linearField") {
            scaledVector <- sapply(decVector, function(x) {min + (max - min) * (x/(vectorLength-1))});


        } else if (scaling == "logLinearProfile") {
            # how to plot this one?
            scaledVector <- sapply(distanceVector, function(x) {min * (max / min) ^ (x / stringLength)});

        } else {
            stop("specify how to scale your fields and profiles!");
        }

    }

temp <- plotScaling(20, "logLinearField", 1e-3, 1e3);

plot(decVector,temp);
