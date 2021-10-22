

###  USing a OTSSC Package functions but modified for issue with 0s being generated from our data

# First function  - modified in line 200

JillOSTSC <- function(sample, label, class, ratio = 1, per = 0.8, r = 1, k = 5, 
                  m = 15, parallel = TRUE, progBar = TRUE) {
  # Oversample a time series sequence imbalance data.
  #
  # Args:
  #   sample:       Univariate sequence data samples.
  #   label:        Labels corresponding to samples.
  #   class:        The number of the classes to be oversampled, starting from the class with the fewest observations, with the default setting to progress to as many classes as possible 
  #   ratio:        The oversampling ratio 
  #              number (>=1) (default = 1)
  #   per:          Ratio of weighting between ESPO and ADASYN (default = 0.8) 
  #   r:            A scalar ratio specifying which level (towards the boundary) we shall push the synthetic data (in EPSO, default = 1)
  #   k:            k Number of nearest neighbours in k-NN (for ADASYN) algorithm (default = 5)
  #   m:            m Seeds from the positive class in m-NN (for ADASYN) algorithm (default = 15)
  #   parallel:     parallel Whether to execute in parallel mode (default = TRUE). 
  #                 (Recommended for datasets with over 30,000 records.)
  #   progBar:      Whether to include progress bars (default = TRUE).
  #
  # Returns:
  #   The oversampled dataset samples data_list$sample and labels data_list$label.
  
  # check if the input sample data had two dimension
  if (is.null(dim(sample)) || length(dim(sample)) != 2) {
    stop ("The input sample data must have two dimensions.")
  }
  
  # check if the numbers of records in label and sample matched
  if (is.null(dim(label))) {
    sizeLabel <- length(label)
  } else {
    sizeLabel <- dim(label)[1]
  }
  
  sizeSample <- dim(sample)[1]
  
  if (sizeLabel != sizeSample) {
    stop ("Number of time series sequences provided in sample do not match the 
          number of classes provided in label. Check dimensions.")
  }
  
  # check if the class input is in the numeric format
  if (!missing(class) && !is.numeric(class)) {
    stop ("The parameter class is not in correct format, which must be a numeric value.")
  }  
  
  # check if the ratio input is in the numeric format
  if (!is.numeric(ratio)) {
    stop ("The parameter ratio is not in correct format, which must be a numeric value.")
  }    
  
  # check if the Percentage input is in the numeric format
  if (!is.numeric(per)) {
    stop ("The parameter per is not in correct format, which must be a numeric value.")
  } 
  
  # check if the r input is in the numeric format
  if (!is.numeric(r)) {
    stop ("The parameter r is not in correct format, which must be a numeric value.")
  }
  
  # check if the k input is in the numeric format
  if (!is.numeric(k)) {
    stop ("The parameter k is not in correct format, which must be a numeric value.")
  }
  
  # check if the m input is in the numeric format
  if (!is.numeric(m)) {
    stop ("The parameter m is not in correct format, which must be a numeric value.")
  }
  
  # check if the class input is only one element
  if (!missing(class) && length(class) != 1) {
    stop ("The parameter class is not in correct format, which must be a single value.")
  }
  
  # check if the ratio input is only one element
  if (length(ratio) != 1) {
    stop ("The parameter ratio is not in correct format, which must be a single value.")
  }
  
  # check if the per input is only one element
  if (length(per) != 1) {
    stop ("The parameter per is not in correct format, which must be a single value.")
  }
  
  # check if the R input is only one element
  if (length(r) != 1) {
    stop ("The parameter r is not in correct format, which must be a single value.")
  }
  
  # check if the k input is only one element
  if (length(k) != 1) {
    stop ("The parameter k is not in correct format, which must be a single value.")
  }
  
  # check if the m input is only one element
  if (length(m) != 1) {
    stop ("The parameter m is not in correct format, which must be a single value.")
  }
  
  # check if the ratio input is in range (0,1]
  if (ratio > 1 || ratio <= 0) {
    stop ("The parameter ratio is not in correct range, which must be betwwen 
          0 to 1, including 1.")
  }
  
  # check if the Percentage input is in range [0,1]
  if (per > 1 || per < 0) {
    stop ("The parameter per is not in correct range, which must be betwwen 
          0 to 1, including 0 and 1.")
  }
  
  # check if the R input is in range [1,+oo)
  if (r < 1) {
    stop ("The parameter r is not in correct range, which must be larger or 
          equal to 1.")
  }
  
  # check if the k input is in range (0,+oo)
  if (k <= 0) {
    stop ("The parameter k is not in correct range, which must be larger than 0.")
  }
  
  # check if the m input is in range (0,+oo)
  if (m <= 0) {
    stop ("The parameter m is not in correct range, which must be larger than 0.")
  }
  
  # check if the k input is an integer
  if (k %% 1 != 0) {
    stop ("The parameter k is not in correct format, which must be an integer.")
  }
  
  # check if the m input is an integer
  if (m%%1 != 0) {
    stop ("The parameter m is not in correct format, which must be an integer.")
  }
  
  # check if the parallel input is a boolean value
  if (!(identical(parallel, FALSE) || identical(parallel, TRUE))) {
    stop ("The parameter parallel is not in correct format, which must be a 
          boolean value.")
  }
  
  # check if the progBar input is a boolean value
  if (!(identical(progBar, FALSE) || identical(progBar, TRUE))) {
    stop ("The parameter progBar is not in correct format, which must be a 
          boolean value.")
  }
  
  # combine labels and features
  fullData <- cbind(label, sample)
  fullData <- matrix(unlist(fullData, use.names = FALSE), 
                     ncol = ncol(fullData))
  
  # clean missing values and non-number values by removing their belonging rows
  fullData <- matrix(suppressWarnings(as.numeric(fullData)), 
                     nrow = nrow(fullData))
  cleanData <- na.omit(fullData)
  
  # determine how many classes need to be oversampled
  Lab <- cleanData[, c(1)]
  claTab <- as.data.frame(table(Lab))  # count frequency of classes
  claTab <- claTab[order(claTab$Freq), ]  # order in ascending
  
  sumFreq <- sum(claTab$Freq)
  
  count <- 0
  for (i in 1:dim(claTab)[1]) {
    if (sumFreq - claTab$Freq[i] > claTab$Freq[i]) {
      count <- count + 1
    }
  }
  
  if (count == 0) {
    stop ("The input dataset is already balanced. No oversampling is necessary.")
  }
  
  if (missing(class)) {
    class <- count
  } 
  
  if (count < class) {
    warning ("Insufficient observations of the minority class. The class number that needs to be oversampled 
             is set to ", count)
  }
  
  ##################  Start of my modifications   ##############################
  #  the function is modified to call the JillReguCovar function rather than the ReguCovar Function
  
  myData <- list()
  for (i in 1:class) {
    targetClass <- as.numeric(as.vector(claTab$Lab[i]))
    newData <- JillReguCovar(cleanData, targetClass, ratio, r, per, k, m, 
                         parallel, progBar)
  
    ##################################################################################  
    
    myData <- rbind(myData, newData)
  }
  
  nData <- list()
  for (i in (class + 1):dim(claTab)[1]) {
    targetClass <- as.numeric(as.vector(claTab$Lab[i]))
    nega <- cleanData[which(cleanData[, c(1)] == targetClass), ]
    nData <- rbind(nData, nega)
  }
  
  # form data
  dataNew <- rbind(myData, nData)
  dataNew <- matrix(unlist(dataNew), ncol=ncol(dataNew))
  
  dataX <- dataNew[, -1]
  dataY <- dataNew[, c(1)]
  dataList <- list("sample" = dataX, "label" = dataY)
  
  return(dataList)
}  


##  Second modified function

JillReguCovar = function(cleanData, targetClass, ratio, r, per, k, m, parallel, progBar) {
  # Generate samples by ESPO and ADASYN.
  #
  # Args:
  #   cleanData:    First column is label data, rest is sample data, without missing values.
  #   targetClass: The class to be oversampled. 
  #   ratio:        The oversampling ratio 
  #              number (>=1) (default = 1)  
  #   r:            A scalar ratio specifying which level (towards the boundary) we shall push the synthetic data (in EPSO, default = 1) 
  #   per:          Ratio of weighting between ESPO and ADASYN (default = 0.8)  
  #   k:            k Number of nearest neighbours in k-NN (for ADASYN) algorithm (default = 5)
  #   m:            m Seeds from the positive class in m-NN (for ADASYN) algorithm (default = 15)
  #   parallel:     parallel Whether to execute in parallel mode (default = TRUE). 
  #                 (Recommended for datasets with over 30,000 records.)
  #   progBar:      Whether to include progress bars (default = TRUE).
  #  Returns:
  #   newData: the oversampled dataset.
  
  # form positive (target class) and negative data
  # The negative data is formed using a one-vs-rest strategy.
  positive <- cleanData[which(cleanData[, c(1)] == targetClass), ]
  
  negative <- cleanData[which(cleanData[, c(1)] != targetClass), ]
  
  p <- positive[, -1]  # remove label column
  n <- negative[, -1]
  
  # Number of sequences to be created
  nTarget <- nrow(n)*ratio
  
  poscnt <- nrow(p)
  if (nTarget > poscnt) { 
    # check if the positive data records have already more than the required number of records to be created
    
    # Compute Regularized Eigen Spectra
    numToGen <- ceiling((nTarget - poscnt)*per)
    numADASYN <- nTarget - poscnt - numToGen
    
    me <- apply(p, 2, mean)  # Mean vector of p
    pCov <- cov(p)  # vector covariance
    v <- eigen(pCov)$vectors  # Eigen axes matrix
    # v <- v[, n:1]
    d <- eigen(pCov)$values  # Eigenvalues
    # d <- d[n:1]
    numD <- ncol(p)  # The feature dimension
    ind <- which(d <= 0.005)  # The unreliable eigenvalues
    if (length(ind) != 0) {
      por <- ind[1]  # [1,por] the portion of reliable
    } else {
      por <- numD
    }
    tCov  <- cov(rbind(p, n))  # The covariance matrix of the total data (column)
    dT <- crossprod(v, tCov) %*% v  # dT = v' * tCov * v
    dT <- diag(dT)  # Turning the diagonal of matrix dT to a vector
    
    # Modify the Eigen spectrum according to a 1-Parameter Model
    # dMod: Modified Eigen Spectrum Value
    dMod <- matrix(0, 1, numD)
    alpha <- d[1]*d[por]*(por-1)/(d[1] - d[por])
    beta  <- (por*d[por] - d[1])/(d[1] - d[por])
    for (i in 1:numD) {
      if (i < por) {
        dMod[i] <- d[i]
      } else {
        dMod[i] = alpha/(i+beta)
        #############  modified here so that any values <0 are replaced#####################
        if (dMod[i] > dT[i] || dMod[i]<0) {
          dMod[i] <- dT[i]
        }
      }
    }
    # Create Oversampled Data by ESPO and ADASYN
    # Users choose if applying in parallel and if adding progress bar
    if (numToGen != 0) {
      if (identical(parallel, FALSE)) {
        if (identical(progBar, FALSE)) {
          sampleESPO <- ESPO(me, v, dMod, p, n, r, por, numToGen)
        } else {
          cat("Oversampling class", targetClass, "... \n")
          sampleESPO <- ESPOBar(me, v, dMod, p, n, r, por, numToGen)
        }
      } else {
        if (identical(progBar, FALSE)) {
          sampleESPO <- ESPOPara(me, v, dMod, p, n, r, por, numToGen)
        } else {
          cat("Oversampling class", targetClass, "... \n")
          sampleESPO <- ESPOParaBar(me, v, dMod, p, n, r, por, numToGen)
        }
      }
    }
    
    if (numADASYN != 0) {
      if (identical(parallel, FALSE)) {
        if (identical(progBar, FALSE)) {
          sampleADA <- ADASYN(t(p), t(n), numADASYN, k, m)
        } else {
          sampleADA <- ADASYNBar(t(p), t(n), numADASYN, k, m)
        }
      } else {
        if (identical(progBar, FALSE)) {
          sampleADA <- ADASYNPara(t(p), t(n), numADASYN, k, m)
        } else {
          sampleADA <- ADASYNParaBar(t(p), t(n), numADASYN, k, m)
        }
      }
    }
    
    # Form new data
    dataTargetClass <- rbind(t(sampleADA), sampleESPO)
    newData <- cbind(matrix(targetClass, nTarget, 1), dataTargetClass)
    return(newData)
  } else {
    return(positive)
  }
}



# remaining functions - not modified

#' Generate samples by ADASYN approach.
#' 
#' @param p minority class samples
#' @param n majority class samples
#' @param nTarget the targeted number of samples to achieve
#' @param k is the number of nearest neighbours in the ADASYN algorithm, with the default value of 5
#' @param m seeds from the positive class in k-NN of the ADASYN algorithm, with the default value of 15
#' @return sampleADA
#' @importFrom fields rdist 
#' @importFrom stats runif 
#' @keywords internal

ADASYN <- function(p, n, nTarget, k, m) {
  # Generate samples by ADASYN.
  #
  # Args:
  #   p:       The minority class samples.
  #   n:       The majority class samples. P and N must have the same feature dimension, greater than one,
  #            with no missing values.
  #   nTarget: The targeted number of samples to achieve.
  #   k:       k-NN search used in the ADASYN algorithm, with a default value of 5.
  #   m:       m-NN search used in the ADASYN, finding seeds from the Positive Class, with the default value of 15.
  #
  # Returns:
  #   The ADASYN oversampled dataset sampleADA.
  nt <- ncol(p)  # number of samples in positive data
  if (nt == 0) {
    stop ("The minority class is empty")
  } else if (nt == 1) {
    sampleADA <- kronecker(matrix(1, 1, nTarget), p)  # duplicate
  } else {
    if (k > nt-1) {
      k <- nt-1  # number of nearest neighbours can not be greater than nt-1
      warning ("The minority class instances is not enough. k is set to ", k)
    } 
    numAtt <- nrow(p)  # Feature dimension
    ratio <- FindRatio(p, n, m)  # the ratio of each positive sample need to be duplicated
    no <- round(nTarget*ratio)  # the number of each positive sample need to be duplicated
    # adjust no to make the total number of new created samples to equal to the number needed
    while (sum(no) != nTarget) {  
      # tmp <- max(no)
      ind <- which.max(no)
      diff <- nTarget - sum(no)
      if (no[ind] + diff > 0) {
        no[ind] <- no[ind] + diff
      } else {
        no[ind] <- 0
      }
    }
    # data generation
    sampleADA <- list()
    for (i in 1:length(no)) {
      if (no[i] == 0) {  # jump the positive samples which don't need to be duplicated
        next
      }
      # k-NN
      d <- rdist(t(p[, i]), t(p))  # the Euclidean distance between each positive sample and other positive data
      d[i] <-Inf  # Set d[i] to infinity manually
      # Find the k indices corresponding to the closest indices
      if (k<log(nt)) {
        minId <- list()
        for (j in 1:k) {
          # tmp <- min(d)
          id <- which.min(d)
          d[id] <-Inf
          minId <- cbind(minId, id)  # sort>=O(n*logn),so we take min: O(n).total time:O(k*n)
        } 
      }else {
        # tmp <- sort(d)
        id <- order(d)
        minId <- id[1:k]
      }
      
      rn <- floor(runif(no[i], min = 0, max = k)) + 1  # random generated No[i] elements integer vector in range 1 to k
      id <- minId[rn]
      weight <- matrix(runif(numAtt * no[i]), nrow = numAtt, ncol = no[i], byrow = TRUE)
      kro <- kronecker(matrix(1, 1, no[i]), p[, i])
      
      # for numeric attributes
      aid <- 1:numAtt
      kro[aid, ] <- kro[aid, ] + weight[aid, ]*(p[aid, unlist(id)] - kro[aid, ])
      
      sampleADA <- cbind(sampleADA, kro)
      sampleADA <- matrix(unlist(sampleADA), ncol = dim(sampleADA)[2])
    }
  }
  return(sampleADA)
}

#' Generate samples by ADASYN approach.
#' 
#' @param p minority class samples
#' @param n majority class samples
#' @param nTarget the targeted number of samples to achieve
#' @param k number of nearest neighbours in k-NN search used by the ADASYN algorithm, with the default value of 5
#' @param m seeds from the positive class in m-NN search of the ADASYN algorithm, with the default value of 15
#' @return sampleADA
#' @importFrom fields rdist 
#' @importFrom stats runif 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal

ADASYNBar <- function(p, n, nTarget, k, m) {
  # Generate samples by ADASYN.
  #
  # Args:
  #   p:       The minority class samples.
  #   n:       The majority class samples. P and N must have the same feature dimention, greater than one,
  #            with no missing values.
  #   nTarget: The targeted number of samples to achieve.
  #   k:       k-NN used in the ADASYN algorithm, with the default value of 5.
  #   m:       m-NN used in ADASYN, finding seeds from the Positive Class, with the default value of 15.
  #
  # Returns:
  #   The ADASYN oversampled dataset sampleADA.
  nt <- ncol(p)  # number of samples in p
  if (nt == 0) {
    stop ("The minority class is empty")
  } else if (nt == 1) {
    sampleADA <- kronecker(matrix(1, 1, nTarget), p)  # duplicate
  } else {
    cat("Oversampling by ADASYN: \n")
    if (k > nt-1) {
      k <- nt-1  # number of nearest neighbours can not be greater than NT-1
      warning ("The minority class instances is not enough. k is set to ", k)
    } 
    numAtt <- nrow(p)  # Feature dimension
    ratio <- FindRatio(p, n, m)  # the ratio of each positive sample need to be duplicated
    no <- round(nTarget * ratio)  # the number of each positive sample need to be duplicated
    # adjust No to make the total number of new created samples to equal to the number needed
    while (sum(no) != nTarget) {
      # tmp <- max(no)
      ind <- which.max(no)
      diff <- nTarget - sum(no)
      if (no[ind] + diff > 0) {
        no[ind] <- no[ind] + diff
      } else {
        no[ind] <- 0
      }
    }
    # data generation
    sampleADA <- list()
    pb <- txtProgressBar(min = 0, max = length(no), style = 3)  # progress bar
    
    for (i in 1:length(no)) {
      if (no[i] == 0) {  # jump the positive samples which don't need to be duplicated
        next
      }
      # k-NN
      d <- rdist(t(p[, i]), t(p))  # the Euclidean distance between each positive sample and other positive data
      d[i] <-Inf  # Set d[i] to infinity manually
      # Find the k indices corresponding to the closest indices
      if (k<log(nt)) {
        minId <- list()
        for (j in 1:k) {
          # tmp <- min(d)
          id <- which.min(d)
          d[id] <-Inf
          minId <- cbind(minId, id)  # sort>=O(n*logn),so we take min: O(n).total time:O(k*n)
        } 
      }else {
        # tmp <- sort(d)
        id <- order(d)
        minId <- id[1:k]
      }
      
      rn <- floor(runif(no[i], min = 0, max = k)) + 1  # random generated No[i] elements integer vector in range 1 to k
      id <- minId[rn]
      weight <- matrix(runif(numAtt * no[i]), nrow = numAtt, ncol = no[i], byrow = TRUE)
      kro <- kronecker(matrix(1, 1, no[i]), p[, i])
      
      # for numeric attributes
      aid <- 1:numAtt
      kro[aid, ] <- kro[aid, ] + weight[aid, ]*(p[aid, unlist(id)] - kro[aid, ])
      sampleADA <- cbind(sampleADA, kro)
      sampleADA <- matrix(unlist(sampleADA), ncol = dim(sampleADA)[2])
      
      setTxtProgressBar(pb, i)  # update progress bar
    }
    close(pb)  # close progress bar
  }
  return(sampleADA)
}

#' Generate samples by ADASYN approach.
#' 
#' @param p minority class samples
#' @param n majority class samples
#' @param nTarget the targeted number of samples to achieve
#' @param k number of the neareat neighbour in k-NN used by the ADASYN algorithm, with the default value of 5
#' @param m seeds from the positive class in m-NN used by the ADASYN algorithm, with the default value of 15
#' @return sampleADA
#' @importFrom fields rdist 
#' @importFrom stats runif 
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @keywords internal

ADASYNPara <- function(p, n, nTarget, k, m) {
  # Generate samples by ADASYN.
  #
  # Args:
  #   p:       The minority class samples.
  #   n:       The majority class samples. P and N must have the same feature dimension, greater than one,
  #            with no missing values.
  #   nTarget: The targeted number of samples to achieve.
  #   k:       k-NN used in the ADASYN algorithm, with the default value of 5.
  #   m:       m-NN used in ADASYN, finding seeds from the Positive Class, with the default value of 15.
  #
  # Returns:
  #   The ADASYN oversampled dataset sampleADA.
  nt <- ncol(p)  # NT is number of samples in P
  if (nt == 0) {
    stop ("The minority class is empty")
  } else if (nt == 1) {
    sampleADA <- kronecker(matrix(1, 1, nTarget), p)  # duplicate
  } else {
    if (k > nt-1) {
      k <- nt-1  # number of nearest neighbours can not be greater than nt-1
      warning ("The minority class instances is not enough. k is set to ", k)
    } 
    numAtt <- nrow(p)  # Feature dimension
    ratio <- FindRatioPara(p, n, m)  # the ratio of each positive sample need to be duplicated
    no <- round(nTarget*ratio)  # the number of each positive sample need to be duplicated
    # adjust No to make the total number of new created samples to equal to the number needed
    while (sum(no) != nTarget) {
      # tmp <- max(no)
      ind <- which.max(no)
      diff <- nTarget - sum(no)
      if (no[ind] + diff > 0) {
        no[ind] <- no[ind] + diff
      } else {
        no[ind] <- 0
      }
    }
    # data generation
    nlen <- length(no)  # number of positive samples    
    i <- 0
    cl <- makeCluster(detectCores(logical = FALSE) - 1)  # start parallel
    registerDoParallel(cl)
    sampleADA <- foreach(i = 1:nlen, .combine = 'cbind') %dopar% {
      if (no[i] != 0) {
        # k-NN
        d <- rdist(t(p[, i]), t(p))  # the Euclidean distance between each positive sample and other positive data
        d[i] <-Inf  # Set d[i] to infinity manually
        # Find the k indices corresponding to the closest indices
        if (k<log(nt)) {
          minId <- list()
          for (j in 1:k) {
            # tmp <- min(d)
            id <- which.min(d)
            d[id] <-Inf
            minId <- cbind(minId, id)  # sort>=O(n*logn),so we take min: O(n).total time:O(k*n)
          } 
        }else {
          # tmp <- sort(d)
          id <- order(d)
          minId <- id[1:k]
        }
        
        rn <- floor(runif(no[i], min = 0, max = k)) + 1  # random generated No[i] elements integer vector in range 1 to k
        id <- minId[rn]
        weight <- matrix(runif(numAtt * no[i]), nrow = numAtt, ncol = no[i], byrow = TRUE)
        kro <- kronecker(matrix(1, 1, no[i]), p[, i])
        
        # for numeric attributes
        aid <- 1:numAtt
        kro[aid, ] <- kro[aid, ] + weight[aid, ]*(p[aid, unlist(id)] - kro[aid, ])
        
        return(kro)
      }
    }
    stopCluster(cl)  # end parallel
  }
  return(sampleADA)
}

#' Generate samples by ADASYN approach.
#' 
#' @param p minority class samples
#' @param n majority class samples
#' @param nTarget the targeted number of samples to achieve
#' @param k number of nearest neighbours in k-NN search used by the ADASYN algorithm, with the default value of 5
#' @param m seeds from the positive class in m-NN of the ADASYN algorithm, with the default value of 15
#' @return sampleADA
#' @importFrom fields rdist 
#' @importFrom stats runif 
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal

ADASYNParaBar <- function(p, n, nTarget, k, m) {
  # Generate samples by ADASYN.
  #
  # Args:
  #   p:       The minority class samples.
  #   n:       The majority class samples. P and N must have the same feature dimention, greater than one,
  #            with no missing values.
  #   nTarget: The targeted number of samples to achieve.
  #   k:       k-NN used in the ADASYN algorithm, with the default value of 5.
  #   m:       m-NN used in ADASYN, finding seeds from the Positive Class, with the default value of 15.
  #
  # Returns:
  #   The ADASYN oversampled dataset sampleADA.
  nt <- ncol(p)  # number of samples in p
  if (nt == 0) {
    stop ("The minority class is empty")
  } else if (nt == 1) {
    sampleADA <- kronecker(matrix(1, 1, nTarget), p)  # duplicate
  } else {
    cat("Oversampling by ADASYN: \n")
    if (k > nt-1) {
      k <- nt-1  # number of nearest neighbours can not be greater than nt-1
      warning ("The minority class instances is not enough. k is set to ", k)
    } 
    
    numAtt <- nrow(p)  # Feature dimension
    ratio <- FindRatioPara(p, n, m)  # the ratio of each positive sample need to be duplicated
    no <- round(nTarget*ratio)  # the number of each positive sample need to be duplicated
    # adjust No to make the total number of new created samples to equal to the number needed
    while (sum(no) != nTarget) {
      # tmp <- max(No)
      ind <- which.max(no)
      diff <- nTarget - sum(no)
      if (no[ind] + diff > 0) {
        no[ind] <- no[ind] + diff
      } else {
        no[ind] <- 0
      }
    }
    # data generation
    nlen <- length(no)  # number of positive samples 
    i <- 0
    cl <- makeCluster(detectCores(logical = FALSE) - 1)  # start parallel
    registerDoSNOW(cl)
    pb <- txtProgressBar(min = 0, max = nlen, style = 3)  # progress bar
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    # registerDoParallel(cl)
    sampleADA <- foreach(i = 1:nlen, .combine = 'cbind', .options.snow = opts) %dopar% {
      if (no[i] != 0) {
        # k-NN
        d <- rdist(t(p[, i]), t(p))  # the Euclidean distance between each positive sample and other positive data
        d[i] <-Inf  # Set d[i] to infinity manually
        # Find the k indices corresponding to the closest indices
        if (k<log(nt)) {
          minId <- list()
          for (j in 1:k) {
            # tmp <- min(d)
            id <- which.min(d)
            d[id] <-Inf
            minId <- cbind(minId, id)  # sort>=O(n*logn),so we take min: O(n).total time:O(k*n)
          } 
        }else {
          # tmp <- sort(d)
          id <- order(d)
          minId <- id[1:k]
        }
        
        rn <- floor(runif(no[i], min = 0, max = k)) + 1  # random generated No[i] elements integer vector in range 1 to k
        id <- minId[rn]
        weight <- matrix(runif(numAtt * no[i]), nrow = numAtt, ncol = no[i], byrow = TRUE)
        kro <- kronecker(matrix(1, 1, no[i]), p[, i])
        
        # for numeric attributes
        aid <- 1:numAtt
        kro[aid, ] <- kro[aid, ] + weight[aid, ]*(p[aid, unlist(id)] - kro[aid, ])
        
        return(kro)
      }
    }
    close(pb)  # close progress bar
    stopCluster(cl)  # end parallel
  }
  return(sampleADA)
}

#' Generate samples by ESPO algorithm.
#' 
#' @param me Mean vector of positive class
#' @param v Eigen axes matrix (Each axis is a column vector)
#' @param dMod Modified Eigen Spectrum Value
#' @param p The minority class samples
#' @param n The majority class samples
#' @param r A scalar ratio specifies which level (towards the boundary) we shall push the synthetic data, 
#'          with the default value 1
#' @param m Scalar specifies the reliable portion of the eigen spectrum
#' @param numToGen The number of samples to be generated
#' @return sampleESPO
#' @importFrom fields rdist 
#' @importFrom MASS mvrnorm
#' @keywords internal

ESPO <- function(me, v, dMod, p, n, r, m, numToGen) {
  # Generate samples by ESPO.
  #
  # Args:
  #   me:       Mean vector of positive class.
  #   v:        Eigen axes matrix (Each axis is a column vector).
  #   dMod:        Modified Eigen Spectrum Value.
  #   p:        The minority class samples.
  #   n:        The majority class samples. P and N must have the same feature dimension, greater than one,
  #             with no missing values.
  #   r:        A scalar ratio specifies which level (towards the boundary) we shall push the synthetic data,
  #             with a default value of 1.
  #   m:        Scalar specifies the reliable portion of the eigen spectrum.
  #   numToGen: The number of samples to be generated.
  #
  # Returns:
  #   The ESPO oversampled dataset sampleESPO.
  rn <- m  # reliable portion of the eigen spectrum
  un <- length(me) - m  # unreliable portion of the eigen spectrum
  
  muR <- matrix(0, 1, rn)  # mean
  sigmaR <- diag(1, rn)  # standard deviation
  
  muU <- matrix(0, 1, un)  # mean
  sigmaU <- diag(1, un)  # standard deviation
  
  sampGen <- matrix(0, numToGen * r, length(me))  # total samples needed 
  sampSel <- matrix(0, numToGen, length(me))  # total samples which should be kept
  prob <- matrix(0, numToGen*r, 1)  # probability of each sample to be kept
  
  cnt <- 0
  dd <- sqrt(dMod)  # square root of modified eigen spectrum value
  
  #  generate new positive data sequence until accepted, using Euclidean distance checking
  while (cnt < r * numToGen) {
    aR <- mvrnorm(1, muR, sigmaR)  # generate random vectors from the multivariate normal distribution
    dens <- exp(-0.5*sum(aR^2) - length(aR)*log(2*pi)/2)  # the density of the multivariate normal distribution
    
    if (un > 0) {
      aU <- mvrnorm(1, muU, sigmaU)
      a <- c(aR, aU)*dd  # the vector in eigen transformed domain
    } else {
      a <- aR*dd
    }
    x <- a %*% t(v) + me  # the modified generated vector
    
    pDist <- rdist(x, p)  # the Euclidean distance between x and positive data
    nDist <- rdist(x, n)  # the Euclidean distance between x and negative data
    
    val <- min(nDist)  # the value of the smallest element in the Euclidean distances between x and negative data
    ind <- which.min(nDist)  # the index of the smallest element in the Euclidean distances between x and negative data
    
    # check whether to keep the generated vector using the Euclidean distance between negative and positive samples
    if (min(pDist) < val) {
      ppDist <- rdist(t(n[ind, ]), p)
      if (val >= min(ppDist) && val <= max(ppDist)) {
        cnt <- cnt + 1
        sampGen[cnt, ] <- x
        prob[cnt, 1] <- dens
      }
    }
  }
  
  # Draw samples from a multivariate normal distribution and discard based on the ratio r
  for (i in 1:numToGen) {
    # tmp <- min(Prob)
    ind <- which.min(prob)
    prob[ind] <- Inf
    sampSel[i, ] <- sampGen[ind, ]
  }
  
  # form new dataset
  sampleESPO <- rbind(sampSel, p)
  return(sampleESPO)
}

#' Generate samples by ESPO algorithm.
#' 
#' @param me Mean vector of positive class
#' @param v Eigen axes matrix (Each axis is a column vector)
#' @param dMod Modified Eigen Spectrum Value
#' @param p The minority class samples
#' @param n The majority class samples
#' @param r scalar specifies which level (towards the boundary) we shall push the synthetic data, 
#'          with the default value =  1
#' @param m Scalar specifies the reliable portion of the eigen spectrum
#' @param numToGen The number of samples to be generated
#' @return sampleESPO
#' @importFrom fields rdist 
#' @importFrom MASS mvrnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal

ESPOBar <- function(me, v, dMod, p, n, r, m, numToGen) {
  # Generate samples by ESPO.
  #
  # Args:
  #   me:       Mean vector of positive class.
  #   v:        Eigen axes matrix (Each axis is a column vector).
  #   dMod:        Modified Eigen Spectrum Value.
  #   p:        The minority class samples.
  #   n:        The majority class samples. P and N must have the same feature dimention, greater than one,
  #             with no missing values.
  #   r:        Scalar ratio specifies which level (towards the boundary) we shall push the synthetic data,
  #             with the default value 1.
  #   m:        Scalar specifies the reliable portion of the eigen spectrum.
  #   NumToGen: The number of samples to be generated.
  #
  # Returns:
  #   The ESPO oversampled dataset sample_espo.
  rn <- m  # reliable portion of the eigen spectrum
  un <- length(me) - m  # unreliable portion of the eigen spectrum
  
  muR <- matrix(0, 1, rn)  # mean
  sigmaR <- diag(1, rn)  # standard deviation
  
  muU <- matrix(0, 1, un)  # mean
  sigmaU <- diag(1, un)  # standard deviation
  
  sampGen <- matrix(0, numToGen * r, length(me))  # total samples needed be created
  sampSel <- matrix(0, numToGen, length(me))  # total samples that should be kept
  prob <- matrix(0, numToGen*r, 1)  # probability of each sample to be kept
  
  dd <- sqrt(dMod)  # square root of modified eigen spectrum value
  cat("Oversampling by ESPO: \n")
  nGener <- r * numToGen  # number of total samples needed be created
  pb <- txtProgressBar(min = 0, max = nGener, style = 3, char = "-")  # progress bar
  for (cnt in 1:nGener) {
    flag <-TRUE
    #  genetare new positive data sequence until accepted upon Euclidean distance checking
    while (flag) {
      aR <- mvrnorm(1, muR, sigmaR)  # generate random vectors from the multivariate normal distribution
      dens <- exp(-0.5*sum(aR^2) - length(aR)*log(2*pi)/2)  # the density of the multivariate normal distribution
      
      if (un > 0) {
        aU <- mvrnorm(1, muU, sigmaU)
        a <- c(aR, aU)*dd  # The vector in Eigen transformed domain
      } else {
        a <- aR*dd
      }
      x <- a %*% t(v) + me  # the modified generated vector
      
      pDist <- rdist(x, p)  # the Euclidean distance between x and positive data
      nDist <- rdist(x, n)  # the Euclidean distance between x and negative data
      
      val <- min(nDist)  # the value of the smallest element in the Euclidean distances between x and negative data
      ind <- which.min(nDist)  # the index of the smallest element in the Euclidean distances between x and negative data
      
      # check whether to keep the generated vector using the Euclidean distance between negative and positive data
      if (min(pDist) < val) {
        ppDist <- rdist(t(n[ind, ]), p)
        if (val >= min(ppDist) && val <= max(ppDist)) {
          flag <- FALSE
          sampGen[cnt, ] <- x
          prob[cnt, 1] <- dens
        }
      }
    }
    setTxtProgressBar(pb, cnt)
  }
  close(pb)  # end progress bar
  
  # Draw samples from the multivariate normal distribution, and discard based on the ratio r
  for (i in 1:numToGen) {
    # tmp <- min(Prob)
    ind <- which.min(prob)
    prob[ind] <- Inf
    sampSel[i, ] <- sampGen[ind, ]
  }
  
  # form new dataset
  sampleESPO <- rbind(sampSel, p)
  return(sampleESPO)
}

#' Generate samples by ESPO algorithm.
#' 
#' @param me Mean vector of positive class
#' @param v Eigen axes matrix (Each axis is a column vector)
#' @param dMod Modified Eigen Spectrum Value
#' @param p The minority class samples
#' @param n The majority class samples
#' @param r A scalar ratio specifies which level (towards the boundary) we shall push the synthetic data, 
#'          with the default value = 1
#' @param m Scalar specifies the reliable portion of the eigen spectrum
#' @param numToGen The number of samples to be generated
#' @return sampleESPO
#' @importFrom fields rdist 
#' @importFrom MASS mvrnorm
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @keywords internal

ESPOPara <- function(me, v, dMod, p, n, r, m, numToGen) {
  # Generate samples by ESPO.
  #
  # Args:
  #   me:       Mean vector of positive class.
  #   v:        Eigen axes matrix (Each axis is a column vector).
  #   dMod:        Modified Eigen Spectrum Value.
  #   p:        The minority class samples.
  #   n:        The majority class samples. p and n must have the same feature dimension, greater than one,
  #             with no missing values.
  #   r:        Scalar ratio specifies which level (towards the boundary) we shall push the synthetic data,
  #             with the default value 1.
  #   m:        Scalar specifies the reliable portion of the eigen spectrum.
  #   numToGen: The number of samples to be generated.
  #
  # Returns:
  #   The ESPO oversampled dataset sampleESPO.
  rn <- m  # reliable portion of the eigen spectrum
  un <- length(me) - m  # unreliable portion of the eigen spectrum
  
  muR <- matrix(0, 1, rn)  # mean
  sigmaR <- diag(1, rn)  # standard deviation
  
  muU <- matrix(0, 1, un)  # mean
  sigmaU <- diag(1, un)  # standard deviation
  
  sampSel <- matrix(0, numToGen, length(me))  # total samples would be kept
  
  dd <- sqrt(dMod)  # square root of modified eigen spectrum value
  
  nGener <- r * numToGen  # number of total samples needed be created
  
  cl <- makeCluster(detectCores(logical = FALSE) - 1)  # start parallel
  registerDoParallel(cl)
  seq <- foreach(cnt = 1:nGener, .combine = 'rbind') %dopar% {
    flag <- TRUE
    #  genetare new positive data sequence until accepted using Euclidean distance checking
    while (flag) {
      aR <- mvrnorm(1, muR, sigmaR)  # generate random vectors from the multivariate normal distribution
      dens <- exp(-0.5*sum(aR^2) - length(aR)*log(2*pi)/2)  # the density of the multivariate normal distribution
      
      if (un > 0) {
        aU <- mvrnorm(1, muU, sigmaU)
        a <- c(aR, aU)*dd  # The vector in Eigen transformed domain
      } else {
        a <- aR*dd
      }
      x <- a %*% t(v) + me  # the modified generated vector
      
      pDist <- rdist(x, p)  # the Euclidean distance between x and positive data
      nDist <- rdist(x, n)  # the Euclidean distance between x and negative data
      
      val <- min(nDist)  # the value of the smallest element in the Euclidean distances between x and negative samples
      ind <- which.min(nDist)  # the index of the smallest element in the Euclidean distances between x and negative samples
      
      # check wheter to keep the generated vector using the Euclidean distance between negative and positive data
      if (min(pDist) < val) {
        ppDist <- rdist(t(n[ind, ]), p)
        if (val >= min(ppDist) && val <= max(ppDist)) {
          res <- list("x" = x, "dens" = dens)
          flag <- FALSE
          return(res)
        }
      }
    }
  }
  stopCluster(cl)  # end parallel
  
  sampGen <- matrix(unlist(seq[, 1]), ncol = length(me), byrow = TRUE)
  prob <- matrix(unlist(seq[, 2]), ncol = 1, byrow = TRUE)
  # Draw samples from a multivariate normal distribution and discard based on the ratio r
  for (i in 1:numToGen) {
    # tmp <- min(Prob)
    ind <- which.min(prob)
    prob[ind] <- Inf
    sampSel[i, ] <- sampGen[ind, ]
  }
  
  # form new dataset
  sampleESPO <- rbind(sampSel, p)
  return(sampleESPO)
}

#' Generate samples by ESPO algorithm.
#' 
#' @param me Mean vector of positive class
#' @param v Eigen axes matrix (Each axis is a column vector)
#' @param dMod Modified Eigen Spectrum Value
#' @param p The minority class samples
#' @param n The majority class samples
#' @param r A scalar specifying which level (towards the boundary) we shall push the synethetic data, 
#'          with the default value = 1
#' @param m Scalar specifies the reliable portion of the eigen spectrum
#' @param numToGen The number of samples to be generated
#' @return sampleESPO
#' @importFrom fields rdist 
#' @importFrom MASS mvrnorm
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal

ESPOParaBar <- function(me, v, dMod, p, n, r, m, numToGen) {
  # Generate samples by ESPO.
  #
  # Args:
  #   me:       Mean vector of positive class.
  #   v:        Eigen axes matrix (Each axis is a column vector).
  #   dMod:        Modified Eigen Spectrum Value.
  #   p:        The minority class samples.
  #   n:        The majority class samples. P and N must have the same feature dimention, greater than one,
  #             with no missing values.
  #   r:        A scalar ratio specifying which level (towards the boundary) we shall push the synthetic data,
  #             with the default value = 1.
  #   m:        Scalar specifies the reliable portion of the eigen spectrum.
  #   numToGen: The number of samples to be generated.
  #
  # Returns:
  #   The ESPO oversampled dataset sampleESPO.
  rn <- m  # reliable portion of the eigen spectrum
  un <- length(me) - m  # unreliable portion of the eigen spectrum
  
  muR <- matrix(0, 1, rn)  # mean
  sigmaR <- diag(1, rn)  # standard deviation
  
  muU <- matrix(0, 1, un)  # mean
  sigmaU <- diag(1, un)  # standard deviation
  
  sampSel <- matrix(0, numToGen, length(me))  # total samples would be kept
  
  dd <- sqrt(dMod)  # square root of modified eigen spectrum value
  
  nGener <- r * numToGen  # number of total samples needed be created
  cat("Oversampling by ESPO: \n")
  cl <- makeCluster(detectCores(logical = FALSE) - 1)  # start parallel
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 0, max = nGener, style = 3, char = "-")  # progress bar
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # registerDoParallel(cl, cores = cores)
  seq <- foreach(cnt = 1:nGener, .combine = 'rbind', .options.snow = opts) %dopar% {
    flag <- TRUE
    #  generate new positive data sequence until accepted, using Euclidean distance checking
    while (flag) {
      aR <- mvrnorm(1, muR, sigmaR)  # draw random vectors from the multivariate normal distribution
      dens <- exp(-0.5*sum(aR^2) - length(aR)*log(2*pi)/2)  # the density of the multivariate normal distribution
      
      if (un > 0) {
        aU <- mvrnorm(1, muU, sigmaU)
        a <- c(aR, aU)*dd  # The vector in Eigen transformed domain
      } else {
        a <- aR*dd
      }
      x <- a %*% t(v) + me  # the modified generated vector
      
      pDist <- rdist(x, p)  # the Euclidean distance between x and positive data
      nDist <- rdist(x, n)  # the Euclidean distance between x and negative data
      
      val <- min(nDist)  # the value of the smallest element of the Euclidean distances between x and the negative samples
      ind <- which.min(nDist)  # the index of the smallest element of the Euclidean distances between x and the negative samples
      
      # check whether to keep the generated vector using the Euclidean distance between negative and positive data
      if (min(pDist) < val) {
        ppDist <- rdist(t(n[ind, ]), p)
        if (val >= min(ppDist) && val <= max(ppDist)) {
          res <- list("x" = x, "dens" = dens)
          flag <- FALSE
          return(res)
        }
      }
    }
  }
  close(pb)  # end progress bar
  stopCluster(cl)  # end parallel
  
  sampGen <- matrix(unlist(seq[, 1]), ncol = length(me), byrow = TRUE)
  prob <- matrix(unlist(seq[, 2]), ncol = 1, byrow = TRUE)
  # Draw samples from a multivariate normal distribution and discard based on the ratio r
  for (i in 1:numToGen) {
    # tmp <- min(Prob)
    ind <- which.min(prob)
    prob[ind] <- Inf
    sampSel[i, ] <- sampGen[ind, ]
  }
  
  # form new dataset
  sampleESPO <- rbind(sampSel, p)
  return(sampleESPO)
}

#' Find the distribution of the positive data.
#' 
#' @param p minority class samples
#' @param n majority class samples
#' @param m For m-NN used in ADASYN, finding seeds from the Positive Class, with the default value = 15
#' @return ratio
#' @importFrom fields rdist 
#' @keywords internal

FindRatio <- function(p, n, m) {
  # Find the distribution of the positive data to determine the ratio of each positive data to be generated.
  #
  # Args:
  #   p: minority class samples.
  #   n: majority class samples. P and N must have the same feature dimension, greater than one,
  #      with no missing values.
  #   m: m-NN used in finding seeds from the Positive Class, with the default value = 15
  #
  # Returns:
  #   ratio: the ratio of each positive sample need to be duplicated.
  c <- cbind(p, n)  # Note that here columns of P and N are the samples
  poscnt <- ncol(p)  # Number of positive records
  
  # The ratio of each positive sample needed to be duplicated
  ratio <- matrix(0, poscnt, 1)
  for (i in 1:poscnt) {
    d <- rdist(t(p[, i]), t(c))  # the Euclidean distance between each positive sample and full data
    d[i] <- Inf
    # find the indices of m number smallest elements using the Euclidean distance between each positive sample and full data
    minId <- matrix(0, m, 1)
    for (j in 1:m) {
      # tmp <- min(d)
      id <- which.min(d)
      d[id] <- Inf
      minId[j] <- id  # sort>=O(n*logn),so we take min: O(n).total time:O(k*n)
    }
    ind <- which(minId > poscnt)  # find all negative samples from this m number closest elements
    ratio[i] <- length(ind)
  }
  ratio <- ratio/sum(ratio)
  return(ratio)
}

#' Find the distribution of the positive data.
#' 
#' @param p minority class samples
#' @param n majority class samples
#' @param m m-NN used in ADASYN, finding seeds from the Positive Class, with the default value = 15
#' @return ratio
#' @importFrom fields rdist 
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @keywords internal

FindRatioPara <- function(p, n, m) {
  # Find the distribution of the positive data to determine the ratio of each positive data to be generated.
  #
  # Args:
  #   p: minority class samples.
  #   n: majority class samples. P and N must have the same feature dimention, greater than one,
  #      with no missing values.
  #   m: m-NN used in finding seeds from the Positive Class, with the default value = 15
  #
  # Returns:
  #   ratio: the ratio of each positive sample needed to be duplicated.
  c <- cbind(p, n)  # Note that here columns of p and n are the samples
  poscnt <- ncol(p)  # Number of positive records
  i <- 0
  
  #   The ratio of each positive sample needed to be duplicated
  cl <- makeCluster(detectCores(logical = FALSE) - 1)  # start parallel
  registerDoParallel(cl)
  seq <- foreach(i = 1:poscnt, .combine = 'rbind') %dopar% {
    d <- rdist(t(p[, i]), t(c))  # the Euclidean distance between each positive sample and full data
    d[i] <- Inf
    # find the indices of m number smallest elements from the Euclidean distance between each positive sample and full data
    minId <- matrix(0, m, 1)
    for (j in 1:m) {
      # tmp <- min(d)
      id <- which.min(d)
      d[id] <- Inf
      minId[j] <- id  # sort>=O(n*logn),so we take min: O(n).total time:O(k*n)
    }
    ind <- which(minId > poscnt)  # find all negative samples from the m number closest elements
    return(length(ind))
  }  
  stopCluster(cl)  # end parallel  
  ratio <- matrix(unlist(seq), ncol = 1, byrow = TRUE)
  ratio <- ratio/sum(ratio)
  
  return(ratio)
}

