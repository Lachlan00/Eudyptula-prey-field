# A customised version of simData in momentuHMM to add shapefile constraints
# to simulations (i.e. keep penguin simulations off land).
# Line 1059 is k loop for point generation
library(momentuHMM)
library(CircStats)

simData_custom <- function (nbAnimals = 1, nbStates = 2, dist, Par, beta = NULL, 
          delta = NULL, formula = ~1, formulaDelta = NULL, mixtures = 1, 
          formulaPi = NULL, covs = NULL, nbCovs = 0, spatialCovs = NULL, 
          zeroInflation = NULL, oneInflation = NULL, circularAngleMean = NULL, 
          centers = NULL, centroids = NULL, angleCovs = NULL, obsPerAnimal = c(500,1500), 
          initialPosition = c(0, 0), DM = NULL, cons = NULL, 
          userBounds = NULL, workBounds = NULL, workcons = NULL, betaRef = NULL, 
          mvnCoords = NULL, stateNames = NULL, model = NULL, states = FALSE, 
          retrySims = 0, lambda = NULL, errorEllipse = NULL,
          # Additions
          landShp=NULL) 
{
  #### Hack - Load internal function #########################################
  mvndists <- momentuHMM:::mvndists
  is.miHMM <- momentuHMM:::is.miHMM
  is.miSum <- momentuHMM:::is.miSum
  delta_bc <- momentuHMM:::delta_bc
  stepdists <- momentuHMM:::stepdists
  angledists <- momentuHMM:::angledists
  nw2w <- momentuHMM:::nw2w
  n2w <- momentuHMM:::n2w
  w2wn <- momentuHMM:::w2wn
  w2n <- momentuHMM:::w2n
  parDef <- momentuHMM:::parDef
  getCovNames <- momentuHMM:::getCovNames
  zeroInflationdists <- momentuHMM:::zeroInflationdists
  oneInflationdists <- momentuHMM:::oneInflationdists
  rwdists <- momentuHMM:::rwdists
  checkInputs <- momentuHMM:::checkInputs
  printMessage <- momentuHMM:::printMessage
  newFormulas <- momentuHMM:::newFormulas
  getWorkBounds <- momentuHMM:::getWorkBounds
  getDM <- momentuHMM:::getDM
  get_ncmean <- momentuHMM:::get_ncmean
  momentuHMMData <- momentuHMM:::momentuHMMData
  ############################################################################
  # Get shp projection
  landShpProj <- proj4string(landShp)
  ############################################################################
  if (!is.null(cons) & is.null(model)) 
    warning("cons argument is deprecated in momentuHMM >= 1.4.0. Please use workBounds instead.")
  if (!is.null(workcons) & is.null(model)) 
    warning("workcons argument is deprecated in momentuHMM >= 1.4.0. Please use workBounds instead.")
  if (!is.null(workBounds) & (!is.null(cons) | !is.null(workcons)) & 
      is.null(model)) 
    stop("workBounds cannot be specified when using deprecated arguments cons or workcons; either workBounds or both cons and workcons must be NULL")
  if (!is.null(model)) {
    if (is.miHMM(model)) {
      model <- model$miSum
    }
    if (inherits(model, "momentuHierHMM") | inherits(model, 
                                                     "hierarchical")) 
      stop("model can not be a 'momentuHierHMM' or 'hierarchical' object; use simHierData instead")
    model <- delta_bc(model)
    nbStates <- length(model$stateNames)
    dist <- model$conditions$dist
    distnames <- names(dist)
    if (is.miSum(model)) {
      model <- formatmiSum(model)
      if (!is.null(model$mle$beta)) 
        model$conditions$workBounds$beta <- matrix(c(-Inf, 
                                                     Inf), length(model$mle$beta), 2, byrow = TRUE)
      if (!is.null(model$Par$beta$pi$est)) 
        model$conditions$workBounds$pi <- matrix(c(-Inf, 
                                                   Inf), length(model$Par$beta$pi$est), 2, byrow = TRUE)
      if (!is.null(model$Par$beta$delta$est)) 
        model$conditions$workBounds$delta <- matrix(c(-Inf, 
                                                      Inf), length(model$Par$beta$delta$est), 2, 
                                                    byrow = TRUE)
      if (!is.null(model$mle$g0)) 
        model$conditions$workBounds$g0 <- matrix(c(-Inf, 
                                                   Inf), length(model$mle$g0), 2, byrow = TRUE)
      if (!is.null(model$mle$theta)) 
        model$conditions$workBounds$theta <- matrix(c(-Inf, 
                                                      Inf), length(model$mle$theta), 2, byrow = TRUE)
    }
    else {
      if (!inherits(model, "momentuHMM")) 
        stop("model must be a 'momentuHMM' object")
    }
    userBounds <- model$conditions$bounds
    workBounds <- model$conditions$workBounds
    mvnCoords <- model$conditions$mvnCoords
    stateNames <- model$stateNames
    estAngleMean <- model$conditions$estAngleMean
    circularAngleMean <- model$conditions$circularAngleMean
    DM <- model$conditions$DM
    cons <- model$conditions$cons
    workcons <- model$conditions$workcons
    betaRef <- model$conditions$betaRef
    zeroInflation <- model$conditions$zeroInflation
    oneInflation <- model$conditions$oneInflation
    formula <- model$conditions$formula
    if (is.null(model$condition$formulaDelta)) {
      formulaDelta <- formDelta <- ~1
    }
    else formulaDelta <- formDelta <- model$condition$formulaDelta
    mixtures <- model$conditions$mixtures
    if (is.null(model$condition$formulaPi)) {
      formulaPi <- formPi <- ~1
    }
    else formulaPi <- formPi <- model$condition$formulaPi
    Par <- model$mle[distnames]
    parCount <- lapply(model$conditions$fullDM, ncol)
    for (i in distnames[!unlist(lapply(circularAngleMean, 
                                       isFALSE))]) {
      parCount[[i]] <- length(unique(gsub("cos", "", gsub("sin", 
                                                          "", colnames(model$conditions$fullDM[[i]])))))
    }
    parindex <- c(0, cumsum(unlist(parCount))[-length(model$conditions$fullDM)])
    names(parindex) <- distnames
    for (i in distnames) {
      if (!is.null(DM[[i]])) {
        Par[[i]] <- model$mod$estimate[parindex[[i]] + 
                                         1:parCount[[i]]]
        if (!isFALSE(circularAngleMean[[i]])) {
          names(Par[[i]]) <- unique(gsub("cos", "", 
                                         gsub("sin", "", colnames(model$conditions$fullDM[[i]]))))
        }
        else names(Par[[i]]) <- colnames(model$conditions$fullDM[[i]])
      }
    }
    for (i in distnames[which(dist %in% angledists)]) {
      if (!estAngleMean[[i]]) {
        estAngleMean[[i]] <- TRUE
        userBounds[[i]] <- rbind(matrix(rep(c(-pi, pi), 
                                            nbStates), nbStates, 2, byrow = TRUE), userBounds[[i]])
        workBounds[[i]] <- rbind(matrix(rep(c(-Inf, 
                                              Inf), nbStates), nbStates, 2, byrow = TRUE), 
                                 workBounds[[i]])
        cons[[i]] <- c(rep(1, nbStates), cons[[i]])
        workcons[[i]] <- c(rep(0, nbStates), workcons[[i]])
        if (!is.null(DM[[i]])) {
          Par[[i]] <- c(rep(0, nbStates), Par[[i]])
          if (is.list(DM[[i]])) {
            DM[[i]]$mean <- ~1
          }
          else {
            tmpDM <- matrix(0, nrow(DM[[i]]) + nbStates, 
                            ncol(DM[[i]]) + nbStates)
            tmpDM[nbStates + 1:nrow(DM[[i]]), nbStates + 
                    1:ncol(DM[[i]])] <- DM[[i]]
            diag(tmpDM)[1:nbStates] <- 1
            DM[[i]] <- tmpDM
          }
        }
      }
    }
    if (!is.null(model$conditions$recharge)) {
      g0 <- nw2w(model$mle$g0, model$conditions$workBounds$g0)
      theta <- nw2w(model$mle$theta, model$conditions$workBounds$theta)
    }
    else g0 <- theta <- NULL
    nbCovsDelta <- ncol(model$covsDelta) - 1
    if (nbStates > 1) {
      beta <- nw2w(model$mle$beta, model$conditions$workBounds$beta)
      foo <- length(model$mod$estimate) - length(g0) - 
        length(theta) - (nbCovsDelta + 1) * (nbStates - 
                                               1) * mixtures + 1
      delta <- matrix(model$mod$estimate[foo:(length(model$mod$estimate) - 
                                                length(g0) - length(theta))], nrow = (nbCovsDelta + 
                                                                                        1) * mixtures)
    }
    else {
      formulaDelta <- NULL
    }
    if (mixtures > 1) {
      nbCovsPi <- ncol(model$covsPi) - 1
      foo <- length(model$mod$estimate) - length(g0) - 
        length(theta) - (nbCovsDelta + 1) * (nbStates - 
                                               1) * mixtures - (nbCovsPi + 1) * (mixtures - 
                                                                                   1) + 1:((nbCovsPi + 1) * (mixtures - 1))
      pie <- matrix(model$mod$estimate[foo], nrow = nbCovsPi + 
                      1, ncol = mixtures - 1)
    }
    else {
      pie <- NULL
      nbCovsPi <- 0
    }
    beta <- list(beta = beta, pi = pie, g0 = g0, theta = theta)
    Par <- lapply(Par, function(x) c(t(x)))
    if (states) 
      model$data$states <- NULL
    if (is.null(covs)) {
      p <- parDef(lapply(dist, function(x) gsub("Consensus", 
                                                "", x)), nbStates, estAngleMean, zeroInflation, 
                  oneInflation, DM, userBounds)
      covNames <- character()
      for (i in distnames) {
        covNames <- c(covNames, getCovNames(model, p, 
                                            i)$DMterms)
      }
      if (!is.null(model$rawCovs)) {
        covNames <- c(colnames(model$rawCovs), covNames)
      }
      covNames <- c(covNames, colnames(model$covsPi)[which(colnames(model$covsPi) != 
                                                             "(Intercept)")], colnames(model$covsDelta)[which(colnames(model$covsDelta) != 
                                                                                                                "(Intercept)")])
      covsCol <- unique(covNames)
      factorterms <- names(model$data)[unlist(lapply(model$data, 
                                                     is.factor))]
      factorcovs <- paste0(rep(factorterms, times = unlist(lapply(model$data[factorterms], 
                                                                  nlevels))), unlist(lapply(model$data[factorterms], 
                                                                                            levels)))
      if (length(covsCol)) {
        for (jj in 1:length(covsCol)) {
          cov <- covsCol[jj]
          form <- formula(paste("~", cov))
          varform <- all.vars(form)
          if (any(varform %in% factorcovs) && !all(varform %in% 
                                                   factorterms)) {
            factorvar <- factorcovs %in% (varform[!(varform %in% 
                                                      factorterms)])
            covsCol[jj] <- rep(factorterms, times = unlist(lapply(model$data[factorterms], 
                                                                  nlevels)))[which(factorvar)]
          }
        }
      }
      covsCol <- unique(covsCol)
      covsCol <- covsCol[!(covsCol %in% "ID")]
      if (length(covsCol)) 
        covs <- model$data[covsCol]
      if (!is.null(mvnCoords) && dist[[mvnCoords]] %in% 
          rwdists) {
        covs[[paste0(mvnCoords, ".x_tm1")]] <- NULL
        covs[[paste0(mvnCoords, ".y_tm1")]] <- NULL
        if (dist[[mvnCoords]] == "rw_mvnorm3") 
          covs[[paste0(mvnCoords, ".z_tm1")]] <- NULL
        if (!ncol(covs)) 
          covs <- NULL
      }
    }
  } 
  else {
    if (!is.list(dist) | is.null(names(dist))) 
      stop("'dist' must be a named list")
    if (!is.list(Par) | is.null(names(Par))) 
      stop("'Par' must be a named list")
    distnames <- names(dist)
    if (!all(distnames %in% names(Par))) 
      stop(distnames[which(!(distnames %in% names(Par)))], 
           " is missing in 'Par'")
    Par <- Par[distnames]
    if (is.null(formulaPi)) {
      formPi <- ~1
    }
    else formPi <- formulaPi
    if (is.null(formulaDelta)) {
      formDelta <- ~1
    }
    else formDelta <- formulaDelta
    mHind <- (is.null(DM) & is.null(userBounds) & is.null(workBounds) & 
                is.null(spatialCovs) & is.null(centers) & is.null(centroids) & 
                ("step" %in% names(dist)) & all(unlist(initialPosition) == 
                                                  c(0, 0)) & is.null(lambda) & is.null(errorEllipse) & 
                !is.list(obsPerAnimal) & is.null(covs) & !nbCovs & 
                !length(attr(terms.formula(formula), "term.labels")) & 
                !length(attr(terms.formula(formDelta), "term.labels")) & 
                is.null(delta) & is.null(betaRef) & is.null(mvnCoords) & 
                mixtures == 1)
    if (all(names(dist) %in% c("step", "angle")) & all(unlist(dist) %in% 
                                                       moveHMMdists) & mHind) {
      zi <- FALSE
      if (!is.null(zeroInflation$step)) 
        zi <- zeroInflation$step
      if (is.null(dist$angle)) 
        dist$angle <- "none"
      data <- moveHMM::simData(nbAnimals, nbStates, dist$step, 
                               dist$angle, Par$step, Par$angle, beta, covs, 
                               nbCovs, zi, obsPerAnimal, model, states)
      attr(data, "class") <- "data.frame"
      data$ID <- as.factor(data$ID)
      return(momentuHMMData(data))
    }
  }
  if (nbAnimals < 1) 
    stop("nbAnimals should be at least 1.")
  if (nbStates < 1) 
    stop("nbStates should be at least 1.")
  if (is.null(zeroInflation)) {
    zeroInflation <- vector("list", length(distnames))
    names(zeroInflation) <- distnames
    for (i in distnames) {
      zeroInflation[[i]] <- FALSE
    }
  } 
  else {
    if (!is.list(zeroInflation) | is.null(names(zeroInflation))) 
      stop("'zeroInflation' must be a named list")
    for (i in distnames) {
      if (is.null(zeroInflation[[i]])) 
        zeroInflation[[i]] <- FALSE
    }
  }
  if (is.null(oneInflation)) {
    oneInflation <- vector("list", length(distnames))
    names(oneInflation) <- distnames
    for (i in distnames) {
      oneInflation[[i]] <- FALSE
    }
  } 
  else {
    if (!is.list(oneInflation) | is.null(names(oneInflation))) 
      stop("'oneInflation' must be a named list")
    for (i in distnames) {
      if (is.null(oneInflation[[i]])) 
        oneInflation[[i]] <- FALSE
    }
  }
  if (!all(unlist(lapply(zeroInflation, is.logical)))) 
    stop("zeroInflation must be a list of logical objects")
  if (!all(unlist(lapply(oneInflation, is.logical)))) 
    stop("oneInflation must be a list of logical objects")
  for (i in distnames) {
    if (!(dist[[i]] %in% zeroInflationdists) & zeroInflation[[i]]) 
      stop(dist[[i]], " distribution cannot be zero inflated")
    if (!(dist[[i]] %in% oneInflationdists) & oneInflation[[i]]) 
      stop(dist[[i]], " distribution cannot be one inflated")
  }
  if (!is.null(mvnCoords)) {
    if (length(mvnCoords) > 1 | !is.character(mvnCoords)) 
      stop("mvnCoords must be a character string")
    if (!(mvnCoords %in% distnames)) 
      stop("mvnCoords not found. Permitted values are: ", 
           paste0(distnames, collapse = ", "))
    if (!(dist[[mvnCoords]] %in% mvndists)) 
      stop("mvnCoords must correspond to a multivariate normal data stream")
    if (any(c("step", "angle") %in% distnames)) 
      stop("step and angle distributions cannot be specified when ", 
           mvnCoords, " ~ ", dist[[mvnCoords]])
    if (mvnCoords %in% c("x", "y", "z")) 
      stop("'x', 'y', and 'z' are reserved and cannot be used for mvnCoords data stream names")
    if (sum(unlist(dist) %in% rwdists) > 1) 
      stop("sorry, simData currently only supports a single multivariate normal random walk distribution (and it must correspond to location data)")
  }
  else if (any(unlist(dist) %in% rwdists)) 
    stop("sorry, simData currently only supports random walk distributions for multivariate location data identified through the mvnCoords argument")
  if (any(unlist(dist) == "rw_norm")) 
    stop("sorry, 'rw_norm' is not currently supported by simData")
  estAngleMean <- vector("list", length(distnames))
  names(estAngleMean) <- distnames
  for (i in distnames) {
    if (dist[[i]] %in% angledists) 
      estAngleMean[[i]] <- TRUE
    else estAngleMean[[i]] <- FALSE
  }
  inputs <- checkInputs(nbStates, dist, Par, estAngleMean, 
                        circularAngleMean, zeroInflation, oneInflation, DM, 
                        userBounds, cons, workcons, stateNames, checkInflation = TRUE)
  p <- inputs$p
  parSize <- p$parSize
  bounds <- p$bounds
  Fun <- lapply(inputs$dist, function(x) paste("r", x, sep = ""))
  spatialcovnames <- NULL
  if (!is.null(spatialCovs)) {
    if (!is.list(spatialCovs)) 
      stop("spatialCovs must be a list")
    spatialcovnames <- names(spatialCovs)
    if (is.null(spatialcovnames)) 
      stop("spatialCovs must be a named list")
    nbSpatialCovs <- length(spatialcovnames)
    if (is.null(mvnCoords)) {
      if (!("step" %in% distnames)) 
        stop("spatialCovs can only be included when 'step' distribution is specified")
      else if (!(inputs$dist[["step"]] %in% stepdists)) 
        stop("spatialCovs can only be included when valid 'step' distributions are specified")
    }
    for (j in 1:nbSpatialCovs) {
      if (!inherits(spatialCovs[[j]], c("RasterLayer", 
                                        "RasterBrick", "RasterStack"))) 
        stop("spatialCovs$", spatialcovnames[j], " must be of class 'RasterLayer', 'RasterStack', or 'RasterBrick'")
      if (any(is.na(raster::getValues(spatialCovs[[j]])))) 
        stop("missing values are not permitted in spatialCovs$", 
             spatialcovnames[j])
      if (inherits(spatialCovs[[j]], c("RasterBrick", 
                                       "RasterStack"))) {
        if (is.null(raster::getZ(spatialCovs[[j]]))) 
          stop("spatialCovs$", spatialcovnames[j], " is a raster stack or brick that must have set z values (see ?raster::setZ)")
        else if (!(names(attributes(spatialCovs[[j]])$z) %in% 
                   names(covs))) {
          if (!is.null(model)) 
            covs[[names(attributes(spatialCovs[[j]])$z)]] <- model$data[[names(attributes(spatialCovs[[j]])$z)]]
          else stop("spatialCovs$", spatialcovnames[j], 
                    " z value '", names(attributes(spatialCovs[[j]])$z), 
                    "' not found in covs")
        }
        zname <- names(attributes(spatialCovs[[j]])$z)
        zvalues <- raster::getZ(spatialCovs[[j]])
        if (!all(unique(covs[[zname]]) %in% zvalues)) 
          stop("data$", zname, " includes z-values with no matching raster layer in spatialCovs$", 
               spatialcovnames[j])
      }
    }
  }
  else nbSpatialCovs <- 0
  if (is.list(obsPerAnimal)) {
    if (length(obsPerAnimal) != nbAnimals) 
      stop("obsPerAnimal must be a list of length ", nbAnimals)
    for (i in 1:length(obsPerAnimal)) {
      if (length(which(obsPerAnimal[[i]] < 1)) > 0) 
        stop("obsPerAnimal elements should have positive values.")
      if (length(obsPerAnimal[[i]]) == 1) 
        obsPerAnimal[[i]] <- rep(obsPerAnimal[[i]], 
                                 2)
      else if (length(obsPerAnimal[[i]]) != 2) 
        stop("obsPerAnimal elements should be of length 1 or 2.")
    }
  } else {
    if (length(which(obsPerAnimal < 1)) > 0) 
      stop("obsPerAnimal should have positive values.")
    if (length(obsPerAnimal) == 1) 
      obsPerAnimal <- rep(obsPerAnimal, 2)
    else if (length(obsPerAnimal) != 2) 
      stop("obsPerAnimal should be of length 1 or 2.")
    tmpObs <- obsPerAnimal
    obsPerAnimal <- vector("list", nbAnimals)
    for (i in 1:nbAnimals) {
      obsPerAnimal[[i]] <- tmpObs
    }
  }
  if (is.list(initialPosition)) {
    if (length(initialPosition) != nbAnimals) 
      stop("initialPosition must be a list of length ", 
           nbAnimals)
    for (i in 1:nbAnimals) {
      if (is.null(mvnCoords) || dist[[mvnCoords]] %in% 
          c("mvnorm2", "rw_mvnorm2")) {
        if (length(initialPosition[[i]]) != 2 | !is.numeric(initialPosition[[i]]) | 
            any(!is.finite(initialPosition[[i]]))) 
          stop("each element of initialPosition must be a finite numeric vector of length 2")
      }
      else if (!is.null(mvnCoords) && dist[[mvnCoords]] %in% 
               c("mvnorm3", "rw_mvnorm3")) {
        if (length(initialPosition[[i]]) != 3 | !is.numeric(initialPosition[[i]]) | 
            any(!is.finite(initialPosition[[i]]))) 
          stop("each element of initialPosition must be a finite numeric vector of length 3")
      }
    }
  }
  else {
    if (is.null(mvnCoords) || dist[[mvnCoords]] %in% c("mvnorm2", 
                                                       "rw_mvnorm2")) {
      if (length(initialPosition) != 2 | !is.numeric(initialPosition) | 
          any(!is.finite(initialPosition))) 
        stop("initialPosition must be a finite numeric vector of length 2")
    }
    else if (!is.null(mvnCoords) && dist[[mvnCoords]] %in% 
             c("mvnorm3", "rw_mvnorm3")) {
      if (all(initialPosition == 0)) 
        initialPosition <- c(0, 0, 0)
      if (length(initialPosition) != 3 | !is.numeric(initialPosition) | 
          any(!is.finite(initialPosition))) 
        stop("initialPosition must be a finite numeric vector of length 3")
    }
    tmpPos <- initialPosition
    initialPosition <- vector("list", nbAnimals)
    for (i in 1:nbAnimals) {
      initialPosition[[i]] <- tmpPos
    }
  }
  if (is.null(model)) {
    if (!is.null(covs) & nbCovs > 0) {
      if (ncol(covs) != nbCovs) 
        warning("covs and nbCovs argument conflicting - nbCovs was set to ncol(covs)")
    }
  }
  if (!is.null(covs)) {
    if (!is.data.frame(covs)) 
      stop("'covs' should be a data.frame")
  }
  if (!is.null(covs)) {
    nbCovs <- ncol(covs)
    if (length(which(is.na(covs))) > 0) 
      warning(paste("There are", length(which(is.na(covs))), 
                    "missing covariate values.", "Each will be replaced by the closest available value."))
    for (i in 1:nbCovs) {
      if (length(which(is.na(covs[, i]))) > 0) {
        if (is.na(covs[1, i])) {
          k <- 1
          while (is.na(covs[k, i])) k <- k + 1
          for (j in k:2) covs[j - 1, i] <- covs[j, i]
        }
        for (j in 2:nrow(covs)) if (is.na(covs[j, i])) 
          covs[j, i] <- covs[j - 1, i]
      }
    }
  }
  if (!is.null(betaRef)) {
    if (length(betaRef) != nbStates) 
      stop("betaRef must be a vector of length ", nbStates)
    if (!is.numeric(betaRef)) 
      stop("betaRef must be a numeric vector")
    if (min(betaRef) < 1 | max(betaRef) > nbStates) 
      stop("betaRef elements must be between 1 and ", 
           nbStates)
  }
  else {
    betaRef <- 1:nbStates
  }
  betaRef <- as.integer(betaRef)
  allNbObs <- rep(NA, nbAnimals)
  for (zoo in 1:nbAnimals) {
    if (obsPerAnimal[[zoo]][1] != obsPerAnimal[[zoo]][2]) 
      allNbObs[zoo] <- sample(obsPerAnimal[[zoo]][1]:obsPerAnimal[[zoo]][2], 
                              size = 1)
    else allNbObs[zoo] <- obsPerAnimal[[zoo]][1]
  }
  cumNbObs <- c(0, cumsum(allNbObs))
  if (!is.null(covs)) {
    covnames <- colnames(covs)
    while (sum(allNbObs) > nrow(covs)) covs <- rbind(covs, 
                                                     covs)
    covs <- data.frame(covs[1:sum(allNbObs), ])
    colnames(covs) <- covnames
    rownames(covs) <- 1:sum(allNbObs)
  }
  allCovs <- NULL
  if (nbCovs > 0) {
    if (is.null(covs)) {
      allCovs <- data.frame(cov1 = rnorm(sum(allNbObs)))
      if (nbCovs > 1) {
        for (j in 2:nbCovs) {
          c <- data.frame(rnorm(sum(allNbObs)))
          colnames(c) <- paste("cov", j, sep = "")
          allCovs <- cbind(allCovs, c)
        }
      }
    }
    else {
      allCovs <- covs
    }
  }
  if (anyDuplicated(colnames(allCovs))) 
    stop("covariates must have unique names")
  if (anyDuplicated(spatialcovnames)) 
    stop("spatialCovs must have unique names")
  if (!is.null(model) & nbSpatialCovs > 0) {
    spInd <- which(!(colnames(allCovs) %in% spatialcovnames))
    if (length(spInd)) {
      allCovs <- allCovs[, spInd, drop = FALSE]
      nbCovs <- ncol(allCovs)
    }
    else {
      allCovs <- NULL
      nbCovs <- 0
    }
  }
  else if (any(colnames(allCovs) %in% spatialcovnames)) 
    stop("spatialCovs name(s) cannot match other covariate name(s)")
  if (!all(angleCovs %in% c(colnames(allCovs), spatialcovnames))) {
    stop("angleCovs ", paste0("'", angleCovs[!(angleCovs %in% 
                                                 c(colnames(allCovs), spatialcovnames))], "'", collapse = ", "), 
         " not found in covs or spatialCovs")
  }
  centerInd <- NULL
  if (!is.null(centers)) {
    if (!is.matrix(centers)) 
      stop("centers must be a matrix")
    if (dim(centers)[2] != 2) 
      stop("centers must be a matrix consisting of 2 columns (i.e., x- and y-coordinates)")
    centerInd <- which(!apply(centers, 1, function(x) any(is.na(x))))
    if (length(centerInd)) {
      if (is.null(rownames(centers))) 
        centerNames <- paste0("center", rep(centerInd, 
                                            each = 2), ".", rep(c("dist", "angle"), length(centerInd)))
      else centerNames <- paste0(rep(rownames(centers), 
                                     each = 2), ".", rep(c("dist", "angle"), length(centerInd)))
      centerCovs <- data.frame(matrix(NA, nrow = sum(allNbObs), 
                                      ncol = length(centerInd) * 2, dimnames = list(NULL, 
                                                                                    centerNames)))
    }
  }
  else centerNames <- NULL
  centroidInd <- NULL
  if (!is.null(centroids)) {
    if (!is.list(centroids)) 
      stop("centroids must be a named list")
    centroidNames <- character()
    for (j in 1:length(centroids)) {
      if (!is.data.frame(centroids[[j]])) 
        stop("each element of centroids must be a data frame")
      if (dim(centroids[[j]])[1] < max(unlist(obsPerAnimal)) | 
          dim(centroids[[j]])[2] != 2) 
        stop("each element of centroids must be a data frame consisting of at least ", 
             max(unlist(obsPerAnimal)), " rows (i.e., the maximum number of observations per animal) and 2 columns (i.e., x- and y-coordinates)")
      if (!all(c("x", "y") %in% colnames(centroids[[j]]))) 
        stop("centroids columns must be named 'x' (x-coordinate) and 'y' (y-coordinate)")
      if (any(is.na(centroids[[j]]))) 
        stop("centroids cannot contain missing values")
      if (is.null(names(centroids[j]))) 
        centroidNames <- c(centroidNames, paste0("centroid", 
                                                 rep(j, each = 2), ".", c("dist", "angle")))
      else centroidNames <- c(centroidNames, paste0(rep(names(centroids[j]), 
                                                        each = 2), ".", c("dist", "angle")))
    }
    centroidCovs <- data.frame(matrix(NA, nrow = sum(allNbObs), 
                                      ncol = length(centroidNames), dimnames = list(NULL, 
                                                                                    centroidNames)))
    centroidInd <- length(centroidNames)/2
  }
  else centroidNames <- NULL
  if (!is.null(model) & length(centerInd)) {
    cInd <- which(!(colnames(allCovs) %in% centerNames))
    if (length(cInd)) {
      allCovs <- allCovs[, cInd, drop = FALSE]
      nbCovs <- ncol(allCovs)
    }
    else {
      allCovs <- NULL
      nbCovs <- 0
    }
  }
  else if (any(colnames(allCovs) %in% centerNames)) 
    stop("centers name(s) cannot match other covariate name(s)")
  if (!is.null(model) & length(centroidInd)) {
    cInd <- which(!(colnames(allCovs) %in% centroidNames))
    if (length(cInd)) {
      allCovs <- allCovs[, cInd, drop = FALSE]
      nbCovs <- ncol(allCovs)
    }
    else {
      allCovs <- NULL
      nbCovs <- 0
    }
  }
  else if (any(colnames(allCovs) %in% centroidNames)) 
    stop("centroids name(s) cannot match other covariate name(s)")
  allNbCovs <- nbCovs + nbSpatialCovs
  if (is.null(delta)) 
    delta0 <- matrix(rep(1, nbStates)/nbStates, mixtures, 
                     nbStates)
  else delta0 <- delta
  zeroMass <- oneMass <- vector("list", length(inputs$dist))
  names(zeroMass) <- names(oneMass) <- distnames
  allStates <- NULL
  allSpatialcovs <- NULL
  if (all(c("step", "angle") %in% distnames)) {
    distnames <- c("step", "angle", distnames[!(distnames %in% 
                                                  c("step", "angle"))])
  }
  data <- data.frame(ID = factor())
  for (i in distnames) {
    if (dist[[i]] %in% mvndists) {
      data[[paste0(i, ".x")]] <- numeric()
      data[[paste0(i, ".y")]] <- numeric()
      if (dist[[i]] %in% c("mvnorm3", "rw_mvnorm3")) 
        data[[paste0(i, ".z")]] <- numeric()
    }
    else data[[i]] <- numeric()
  }
  if ("angle" %in% distnames) {
    if (inputs$dist[["angle"]] %in% angledists & ("step" %in% 
                                                  distnames)) 
      if (inputs$dist[["step"]] %in% stepdists) {
        data$x <- numeric()
        data$y <- numeric()
      }
  }
  else if ("step" %in% distnames) {
    if (inputs$dist[["step"]] %in% stepdists) {
      data$x <- numeric()
      data$y <- numeric()
    }
  }
  else if (!is.null(mvnCoords)) {
    data[[paste0(mvnCoords, ".x")]] <- numeric()
    data[[paste0(mvnCoords, ".y")]] <- numeric()
    if (dist[[mvnCoords]] %in% c("mvnorm3", "rw_mvnorm3")) 
      data[[paste0(mvnCoords, ".z")]] <- numeric()
  }
  else if (nbSpatialCovs | length(centerInd) | length(centroidInd) | 
           length(angleCovs)) 
    stop("spatialCovs, angleCovs, centers, and/or centroids cannot be specified without valid step length and turning angle distributions")
  rwInd <- any(unlist(dist) %in% rwdists)
  printMessage(nbStates, dist, p, DM, formula, formDelta, 
               formPi, mixtures, "Simulating", FALSE)
  if (length(all.vars(formula))) 
    if (!all(all.vars(formula) %in% c("ID", names(allCovs), 
                                      centerNames, centroidNames, spatialcovnames))) 
      stop("'formula' covariate(s) not found")
  if (length(all.vars(formPi))) 
    if (!all(all.vars(formPi) %in% c("ID", names(allCovs), 
                                     centerNames, centroidNames, spatialcovnames))) 
      stop("'formulaPi' covariate(s) not found")
  if (length(all.vars(formDelta))) 
    if (!all(all.vars(formDelta) %in% c("ID", names(allCovs), 
                                        centerNames, centroidNames, spatialcovnames))) 
      stop("'formulaDelta' covariate(s) not found")
  if (("ID" %in% all.vars(formula) | "ID" %in% all.vars(formPi) | 
       "ID" %in% all.vars(formDelta)) & nbAnimals < 2) 
    stop("ID cannot be a covariate when nbAnimals=1")
  newForm <- newFormulas(formula, nbStates, betaRef)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  tmpCovs <- data.frame(ID = factor(1, levels = 1:nbAnimals))
  if (!is.null(allCovs)) 
    tmpCovs <- cbind(tmpCovs, allCovs[1, , drop = FALSE])
  if (nbSpatialCovs) {
    for (j in 1:nbSpatialCovs) {
      for (i in 1:length(initialPosition)) {
        if (is.na(raster::cellFromXY(spatialCovs[[j]], 
                                     initialPosition[[i]]))) 
          stop("initialPosition for individual ", i, 
               " is not within the spatial extent of the ", 
               spatialcovnames[j], " raster")
      }
      getCell <- raster::cellFromXY(spatialCovs[[j]], 
                                    initialPosition[[1]])
      spCov <- spatialCovs[[j]][getCell]
      if (inherits(spatialCovs[[j]], c("RasterStack", 
                                       "RasterBrick"))) {
        zname <- names(attributes(spatialCovs[[j]])$z)
        zvalues <- raster::getZ(spatialCovs[[j]])
        spCov <- spCov[1, which(zvalues == tmpCovs[[zname]][1])]
      }
      tmpCovs[[spatialcovnames[j]]] <- spCov
    }
  }
  if (length(centerInd)) {
    for (j in 1:length(centerInd)) {
      tmpDistAngle <- distAngle(initialPosition[[1]], 
                                initialPosition[[1]], centers[centerInd[j], 
                                                              ])
      tmpCovs[[centerNames[(j - 1) * 2 + 1]]] <- tmpDistAngle[1]
      tmpCovs[[centerNames[(j - 1) * 2 + 2]]] <- tmpDistAngle[2]
    }
  }
  if (length(centroidInd)) {
    for (j in 1:centroidInd) {
      tmpDistAngle <- distAngle(initialPosition[[1]], 
                                initialPosition[[1]], as.numeric(centroids[[j]][1, 
                                                                                ]))
      tmpCovs[[centroidNames[(j - 1) * 2 + 1]]] <- tmpDistAngle[1]
      tmpCovs[[centroidNames[(j - 1) * 2 + 2]]] <- tmpDistAngle[2]
    }
  }
  if (mixtures > 1) {
    if (!is.null(beta)) {
      if (!is.list(beta)) 
        stop("beta must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
    }
  }
  if (!is.null(recharge)) {
    g0covs <- model.matrix(recharge$g0, tmpCovs[1, , drop = FALSE])
    nbG0covs <- ncol(g0covs) - 1
    recovs <- model.matrix(recharge$theta, tmpCovs)
    nbRecovs <- ncol(recovs) - 1
    if (!nbRecovs) 
      stop("invalid recharge model -- theta must include an intercept and at least 1 covariate")
    tmpcovs <- cbind(tmpCovs, rep(0, nrow(tmpCovs)))
    colnames(tmpcovs) <- c(colnames(tmpCovs), "recharge")
    tmpCovs <- tmpcovs
    if (!is.null(beta)) {
      if (!is.list(beta)) 
        stop("beta must be a list with elements named 'beta', 'g0', and/or 'theta' when a recharge model is specified")
    }
    for (j in 1:nbStates) {
      formulaStates[[j]] <- as.formula(paste0(Reduce(paste, 
                                                     deparse(formulaStates[[j]])), "+recharge"))
    }
    formterms <- c(formterms, "recharge")
    newformula <- as.formula(paste0(Reduce(paste, deparse(newformula)), 
                                    "+recharge"))
    beta0 <- beta
  }
  else {
    nbRecovs <- 0
    nbG0covs <- 0
    g0covs <- NULL
    recovs <- NULL
    if (is.null(model) & !is.list(beta)) {
      beta0 <- list(beta = beta)
    }
    else beta0 <- beta
  }
  nbBetaCovs <- ncol(model.matrix(newformula, tmpCovs))
  if (is.null(beta0$beta)) 
    beta0$beta <- matrix(rnorm(nbStates * (nbStates - 1) * 
                                 nbBetaCovs * mixtures), nrow = nbBetaCovs * mixtures)
  if (nbRecovs) {
    if (is.null(beta0$g0) | is.null(beta0$theta)) 
      stop("beta$g0 and beta$theta must be specified for recharge model")
    if (length(beta0$g0) != (nbG0covs + 1) | any(!is.numeric(beta0$g0))) 
      stop("beta$g0 must be a numeric vector of length ", 
           nbG0covs + 1)
    if (length(beta0$theta) != (nbRecovs + 1) | any(!is.numeric(beta0$theta))) 
      stop("beta$theta must be a numeric vector of length ", 
           nbRecovs + 1)
  }
  if (ncol(beta0$beta) != nbStates * (nbStates - 1) | (nrow(beta0$beta)/mixtures) != 
      nbBetaCovs) {
    error <- paste("beta has wrong dimensions: it should have", 
                   nbBetaCovs * mixtures, "rows and", nbStates * (nbStates - 
                                                                    1), "columns.")
    stop(error)
  }
  if (nbStates > 1) {
    for (state in 1:(nbStates * (nbStates - 1))) {
      noBeta <- which(match(colnames(model.matrix(newformula, 
                                                  tmpCovs)), colnames(model.matrix(formulaStates[[state]], 
                                                                                   tmpCovs)), nomatch = 0) == 0)
      if (length(noBeta)) 
        beta0$beta[noBeta, state] <- 0
    }
  }
  covsPi <- model.matrix(formPi, tmpCovs)
  nbCovsPi <- ncol(covsPi) - 1
  if (!nbCovsPi & is.null(formulaPi)) {
    if (is.null(beta0$pi)) {
      beta0$pi <- matrix(1/mixtures, (nbCovsPi + 1), mixtures)
    }
    else {
      beta0$pi <- matrix(beta0$pi, (nbCovsPi + 1), mixtures)
    }
    if (length(beta0$pi) != (nbCovsPi + 1) * mixtures) 
      stop(paste("beta$pi has the wrong length: it should have", 
                 mixtures, "elements."))
    beta0$pi <- matrix(log(beta0$pi[-1]/beta0$pi[1]), nbCovsPi + 
                         1, mixtures - 1)
  }
  else {
    if (is.null(beta0$pi)) 
      beta0$pi <- matrix(0, nrow = (nbCovsPi + 1), ncol = mixtures - 
                           1)
    if (is.null(dim(beta0$pi)) || (ncol(beta0$pi) != mixtures - 
                                   1 | nrow(beta0$pi) != (nbCovsPi + 1))) 
      stop(paste("beta$pi has wrong dimensions: it should have", 
                 (nbCovsPi + 1), "rows and", mixtures - 1, "columns."))
  }
  covsDelta <- model.matrix(formDelta, tmpCovs)
  nbCovsDelta <- ncol(covsDelta) - 1
  if (!nbCovsDelta & is.null(formulaDelta)) {
    if (mixtures == 1) {
      if (length(delta0) != (nbCovsDelta + 1) * nbStates) 
        stop(paste("delta has the wrong length: it should have", 
                   nbStates * mixtures, "elements."))
      deltaB <- matrix(log(delta0[-1]/delta0[1]), 1)
    }
    else {
      if (is.null(dim(delta0)) || (ncol(delta0) != nbStates | 
                                   nrow(delta0) != mixtures)) 
        stop(paste("delta has wrong dimensions: it should have", 
                   mixtures, "rows and", nbStates, "columns."))
      deltaB <- apply(delta0, 1, function(x) log(x[-1]/x[1]))
    }
  }
  else {
    if (is.null(dim(delta0)) || (ncol(delta0) != nbStates - 
                                 1 | nrow(delta0) != (nbCovsDelta + 1) * mixtures)) 
      stop(paste("delta has wrong dimensions: it should have", 
                 (nbCovsDelta + 1) * mixtures, "rows and", nbStates - 
                   1, "columns."))
    deltaB <- delta0
  }
  parCount <- lapply(Par[distnames], length)
  parindex <- c(0, cumsum(unlist(parCount))[-length(distnames)])
  names(parindex) <- distnames
  if (is.null(workBounds)) {
    workBounds <- vector("list", length(distnames))
    names(workBounds) <- distnames
  }
  workBounds <- getWorkBounds(workBounds, distnames, unlist(Par[distnames]), 
                              parindex, parCount, inputs$DM, beta0, deltaB)
  wnbeta <- w2wn(beta0$beta, workBounds$beta)
  wnpi <- w2wn(beta0$pi, workBounds$pi)
  if (!is.null(recharge)) {
    wng0 <- w2wn(beta0$g0, workBounds$g0)
    wntheta <- w2wn(beta0$theta, workBounds$theta)
  }
  mix <- rep(1, nbAnimals)
  if (!nbSpatialCovs | !retrySims) {
    ######### Run through animal simulations #########
    # Add progress bar
    pb <- txtProgressBar(max=nbAnimals, style = 3)
    setTxtProgressBar(pb, 0)
    ##################################################
    for (zoo in 1:nbAnimals) {
      nbObs <- allNbObs[zoo]
      d <- data.frame(ID = factor(rep(zoo, nbObs)))
      subCovs <- data.frame(ID = rep(factor(zoo, levels = 1:nbAnimals), 
                                     nbObs))
      if (nbCovs > 0) {
        if (zoo < 2) 
          ind1 <- 1
        else ind1 <- sum(allNbObs[1:(zoo - 1)]) + 1
        ind2 <- sum(allNbObs[1:zoo])
        subCovs <- cbind(subCovs, data.frame(allCovs[ind1:ind2, 
                                                     , drop = FALSE]))
      }
      if (length(centerInd)) 
        subCovs <- cbind(subCovs, centerCovs[cumNbObs[zoo] + 
                                               1:nbObs, ])
      if (length(centroidInd)) 
        subCovs <- cbind(subCovs, centroidCovs[cumNbObs[zoo] + 
                                                 1:nbObs, ])
      subSpatialcovs <- as.data.frame(matrix(NA, nrow = nbObs, 
                                             ncol = nbSpatialCovs))
      colnames(subSpatialcovs) <- spatialcovnames
      subAnglecovs <- as.data.frame(matrix(NA, nrow = nbObs, 
                                           ncol = length(angleCovs)))
      colnames(subAnglecovs) <- angleCovs
      X <- matrix(0, nrow = nbObs, ncol = 2)
      if (!is.null(mvnCoords) && dist[[mvnCoords]] %in% 
          c("mvnorm3", "rw_mvnorm3")) 
        X <- matrix(0, nrow = nbObs, ncol = 3)
      X[1, ] <- initialPosition[[zoo]]
      phi <- 0
      genData <- genArgs <- vector("list", length(distnames))
      names(genData) <- names(genArgs) <- distnames
      for (i in distnames) {
        if (inputs$dist[[i]] %in% mvndists) {
          genData[[i]] <- matrix(NA, nbObs, ncol(X))
        }
        else genData[[i]] <- rep(NA, nbObs)
        genArgs[[i]] <- list(1)
      }
      gamma <- matrix(0, nbStates, nbStates)
      gamma[cbind(1:nbStates, betaRef)] <- 1
      gamma <- t(gamma)
      if (!nbSpatialCovs & !length(centerInd) & !length(centroidInd) & 
          !length(angleCovs) & !rwInd) {
        if (!is.null(recharge)) {
          g0 <- model.matrix(recharge$g0, subCovs[1, 
                                                  , drop = FALSE]) %*% wng0
          subCovs[, "recharge"] <- cumsum(c(g0, model.matrix(recharge$theta, 
                                                             subCovs[-nrow(subCovs), ]) %*% wntheta))
        }
        DMcov <- model.matrix(newformula, subCovs)
        DMinputs <- getDM(subCovs, inputs$DM, inputs$dist, 
                          nbStates, p$parNames, p$bounds, Par, cons, 
                          workcons, zeroInflation, oneInflation, inputs$circularAngleMean)
        fullDM <- DMinputs$fullDM
        DMind <- DMinputs$DMind
        wpar <- n2w(Par, bounds, beta0, deltaB, nbStates, 
                    inputs$estAngleMean, inputs$DM, DMinputs$cons, 
                    DMinputs$workcons, p$Bndind, inputs$dist)
        if (any(!is.finite(wpar))) 
          stop("Scaling error. Check initial parameter values and bounds.")
        ncmean <- get_ncmean(distnames, fullDM, inputs$circularAngleMean, 
                             nbStates)
        nc <- ncmean$nc
        meanind <- ncmean$meanind
        covsDelta <- model.matrix(formDelta, subCovs[1, 
                                                     , drop = FALSE])
        covsPi <- model.matrix(formPi, subCovs[1, , 
                                               drop = FALSE])
        fullsubPar <- w2n(wpar, bounds, parSize, nbStates, 
                          nbBetaCovs - 1, inputs$estAngleMean, inputs$circularAngleMean, 
                          inputs$consensus, stationary = FALSE, DMinputs$cons, 
                          fullDM, DMind, DMinputs$workcons, nbObs, inputs$dist, 
                          p$Bndind, nc, meanind, covsDelta, workBounds, 
                          covsPi)
        pie <- fullsubPar$pi
        if (mixtures > 1) 
          mix[zoo] <- sample.int(mixtures, 1, prob = pie)
        gFull <- DMcov %*% wnbeta[(mix[zoo] - 1) * nbBetaCovs + 
                                    1:nbBetaCovs, ]
        g <- gFull[1, , drop = FALSE]
        delta0 <- fullsubPar$delta[mix[zoo], ]
      }
      else {
        if (nbSpatialCovs) {
          for (j in 1:nbSpatialCovs) {
            getCell <- raster::cellFromXY(spatialCovs[[j]], 
                                          c(X[1, 1], X[1, 2]))
            if (is.na(getCell)) 
              stop("Movement is beyond the spatial extent of the ", 
                   spatialcovnames[j], " raster. Try expanding the extent of the raster.")
            spCov <- spatialCovs[[j]][getCell]
            if (inherits(spatialCovs[[j]], c("RasterStack", 
                                             "RasterBrick"))) {
              zname <- names(attributes(spatialCovs[[j]])$z)
              zvalues <- raster::getZ(spatialCovs[[j]])
              spCov <- spCov[1, which(zvalues == subCovs[1, 
                                                         zname])]
            }
            subSpatialcovs[1, j] <- spCov
            if (spatialcovnames[j] %in% angleCovs) {
              subAnglecovs[1, spatialcovnames[j]] <- subSpatialcovs[1, 
                                                                    j]
              subSpatialcovs[1, j] <- 0
            }
          }
        }
        for (j in angleCovs[which(angleCovs %in% names(subCovs))]) {
          subAnglecovs[1, j] <- subCovs[1, j]
          subCovs[1, j] <- 0
        }
        if (length(centerInd)) {
          for (j in 1:length(centerInd)) {
            subCovs[1, centerNames[(j - 1) * 2 + 1:2]] <- distAngle(X[1, 
                                                                      ], X[1, ], centers[centerInd[j], ])
          }
        }
        if (length(centroidInd)) {
          for (j in 1:centroidInd) {
            subCovs[1, centroidNames[(j - 1) * 2 + 1:2]] <- distAngle(X[1, 
                                                                        ], X[1, ], as.numeric(centroids[[j]][1, 
                                                                                                             ]))
          }
        }
        if (!is.null(recharge)) {
          g0 <- model.matrix(recharge$g0, cbind(subCovs[1, 
                                                        , drop = FALSE], subSpatialcovs[1, , drop = FALSE])) %*% 
            wng0
          subCovs[1, "recharge"] <- g0
        }
        for (i in distnames) {
          if (dist[[i]] %in% rwdists) {
            if (dist[[i]] %in% c("rw_mvnorm2")) {
              subCovs[1, paste0(i, ".x_tm1")] <- X[1, 
                                                   1]
              subCovs[1, paste0(i, ".y_tm1")] <- X[1, 
                                                   2]
            }
            else if (dist[[i]] %in% c("rw_mvnorm3")) {
              subCovs[1, paste0(i, ".x_tm1")] <- X[1, 
                                                   1]
              subCovs[1, paste0(i, ".y_tm1")] <- X[1, 
                                                   2]
              subCovs[1, paste0(i, ".z_tm1")] <- X[1, 
                                                   3]
            }
          }
        }
        maxlag <- getDM(cbind(subCovs[1, , drop = FALSE], 
                              subSpatialcovs[1, , drop = FALSE]), inputs$DM, 
                        inputs$dist, nbStates, p$parNames, p$bounds, 
                        Par, cons, workcons, zeroInflation, oneInflation, 
                        inputs$circularAngleMean, wlag = TRUE)$lag
        covsPi <- model.matrix(formPi, cbind(subCovs[1, 
                                                     , drop = FALSE], subSpatialcovs[1, , drop = FALSE]))
        pie <- mlogit(wnpi, covsPi, nbCovsPi, 1, mixtures)
        if (mixtures > 1) 
          mix[zoo] <- sample.int(mixtures, 1, prob = pie)
        g <- model.matrix(newformula, cbind(subCovs[1, 
                                                    , drop = FALSE], subSpatialcovs[1, , drop = FALSE])) %*% 
          wnbeta[(mix[zoo] - 1) * nbBetaCovs + 1:nbBetaCovs, 
                 ]
        covsDelta <- model.matrix(formDelta, cbind(subCovs[1, 
                                                           , drop = FALSE], subSpatialcovs[1, , drop = FALSE]))
        delta0 <- mlogit(deltaB[(mix[zoo] - 1) * (nbCovsDelta + 
                                                    1) + 1:(nbCovsDelta + 1), , drop = FALSE], 
                         covsDelta, nbCovsDelta, 1, nbStates)
      }
      gamma[!gamma] <- exp(g)
      gamma <- t(gamma)
      gamma <- gamma/apply(gamma, 1, sum)
      if (nbStates > 1) {
        Z <- rep(NA, nbObs)
        Z[1] <- sample(1:nbStates, size = 1, prob = delta0 %*% 
                         gamma)
      }
      else Z <- rep(1, nbObs)
      ######
      # Data generation loop simulations
      ######
      for (k in 1:(nbObs - 1)) {
        ##############################
        # Nesting loop in While loop #
        ##############################
        onLand <- TRUE
        while(onLand){
        ##############################
          if (nbSpatialCovs | length(centerInd) | length(centroidInd) | 
              length(angleCovs) | rwInd) {
            DMinputs <- getDM(cbind(subCovs[k - maxlag:0, 
                                            , drop = FALSE], subSpatialcovs[k - maxlag:0, 
                                                                            , drop = FALSE]), inputs$DM, inputs$dist, 
                              nbStates, p$parNames, p$bounds, Par, cons, 
                              workcons, zeroInflation, oneInflation, inputs$circularAngleMean, 
                              wlag = TRUE)
            fullDM <- DMinputs$fullDM
            DMind <- DMinputs$DMind
            wpar <- n2w(Par, bounds, beta0, deltaB, nbStates, 
                        inputs$estAngleMean, inputs$DM, DMinputs$cons, 
                        DMinputs$workcons, p$Bndind, inputs$dist)
            if (any(!is.finite(wpar))) 
              stop("Scaling error. Check initial parameter values and bounds.")
            nc <- meanind <- vector("list", length(distnames))
            names(nc) <- names(meanind) <- distnames
            for (i in distnames) {
              nc[[i]] <- apply(fullDM[[i]], 1:2, function(x) !all(unlist(x) == 
                                                                    0))
              if (!isFALSE(inputs$circularAngleMean[[i]])) {
                meanind[[i]] <- which((apply(fullDM[[i]][1:nbStates, 
                                                         , drop = FALSE], 1, function(x) !all(unlist(x) == 
                                                                                                0))))
                if (length(meanind[[i]])) {
                  angInd <- which(is.na(match(gsub("cos", 
                                                   "", gsub("sin", "", colnames(nc[[i]]))), 
                                              colnames(nc[[i]]), nomatch = NA)))
                  sinInd <- colnames(nc[[i]])[which(grepl("sin", 
                                                          colnames(nc[[i]])[angInd]))]
                  nc[[i]][meanind[[i]], sinInd] <- ifelse(nc[[i]][meanind[[i]], 
                                                                  sinInd], nc[[i]][meanind[[i]], sinInd], 
                                                          nc[[i]][meanind[[i]], gsub("sin", 
                                                                                     "cos", sinInd)])
                  nc[[i]][meanind[[i]], gsub("sin", "cos", 
                                             sinInd)] <- ifelse(nc[[i]][meanind[[i]], 
                                                                        gsub("sin", "cos", sinInd)], nc[[i]][meanind[[i]], 
                                                                                                             gsub("sin", "cos", sinInd)], nc[[i]][meanind[[i]], 
                                                                                                                                                  sinInd])
                }
              }
            }
            subPar <- w2n(wpar, bounds, parSize, nbStates, 
                          nbBetaCovs - 1, inputs$estAngleMean, inputs$circularAngleMean, 
                          inputs$consensus, stationary = FALSE, DMinputs$cons, 
                          fullDM, DMind, DMinputs$workcons, 1, inputs$dist, 
                          p$Bndind, nc, meanind, covsDelta, workBounds, 
                          covsPi)
          }
          else {
            subPar <- lapply(fullsubPar[distnames], function(x) x[,k, drop = FALSE])
          }
          ####### Loop over step and angles to generate values
          for (i in distnames) {
            zeroMass[[i]] <- rep(0, nbStates)
            oneMass[[i]] <- rep(0, nbStates)
            if (zeroInflation[[i]] | oneInflation[[i]]) {
              if (zeroInflation[[i]]) 
                zeroMass[[i]] <- subPar[[i]][parSize[[i]] * 
                                               nbStates - nbStates * oneInflation[[i]] - 
                                               (nbStates - 1):0]
              if (oneInflation[[i]]) 
                oneMass[[i]] <- subPar[[i]][parSize[[i]] * 
                                              nbStates - (nbStates - 1):0]
              subPar[[i]] <- subPar[[i]][-(parSize[[i]] * 
                                             nbStates - (nbStates * oneInflation[[i]] - 
                                                           nbStates * zeroInflation[[i]] - 1):0)]
            } # F
            if (inputs$dist[[i]] %in% mvndists) {
              if (inputs$dist[[i]] == "mvnorm2" || inputs$dist[[i]] == 
                  "rw_mvnorm2") {
                genArgs[[i]][[2]] <- c(subPar[[i]][Z[k]], 
                                       subPar[[i]][nbStates + Z[k]])
                genArgs[[i]][[3]] <- matrix(c(subPar[[i]][nbStates * 
                                                            2 + Z[k]], subPar[[i]][nbStates * 3 + 
                                                                                     Z[k]], subPar[[i]][nbStates * 3 + Z[k]], 
                                              subPar[[i]][nbStates * 4 + Z[k]]), 2, 
                                            2)
              }
              else if (inputs$dist[[i]] == "mvnorm3" || 
                       inputs$dist[[i]] == "rw_mvnorm3") {
                genArgs[[i]][[2]] <- c(subPar[[i]][Z[k]], 
                                       subPar[[i]][nbStates + Z[k]], subPar[[i]][2 * 
                                                                                   nbStates + Z[k]])
                genArgs[[i]][[3]] <- matrix(c(subPar[[i]][nbStates * 
                                                            3 + Z[k]], subPar[[i]][nbStates * 4 + 
                                                                                     Z[k]], subPar[[i]][nbStates * 5 + Z[k]], 
                                              subPar[[i]][nbStates * 4 + Z[k]], subPar[[i]][nbStates * 
                                                                                              6 + Z[k]], subPar[[i]][nbStates * 
                                                                                                                       7 + Z[k]], subPar[[i]][nbStates * 
                                                                                                                                                5 + Z[k]], subPar[[i]][nbStates * 
                                                                                                                                                                         7 + Z[k]], subPar[[i]][nbStates * 
                                                                                                                                                                                                  8 + Z[k]]), 3, 3)
              }
            } # F
            else if (inputs$dist[[i]] == "cat") {
              genArgs[[i]][[2]] <- subPar[[i]][seq(Z[k], 
                                                   (parSize[[i]] + 1) * nbStates, nbStates)]
            } # F
            else {
              for (j in 1:(parSize[[i]] - zeroInflation[[i]] - oneInflation[[i]])) 
                genArgs[[i]][[j + 1]] <- subPar[[i]][(j - 1) * nbStates + Z[k]]
            }
            # Angle is generated here
            if (inputs$dist[[i]] %in% angledists) {
              genData[[i]][k] <- do.call(Fun[[i]], genArgs[[i]])
              # correction if angle > or < than pi
              if (genData[[i]][k] > pi) 
                genData[[i]][k] <- genData[[i]][k] - 2 * pi
              if (genData[[i]][k] < -pi) 
                genData[[i]][k] <- genData[[i]][k] + 2 * pi
              if (i == "angle" & ("step" %in% distnames)) {
                if (inputs$dist[["step"]] %in% stepdists) {
                  if (genData$step[k] > 0) {
                    phi <- phi + genData[[i]][k]
                  }
                  # Here is where the next coordinate is generated
                  # m is the new pos coord units to add
                  m <- genData$step[k] * c(Re(exp((0+1i) * phi)), Im(exp((0+1i) * phi)))
                  X[k + 1, ] <- X[k, ] + m
                }
              }
            }
            else {
              ###### Step is generated from gamma distributions here
              if (inputs$dist[[i]] == "gamma") {
                shape <- genArgs[[i]][[2]]^2/genArgs[[i]][[3]]^2
                scale <- genArgs[[i]][[3]]^2/genArgs[[i]][[2]]
                genArgs[[i]][[2]] <- shape
                genArgs[[i]][[3]] <- 1/scale
              }
              probs <- c(1 - zeroMass[[i]][Z[k]] - oneMass[[i]][Z[k]], 
                         zeroMass[[i]][Z[k]], oneMass[[i]][Z[k]])
              rU <- which(rmultinom(1, 1, prob = probs) == 1)
              if (rU == 1) {
                if (inputs$dist[[i]] %in% mvndists) {
                  genData[[i]][k, ] <- do.call(Fun[[i]], 
                                               genArgs[[i]])
                }
                else genData[[i]][k] <- do.call(Fun[[i]], genArgs[[i]])
              }
              else if (rU == 2) {
                genData[[i]][k] <- 0
              }
              else {
                genData[[i]][k] <- 1
              }
            }
            if (!is.null(mvnCoords) && i == mvnCoords) {
              X[k + 1, ] <- genData[[i]][k, ]
              d[[i]] <- X
            }
            #######
            # On this line step and or angle in (d) is replaced by generated data
            #######
            else d[[i]] <- genData[[i]]
          }
          #### ^^^ Here ends the loop to generate step and angles
          gamma <- matrix(0, nbStates, nbStates)
          gamma[cbind(1:nbStates, betaRef)] <- 1
          gamma <- t(gamma)
          if (nbSpatialCovs | length(centerInd) | length(centroidInd) | 
              length(angleCovs) | rwInd) {
            if (nbSpatialCovs) {
              for (j in 1:nbSpatialCovs) {
                getCell <- raster::cellFromXY(spatialCovs[[j]], 
                                              c(X[k + 1, 1], X[k + 1, 2]))
                if (is.na(getCell)) 
                  stop("Movement is beyond the spatial extent of the ", 
                       spatialcovnames[j], " raster. Try expanding the extent of the raster.")
                spCov <- spatialCovs[[j]][getCell]
                if (inherits(spatialCovs[[j]], c("RasterStack", 
                                                 "RasterBrick"))) {
                  zname <- names(attributes(spatialCovs[[j]])$z)
                  zvalues <- raster::getZ(spatialCovs[[j]])
                  spCov <- spCov[1, which(zvalues == subCovs[k + 
                                                               1, zname])]
                }
                subSpatialcovs[k + 1, j] <- spCov
                if (spatialcovnames[j] %in% angleCovs) {
                  subAnglecovs[k + 1, spatialcovnames[j]] <- subSpatialcovs[k + 1, j]
                  subSpatialcovs[k + 1, j] <- 
                    circAngles(subAnglecovs[k:(k + 1), spatialcovnames[j]], 
                               data.frame(x = X[k:(k + 1), 1], y = X[k:(k + 1), 2]))[2]
                }
              }
            }
            for (j in angleCovs[which(angleCovs %in% names(subCovs))]) {
              subAnglecovs[k + 1, j] <- subCovs[k + 1,j]
              subCovs[k + 1, j] <- 
                circAngles(subAnglecovs[k:(k + 1), j], 
                           data.frame(x = X[k:(k + 1), 1], y = X[k:(k + 1), 2]))[2]
            }
            if (length(centerInd)) {
              for (j in 1:length(centerInd)) {
                subCovs[k + 1, centerNames[(j - 1) * 2 + 1:2]] <- 
                  distAngle(X[k, ], X[k + 1,], centers[centerInd[j], ])
              }
            }
            if (length(centroidInd)) {
              for (j in 1:centroidInd) {
                subCovs[k + 1, centroidNames[(j - 1) * 2 + 1:2]] <- 
                  distAngle(X[k, ], X[k + 1, ], as.numeric(centroids[[j]][k + 1, ]))
              }
            }
            if (!is.null(recharge)) {
              subCovs[k + 1, "recharge"] <- 
                subCovs[k,"recharge"] + model.matrix(recharge$theta, 
                                                     cbind(subCovs[k, , drop = FALSE], 
                                                           subSpatialcovs[k,, drop = FALSE])) %*% wntheta
            }
            for (i in distnames) {
              if (dist[[i]] %in% rwdists) {
                if (dist[[i]] %in% c("rw_mvnorm2")) 
                  subCovs[k + 1, paste0(i, c(".x_tm1", 
                                             ".y_tm1"))] <- X[k + 1, ]
                else if (dist[[i]] %in% c("rw_mvnorm3")) 
                  subCovs[k + 1, paste0(i, c(".x_tm1", 
                                             ".y_tm1", ".z_tm1"))] <- X[k + 1,]
              }
            }
            g <- model.matrix(newformula, cbind(subCovs[k + 1, , drop = FALSE], 
                                                subSpatialcovs[k + 1, , drop = FALSE])) %*% wnbeta[(mix[zoo] - 1) * nbBetaCovs + 1:nbBetaCovs, ]
          }
          else {
            g <- gFull[k + 1, , drop = FALSE]
          }
          gamma[!gamma] <- exp(g)
          gamma <- t(gamma)
          gamma <- gamma/apply(gamma, 1, sum)
          # Z I believe is the state
          Z[k + 1] <- sample(1:nbStates, size = 1, prob = gamma[Z[k],])
          #####################################################################################
          # CODE ADDITIONS
          # USE `X` to check position on land
          k_location <- data.frame(lon=X[k+1,1]*1000, lat=X[k+1,2]*1000)
          coordinates(k_location) <- ~ lon + lat
          proj4string(k_location) <- landShpProj
          # if onLand FALSE, while loop ends and next interation begins
          onLand <- !is.na(over(k_location, landShp)$area)
        }
        #####################################################################################
      }
      allStates <- c(allStates, Z)
      if (nbSpatialCovs > 0) {
        allSpatialcovs <- rbind(allSpatialcovs, subSpatialcovs)
      }
      if ("angle" %in% distnames) {
        if (inputs$dist[["angle"]] %in% angledists & 
            ("step" %in% distnames)) 
          if (inputs$dist[["step"]] %in% stepdists) {
            d$angle[1] <- NA
            step0 <- which(d$step == 0)
            d$angle[c(step0, step0 + 1)] <- NA
            d$x = X[, 1]
            d$y = X[, 2]
          }
      }
      else if ("step" %in% distnames) {
        if (inputs$dist[["step"]] %in% stepdists) {
          d$x = c(initialPosition[[zoo]][1], initialPosition[[zoo]][1] + 
                    cumsum(d$step)[-nrow(d)])
          d$y = rep(initialPosition[[zoo]][2], nrow(d))
        }
      }
      else if (!is.null(mvnCoords)) {
        d[[paste0(mvnCoords, ".x")]] <- d[[mvnCoords]][, 
                                                       1]
        d[[paste0(mvnCoords, ".y")]] <- d[[mvnCoords]][, 
                                                       2]
        if (dist[[mvnCoords]] %in% c("mvnorm3", "rw_mvnorm3")) 
          d[[paste0(mvnCoords, ".z")]] <- d[[mvnCoords]][, 
                                                         3]
        d[[mvnCoords]] <- NULL
      }
      for (j in angleCovs[which(angleCovs %in% names(subCovs))]) allCovs[cumNbObs[zoo] + 
                                                                           1:nbObs, j] <- subCovs[, j]
      if (length(centerInd)) 
        centerCovs[cumNbObs[zoo] + 1:nbObs, ] <- subCovs[, 
                                                         centerNames]
      if (length(centroidInd)) 
        centroidCovs[cumNbObs[zoo] + 1:nbObs, ] <- subCovs[, 
                                                           centroidNames]
      data <- rbind(data, d)
      ##########################
      setTxtProgressBar(pb, zoo)
      ##########################
    }
    ##########################
    close(pb)
    ######## END OF SIMULATION GENRATION LOOP ########
    if (nbCovs > 0) 
      data <- cbind(data, allCovs)
    if (nbSpatialCovs > 0) {
      colnames(allSpatialcovs) <- spatialcovnames
      for (j in spatialcovnames) {
        if (any(raster::is.factor(spatialCovs[[j]]))) {
          allSpatialcovs[[j]] <- factor(allSpatialcovs[[j]], 
                                        levels = unique(unlist(raster::levels(spatialCovs[[j]]))))
        }
      }
      data <- cbind(data, allSpatialcovs)
    }
    if (length(centerInd)) {
      data <- cbind(data, centerCovs)
      for (j in which(grepl(".angle", names(data)))) {
        if (names(data[j]) %in% centerNames) 
          class(data[[j]]) <- c(class(data[[j]]), "angle")
      }
    }
    if (length(centroidInd)) {
      data <- cbind(data, centroidCovs)
      for (j in which(grepl(".angle", names(data)))) {
        if (names(data[j]) %in% centroidNames) 
          class(data[[j]]) <- c(class(data[[j]]), "angle")
      }
    }
    if (states) 
      data <- cbind(data, states = allStates)
    for (i in distnames) {
      if (inputs$dist[[i]] %in% angledists) 
        class(data[[i]]) <- c(class(data[[i]]), "angle")
    }
    for (i in angleCovs) {
      class(data[[i]]) <- c(class(data[[i]]), "angle")
    }
    if (!is.null(mvnCoords)) {
      attr(data, "coords") <- paste0(mvnCoords, c(".x", 
                                                  ".y"))
    }
    out <- simObsData(momentuHMMData(data), lambda, errorEllipse)
    message("DONE")
    return(out)
  }
  else {
    simCount <- 0
    cat("Attempting to simulate tracks within spatial extent(s) of raster layers(s). Press 'esc' to force exit from 'simData'\n", 
        sep = "")
    while (simCount < retrySims) {
      cat("\r    Attempt ", simCount + 1, " of ", retrySims, 
          "...", sep = "")
      tmp <- suppressMessages(tryCatch(simData(nbAnimals, 
                                               nbStates, dist, Par, beta, delta, formula, formulaDelta, 
                                               mixtures, formulaPi, covs, nbCovs, spatialCovs, 
                                               zeroInflation, oneInflation, circularAngleMean, 
                                               centers, centroids, angleCovs, obsPerAnimal, 
                                               initialPosition, DM, cons, userBounds, workBounds, 
                                               workcons, betaRef, mvnCoords, stateNames, model, 
                                               states, retrySims = 0, lambda, errorEllipse), 
                                       error = function(e) e))
      if (inherits(tmp, "error")) {
        if (grepl("Try expanding the extent of the raster", 
                  tmp)) 
          simCount <- simCount + 1
        else stop(tmp)
      }
      else {
        simCount <- retrySims
        cat("DONE\n")
        return(tmp)
      }
    }
    cat("FAILED\n")
    stop(tmp)
  }
}
