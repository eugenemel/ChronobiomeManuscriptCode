# Supporting functions for Chronobiome analyses

# Covariates adjusted for in every model of a biomarker. Dropped only when a variable is
# the stratifier or exposure: sex in sex-stratified fits, chronotype in chronotype analyses.
COVARIATES <- c("age", "sex", "smoking", "town_depr", "fastingtime", "storagetime", "chronotype", "sleepduration")

# Survival models (biomarker as predictor of an event) also adjust for draw time and self-rated health
SURVIVAL_COVARIATES <- c(COVARIATES, "hourminute", "health_rating")

# Right-hand side of a model formula from a covariate set, e.g. "age + sex + smoking + ..."
covariateTerms <- function(covariates = COVARIATES, drop = NULL) {
  paste(setdiff(covariates, drop), collapse = " + ")
}

# Fit a GAM (Generalized Additive Model) to examine diurnal rhythms of a biomarker
# dropCovariates: covariates to leave out, e.g. "chronotype" when fitting within one chronotype group
# returnFit: return the fitted bam object itself
getGamFit <- function(biomarker, rowsubset = c(), adjustforSex = FALSE, plot = FALSE, plotdata = FALSE, verbose = FALSE, prcChange = FALSE, zscale = FALSE, predict = TRUE, dropCovariates = NULL, returnFit = FALSE) {
  hasBiomarkerData <- !is.na(X[, biomarker])

  # too few observations, return empty data frame
  if (sum(hasBiomarkerData) < 10) {
    return(data.frame())
  }

  # construct model data frame
  tmpdf <- X %>% select(bm = all_of(biomarker), "hourminute", all_of(COVARIATES))

  # Filter data based on provided rowsubset
  if (length(rowsubset) > 0) {
    tmpdf <- tmpdf[rowsubset & hasBiomarkerData, ]
  } else {
    tmpdf <- tmpdf[hasBiomarkerData, ]
  }

  if (prcChange) {
    tmpdf$bm <- tmpdf$bm / mean(tmpdf$bm, na.rm = TRUE) * 100
  }
  if (zscale) {
    tmpdf$bm <- as.double(scale(tmpdf$bm, scale = TRUE))
  }

  # Define model formula
  dropCovariates <- c(dropCovariates, if (!adjustforSex) "sex")
  formula_str <- paste("bm ~ s(hourminute, k=24, bs='ps') +", covariateTerms(drop = dropCovariates))

  tryCatch(
    expr = {
      fit <- mgcv::bam(as.formula(formula_str), data = tmpdf, gamma = 1)

      if (returnFit) {
        return(fit)
      }
      if (plot) {
        mgcv::plot.gam(fit, pages = 1)
      }
      if (verbose) {
        print(summary(fit))
      }

      if (plotdata) {
        pdf(NULL)
        plotdata <- mgcv::plot.gam(fit)[[1]]
        plotdf <- data.frame(x = plotdata$x, fit = plotdata$fit, se = plotdata$se)
        dev.off()
        return(plotdf)
      }
      
      if (predict) { 
       return(predict(fit))
      }

      tmp1 <- broom::tidy(fit)
      tmp1$n <- fit$df.null
      betas <- as.vector(fit$coefficients[grepl("hour", names(fit$coefficients))])
      names(betas) <- paste0("beta", seq_along(betas))
      tmp1 <- as.data.frame(c(as.list(tmp1), as.list(betas)))
      tmp1$range <- max(betas) - min(betas)
      tmp1$biomarker <- biomarker

      return(tmp1)
    },
    error = function(e) {
      print(e)
      return(data.frame())
    }
  ) # end tryCatch
}

# Fit the GAM for one population subset and label the coefficient table
fitPopulation <- function(biomarker, rowsubset, population, adjustforSex = FALSE) {
  tmp <- getGamFit(biomarker, rowsubset = rowsubset, adjustforSex = adjustforSex, predict = FALSE)
  if (nrow(tmp) > 0) {
    tmp$population <- population
  }
  return(tmp)
}

# Aggregate GAM fits across all demographic/health subgroups
# good = Good/Excellent, poor = Poor/Fair self-rated health; "_all" strata are adjusted for sex
gamwrap <- function(biomarker) {
  goodHealth <- X$health_rating %in% c("Good", "Excellent")
  poorHealth <- X$health_rating %in% c("Poor", "Fair")
  return(dplyr::bind_rows(
    fitPopulation(biomarker, X$sex == "Male" & goodHealth, "good_health_males"),
    fitPopulation(biomarker, X$sex == "Female" & goodHealth, "good_health_females"),
    fitPopulation(biomarker, X$sex == "Female" & poorHealth, "poor_health_females"),
    fitPopulation(biomarker, X$sex == "Male" & poorHealth, "poor_health_males"),
    fitPopulation(biomarker, goodHealth, "good_health_all", adjustforSex = TRUE),
    fitPopulation(biomarker, poorHealth, "poor_health_all", adjustforSex = TRUE)
  ))
}

# Plot GAM fits comparing good vs poor health or males vs females
plotGamFit <- function(biomarker, healthStratifed = FALSE, prcChange = FALSE, zscale = FALSE, limitSex = NA) {
  mycolors <- c("#E41A1C", "#377EB8")
  if (healthStratifed) {
    adjustforSex <- TRUE
    rowsubsetH <- X$health_rating %in% c("Good", "Excellent")
    if (!is.na(limitSex)) {
      rowsubsetH <- rowsubsetH & X$sex == limitSex
      adjustforSex <- FALSE
    }

    tmp1 <- getGamFit(biomarker,
      rowsubset = rowsubsetH, adjustforSex = adjustforSex,
      plotdata = TRUE, prcChange = prcChange, zscale = zscale
    )
    tmp1$population <- "good_health"
    
    rowsubsetF <- X$health_rating %in% c("Fair", "Poor")
    if (!is.na(limitSex)) {
      rowsubsetF <- rowsubsetF & X$sex == limitSex # Corrected from rowsubsetH
      adjustforSex <- FALSE
    }
    tmp2 <- getGamFit(biomarker,
      rowsubset = rowsubsetF, adjustforSex = adjustforSex,
      plotdata = TRUE, prcChange = prcChange, zscale = zscale
    )
    tmp2$population <- "poor_health"
    mycolors <- c("#4DAF4A", "#E41A1C")
  } else {
    rowsubsetM <- X$sex == "Male" & X$health_rating %in% c("Good", "Excellent")
    rowsubsetF <- X$sex == "Female" & X$health_rating %in% c("Good", "Excellent")
    tmp1 <- getGamFit(biomarker, rowsubset = rowsubsetF, plotdata = TRUE, prcChange = prcChange, zscale = zscale)
    tmp1$population <- "good_health_females"
    tmp2 <- getGamFit(biomarker, rowsubset = rowsubsetM, plotdata = TRUE, prcChange = prcChange, zscale = zscale)
    tmp2$population <- "good_health_males"
  }
  tmpdf <- bind_rows(tmp1, tmp2)

  ggplot(tmpdf, aes(x, fit, fill = population)) +
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.8) +
    scale_fill_manual(values = mycolors) +
    geom_hline(yintercept = 0) +
    theme_bw(base_size = 16) +
    scale_x_continuous(limits = c(9, 20), breaks = c(9, 12, 15, 18, 20)) +
    xlab("Hour") +
    ylab("Biomarker Level (z-scaled)") +
    theme(
      legend.position = "inside", legend.position.inside = c(0.8, 0.95),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    ) +
    ggtitle(biomarker)
}

# Plot simple smoothed trends unadjusted for covariates (quick exploration)

plotGamUnadjusted <- function(biomarker, rowsetDisease = NA, sexstratifed = FALSE) {
  tmpdf <- X %>% select(age, hourminute, sex, .data[[biomarker]])

  if (length(rowsetDisease) > 1) {
    mycolors <- c("#19167f", "#E41A1C", "#19167f", "#E41A1C")
    tmpdf$hadEvent <- ifelse(rowsetDisease, "sick", "healthy")
    p1 <- ggplot(tmpdf, aes(hourminute, .data[[biomarker]], fill = hadEvent)) +
      geom_smooth(span = 1.2) +
      scale_fill_manual(values = mycolors) +
      scale_color_manual(values = mycolors) +
      theme_bw(base_size = 16) +
      geom_hline(yintercept = 0) +
      scale_x_continuous(limits = c(9, 20), breaks = c(9, 12, 15, 18, 20)) +
      xlab("Hour") + ylab("Biomarker level (z-scaled)") +
      theme(legend.position = "none") +
      ggtitle(biomarker)
  } else {
    mycolors <- c("#E41A1C", "#377EB8")
    p1 <- ggplot(tmpdf, aes(hourminute, .data[[biomarker]], fill = sex)) +
      geom_smooth(span = 1.2) +
      scale_color_manual(values = mycolors) +
      scale_fill_manual(values = rev(mycolors)) +
      theme_bw(base_size = 16) +
      geom_hline(yintercept = 0) +
      scale_x_continuous(limits = c(9, 20), breaks = c(9, 12, 15, 18, 20)) +
      theme(legend.position = "none") +
      xlab("Hour") + ylab("Biomarker level (z-scaled)") +
      ggtitle(biomarker)
  }
  if (sexstratifed) {
    p1 <- p1 + facet_wrap(~sex)
  }
  print(p1)
}

# Compare GAM fits between a disease group and healthy controls

plotGamFitDisease <- function(biomarker, icd10, rowsetDisease, adjustforSex = TRUE) {
  if (adjustforSex == FALSE) {
    tmp1 <- getGamFit(biomarker, rowsubset = rowsetDisease & X$sex == "Male", plotdata = TRUE, adjustforSex = FALSE)
    tmp1$population <- paste("sick")
    tmp1$sex <- "Males"

    tmp2 <- getGamFit(biomarker, rowsubset = !rowsetDisease & X$sex == "Male", plotdata = TRUE, adjustforSex = FALSE)
    tmp2$population <- paste("healthy")
    tmp2$sex <- "Males"

    tmp3 <- getGamFit(biomarker, rowsubset = rowsetDisease & X$sex == "Female", plotdata = TRUE, adjustforSex = FALSE)
    tmp3$population <- paste("sick")
    tmp3$sex <- "Females"

    tmp4 <- getGamFit(biomarker, rowsubset = !rowsetDisease & X$sex == "Female", plotdata = TRUE, adjustforSex = FALSE)
    tmp4$population <- paste("healthy")
    tmp4$sex <- "Females"
    tmpdf <- bind_rows(tmp1, tmp2, tmp3, tmp4)
  } else {
    tmp1 <- getGamFit(biomarker, rowsubset = rowsetDisease, adjustforSex = TRUE, plotdata = TRUE)
    tmp1$population <- "sick"
    tmp2 <- getGamFit(biomarker, rowsubset = !rowsetDisease, adjustforSex = TRUE, plotdata = TRUE)
    tmp2$population <- "healthy"
    tmpdf <- bind_rows(tmp1, tmp2)
    tmpdf$sex <- "all"

    # View(tmpdf);
  }

  mycolors <- c("#19167f", "#E41A1C", "#19167f", "#E41A1C")

  p0 <- ggplot(tmpdf, aes(x, fit, fill = population)) +
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.8) +
    scale_fill_manual(values = mycolors) +
    geom_hline(yintercept = 0) +
    theme_bw(base_size = 16) +
    xlab("Hour") +
    ylab("Biomarker level (z-scaled)") +
    scale_x_continuous(limits = c(9, 20), breaks = c(9, 12, 15, 18, 20)) +
    theme(
      legend.position = "inside", legend.position.inside = c(0.8, 0.95),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    ) +
    ggtitle(paste(biomarker, icd10))
  print(p0)
}

# Generate a list of GAM plots for a set of biomarkers
plotBiomarkerSet <- function(biomarkers, healthStratifed = FALSE, ylimits = NA, diseaseRows = NA, icd10 = "", unadjusted = FALSE) {
  TINYTHEME <- theme_minimal() + theme(
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # panel.grid = element_blank(),
    # axis.text.y = element_blank(),
    # axis.text.x = element_blank()
  )

  PLOTS <- list()
  for (bm in biomarkers) {
    if (length(diseaseRows) > 1 & unadjusted == TRUE) {
      p0 <- plotGamUnadjusted(bm, rowsetDisease = diseaseRows, sexstratifed = FALSE)
    } else if (length(diseaseRows) > 1) {
      p0 <- plotGamFitDisease(bm, rowsetDisease = diseaseRows, icd10 = "")
    } else {
      p0 <- plotGamFit(bm, healthStratifed) + TINYTHEME
    }
    p0 <- p0 + TINYTHEME
    if (length(ylimits) == 2) p0 <- p0 + ylim(ylimits)
    PLOTS[[bm]] <- p0
  }
  return(PLOTS)
}

## Function to test for differences in spline shapes (permutation/bootstrap approach)
testDiffSpline <- function(biomarker, popname, rowsA, rowsB = NA, showplot = FALSE, adjustforSex = TRUE) {
  
  if (length(rowsB) <= 1) { rowsB <- !rowsA } # if rowsB not provided, take the complement of rowsA
  
  hasBiomarkerData <- !is.na(X[, biomarker]) # check we have data for this biomarker
  
  # make sure we have some data in both groups
  rowsA <- rowsA & hasBiomarkerData
  rowsB <- rowsB & hasBiomarkerData
  if (sum(rowsA) < 10 | sum(rowsB) < 10) { return(data.frame()) } # too few observations
 
  selected <- data.frame()
 # tmp1 <- getGamFit(biomarker, rowsubset = rowsA, plotdata = TRUE, adjustforSex = adjustforSex)
 # tmp2 <- getGamFit(biomarker, rowsubset = rowsB, plotdata = TRUE, adjustforSex = adjustforSex)
  idxA <- which(rowsA)
  idxB <- which(rowsB)
  
  for (itr in 1:100) {
    subA <- seq_len(nrow(X)) %in% sample(idxA, round(0.8 * length(idxA)))
    subB <- seq_len(nrow(X)) %in% sample(idxB, round(0.8 * length(idxB)))
    tmp1 <- getGamFit(biomarker, rowsubset = subA, plotdata = TRUE, adjustforSex = adjustforSex)
    tmp2 <- getGamFit(biomarker, rowsubset = subB, plotdata = TRUE, adjustforSex = adjustforSex)
                      
    # Simulate data based on fit variability
    for (i in 1:nrow(tmp1)) {
      tmp1$y[i] <- rnorm(n = 1, mean = tmp1$fit[i], sd = tmp1$se[i] * 2)
    }
    for (i in 1:nrow(tmp2)) {
      tmp2$y[i] <- rnorm(n = 1, mean = tmp2$fit[i], sd = tmp2$se[i] * 2)
    }
    tmp1$Y <- tmp1$y - tmp2$y
    lmfit <- broom::tidy(mgcv::gam(Y ~ s(x), data = tmp1))
    lmfit$sumDelta <- sum(abs(tmp1$fit - tmp2$fit))
    lmfit$noDiseaseRange <- max(tmp2$fit) - min(tmp2$fit) # popB = Controls
    lmfit$DiseaseRange <- max(tmp1$fit) - min(tmp1$fit)   # popA = Cases
    lmfit$maxDiff <- max(abs(tmp1$fit - tmp2$fit))
    selected <- bind_rows(selected, lmfit)
  }
  
  if (showplot) {
    plotdf <- tmp2
    plotdf$delta <- tmp2$fit - tmp1$fit   # difference in fit between B and A
    p0 <- ggplot(plotdf, aes(x, delta)) + geom_line() + geom_point() + theme_bw() +
      geom_hline(yintercept = 0) + ggtitle(paste("Spline Difference", biomarker, popname))
    print(p0)
  }
  
  robust <- selected %>% group_by(term) %>% summarise_all(median)
  robust$nA <- sum(rowsA)
  robust$nB <- sum(rowsB)
  robust$biomarker <- biomarker
  robust$model <- "prior_disease"
  return(robust)
}

## Faster version of testDiffSpline: same inputs, statistics and output, 2 GAM fits instead of 200.
# Each group is fitted once on all of its rows. The refits on random subsamples (fraction f = subsample)
# are replaced by draws from the fit's posterior: given the data, a subsample estimate varies around the
# full-data estimate with covariance (1/f - 1) * Vp, and its standard errors are 1/sqrt(f) times larger.
# The simulated noise, the GAM on the difference curve and the median summaries follow testDiffSpline.
# fitB: optional control-group fit from getGamFit(biomarker, rowsubset = rowsB & <has data>, returnFit = TRUE),
# so a caller comparing several case groups against the same controls fits the controls only once.
testDiffSplineFast <- function(biomarker, popname, rowsA, rowsB = NA, showplot = FALSE, adjustforSex = TRUE,
                               nIter = 100, subsample = 0.8, fitB = NULL) {

  if (length(rowsB) <= 1) { rowsB <- !rowsA } # if rowsB not provided, take the complement of rowsA

  hasBiomarkerData <- !is.na(X[, biomarker]) # check we have data for this biomarker

  # make sure we have some data in both groups
  rowsA <- rowsA & hasBiomarkerData
  rowsB <- rowsB & hasBiomarkerData
  if (sum(rowsA) < 10 | sum(rowsB) < 10) { return(data.frame()) } # too few observations

  fitA <- getGamFit(biomarker, rowsubset = rowsA, adjustforSex = adjustforSex, returnFit = TRUE)
  if (is.null(fitB)) {
    fitB <- getGamFit(biomarker, rowsubset = rowsB, adjustforSex = adjustforSex, returnFit = TRUE)
  }
  if (!inherits(fitA, "gam") || !inherits(fitB, "gam")) { return(data.frame()) } # a fit failed

  # 100 time points over the hours both groups cover (plot.gam's grid in testDiffSpline)
  x <- seq(max(min(fitA$model$hourminute), min(fitB$model$hourminute)),
           min(max(fitA$model$hourminute), max(fitB$model$hourminute)), length.out = 100)

  # spline curve on the grid: point estimate, nIter subsample-like draws, and subsample-sized se.
  # se is 2 standard errors, matching the se that plot.gam returns to testDiffSpline.
  drawCurves <- function(fit) {
    idx <- fit$smooth[[1]]$first.para:fit$smooth[[1]]$last.para
    newdata <- fit$model[rep(1, length(x)), ]
    newdata$hourminute <- x
    Xs <- predict(fit, newdata = newdata, type = "lpmatrix")[, idx]
    Vs <- fit$Vp[idx, idx]
    draws <- mgcv::rmvn(nIter, fit$coefficients[idx], (1 / subsample - 1) * Vs)
    list(est = drop(Xs %*% fit$coefficients[idx]),
         fit = Xs %*% t(draws),
         se = 2 * sqrt(rowSums((Xs %*% Vs) * Xs) / subsample))
  }
  curvesA <- drawCurves(fitA)
  curvesB <- drawCurves(fitB)

  # the difference-curve GAM always uses the same time points, so its basis is set up once
  # and refitted to each simulated Y (identical results to fitting Y ~ s(x) from scratch)
  G <- mgcv::gam(Y ~ s(x), data = data.frame(x, Y = numeric(length(x))), fit = FALSE)

  selected <- data.frame()
  for (itr in 1:nIter) {
    fitAi <- curvesA$fit[, itr]
    fitBi <- curvesB$fit[, itr]

    # Simulate data based on fit variability
    Y <- rnorm(length(x), mean = fitAi, sd = curvesA$se * 2) - rnorm(length(x), mean = fitBi, sd = curvesB$se * 2)
    G$y <- Y
    G$mf$Y <- Y
    lmfit <- broom::tidy(mgcv::gam(G = G))
    lmfit$sumDelta <- sum(abs(fitAi - fitBi))
    lmfit$noDiseaseRange <- max(fitBi) - min(fitBi) # popB = Controls
    lmfit$DiseaseRange <- max(fitAi) - min(fitAi)   # popA = Cases
    lmfit$maxDiff <- max(abs(fitAi - fitBi))
    selected <- bind_rows(selected, lmfit)
  }

  if (showplot) {
    plotdf <- data.frame(x = x, delta = curvesB$est - curvesA$est)   # difference in fit between B and A
    p0 <- ggplot(plotdf, aes(x, delta)) + geom_line() + geom_point() + theme_bw() +
      geom_hline(yintercept = 0) + ggtitle(paste("Spline Difference", biomarker, popname))
    print(p0)
  }

  robust <- selected %>% group_by(term) %>% summarise_all(median)
  robust$nA <- sum(rowsA)
  robust$nB <- sum(rowsB)
  robust$biomarker <- biomarker
  robust$model <- "prior_disease"
  return(robust)
}

# Wrappers below take fast = TRUE to use testDiffSplineFast; fast = FALSE runs the original testDiffSpline

# Wrapper to test spline differences between sexes in healthy individuals
testDiffSplineSex <- function(biomarker, fast = TRUE) {
  rowsubsetF <- X$sex == "Female" & X$health_rating %in% c("Good", "Excellent")
  rowsubsetM <- X$sex == "Male" & X$health_rating %in% c("Good", "Excellent")
  diffSpline <- if (fast) testDiffSplineFast else testDiffSpline
  out <- diffSpline(biomarker,
    popname = "Female",
    rowsA = rowsubsetF, rowsB = rowsubsetM, adjustforSex = FALSE
  )
  out$model <- "sexcheck"
  return(out)
}

# Wrapper to test spline differences between health ratings (Poor vs Good)
testDiffSplineHealth <- function(biomarker, fast = TRUE) {
  rowsubsetP <- X$health_rating %in% c("Fair", "Poor")
  rowsubsetG <- X$health_rating %in% c("Good", "Excellent")
  diffSpline <- if (fast) testDiffSplineFast else testDiffSpline
  out <- diffSpline(biomarker, popname = "Poor", rowsubsetP, rowsubsetG)
  out$model <- "health_check"
  return(out)
}

# Wrapper to test spline differences between sexes in healthy individuals
testDiffSplineWomen <- function(biomarker, fast = TRUE) {
  rowsubsetA <- X$sex == "Female" & X$health_rating %in% c("Good", "Excellent") & X$post_menopause == FALSE
  rowsubsetB <- X$sex == "Female" & X$health_rating %in% c("Good", "Excellent") & X$post_menopause == TRUE
  diffSpline <- if (fast) testDiffSplineFast else testDiffSpline
  out <- diffSpline(biomarker,
                    popname = "FemalePrePostMenopause",
                    rowsA = rowsubsetA, rowsB = rowsubsetB, adjustforSex = FALSE
  )
  out$model <- "FemalePrePostMenopause"
  return(out)
}
#testDiffSplineWomen("oestradiol")
#testDiffSplineWomen("testosterone")

# Fit a sex-adjusted diurnal spline for one chronotype group
# (morning, more_morning, more_evening, evening).
# Requires SLEEP (with recoded chronotype factor) and X in globalenv.
# Chronotype is constant within a group, so it is dropped from the covariates.
fitChronotype <- function(biomarker, chronotypeLevel) {
  ids <- SLEEP %>% filter(chronotype == chronotypeLevel) %>% pull(individual_id)
  tmp <- getGamFit(biomarker, rowsubset = X$individual_id %in% ids, adjustforSex = TRUE, predict = FALSE, dropCovariates = "chronotype")
  tmp$population <- chronotypeLevel
  return(tmp)
}

# Aggregate GAM fits across all four chronotype groups (mirrors gamwrap structure)
gamwrapChronotype <- function(biomarker) {
  chronotypes <- c("morning", "more_morning", "more_evening", "evening")
  dplyr::bind_rows(lapply(chronotypes, function(level) fitChronotype(biomarker, level)))
}

# Test whether chronotype modifies the diurnal spline shape of a biomarker.
# Joins SLEEP (must exist in globalenv) to X on individual_id, then fits:
#   bm ~ s(hourminute, k=24) + s(hourminute, by=chronotype_num, k=24) +
#        chronotype_num + COVARIATES (without the covariate coding of chronotype)
# The by-variable smooth captures the chronotype-linear gradient in spline shape.
# Returns broom::tidy() of the fit with biomarker column appended.
testChronotypeSpline <- function(biomarker) {
  hasBiomarkerData <- !is.na(X[, biomarker])
  if (sum(hasBiomarkerData) < 10) return(data.frame())

  tmpdf <- X %>%
    select(individual_id, bm = all_of(biomarker), hourminute, all_of(setdiff(COVARIATES, "chronotype"))) %>%
    left_join(SLEEP %>% select(individual_id, chronotype), by = "individual_id") %>%
    filter(!is.na(chronotype) & hasBiomarkerData) %>%
    mutate(chronotype_num = as.numeric(chronotype))  # 1=morning … 4=evening

  if (nrow(tmpdf) < 10) return(data.frame())

  formula_str <- paste(
    "bm ~ s(hourminute, k = 24, bs = 'ps') + s(hourminute, by = chronotype_num, k = 24, bs = 'ps') + chronotype_num +",
    covariateTerms(drop = "chronotype")
  )

  tryCatch({
    fit <- mgcv::bam(as.formula(formula_str), data = tmpdf, gamma = 1)
    out <- broom::tidy(fit)
    out$n <- fit$df.null
    out$biomarker <- biomarker
    return(out)
  }, error = function(e) {
    print(e)
    return(data.frame())
  })
}

# Plot diurnal spline for morning vs evening chronotypes for a single biomarker
plotChronotypeSpline <- function(biomarker) {
  hasBiomarkerData <- !is.na(X[, biomarker])
  tmpdf <- X %>%
    select(individual_id, bm = all_of(biomarker), hourminute, all_of(setdiff(COVARIATES, "chronotype"))) %>%
    left_join(SLEEP %>% select(individual_id, chronotype), by = "individual_id") %>%
    filter(!is.na(chronotype) & hasBiomarkerData)

  morning <- tmpdf %>% filter(chronotype %in% c("morning"))
  evening <- tmpdf %>% filter(chronotype %in% c("evening"))

  fit_m <- getGamFit(biomarker, rowsubset = tmpdf$individual_id %in% morning$individual_id, plotdata = TRUE, adjustforSex = TRUE, dropCovariates = "chronotype")
  fit_e <- getGamFit(biomarker, rowsubset = tmpdf$individual_id %in% evening$individual_id, plotdata = TRUE, adjustforSex = TRUE, dropCovariates = "chronotype")
  fit_m$chronotype <- "morning"
  fit_e$chronotype <- "evening"

  mycolors <- c("#E41A1C", "#377EB8")
  ggplot(bind_rows(fit_m, fit_e), aes(x, fit, fill = chronotype)) +
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.5) +
    scale_fill_manual(values = mycolors) +
    geom_hline(yintercept = 0) +
    theme_bw(base_size = 14) +
    xlab("Hour") + ylab("Relative Abundance") +
    ggtitle(biomarker) +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.95),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))
}

## utility functions for pairs plots
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if (missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.lm <- function(x, y, col = par("col"), bg = NA, pch = par("pch"),
                     cex = 1, col.smooth = 2, ...) {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) {
    abline(MASS::rlm(y[ok] ~ x[ok]), col = col.smooth, ...)
  }
}

## utility functions for comparing fits

# Correlation between row i of two tables (first column = biomarker, remaining columns = betas)
pair_cor <- function(df1, df2, i) {
  cor(as.double(df1[i, -1]), as.double(df2[i, -1]))
}

# Fit y ~ x (per group, if df is grouped) and add each row's relative deviation:
#   reference = "fit":      |y - predicted| / predicted   (distance from the fitted line)
#   reference = "identity": |y - x| / predicted           (distance from the y = x line)
# Also adds the fitted slope and R-squared (uncentered when through_origin = TRUE).
fit_deviation <- function(df, x, y, through_origin = FALSE, reference = c("fit", "identity")) {
  reference <- match.arg(reference)
  formula_str <- paste(y, "~", x, if (through_origin) "- 1" else "")
  df %>%
    filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>%
    group_modify(~ {
      fit <- lm(as.formula(formula_str), data = .x)
      ref <- if (reference == "fit") fitted(fit) else .x[[x]]
      .x %>% mutate(
        predicted = fitted(fit),
        rel_dev   = abs(.data[[y]] - ref) / predicted,
        slope     = coef(fit)[[x]],
        r2        = summary(fit)$r.squared
      )
    }) %>%
    ungroup()
}

### Disease assocaitons functions

# Time to first event (years from first visit) for the given ICD10 codes.
# Requires UKBB, PREV_CONDITIONS, EVENTS and CENSOR_DATE in globalenv.
# Individuals without an event are censored at death or censor_date.
getAllDiseaseEvents <- function(icd10str, censor_date = CENSOR_DATE) {
  tmp <- UKBB %>%
    left_join(PREV_CONDITIONS, by = "individual_id") %>%
    left_join(EVENTS %>% filter(icd10 %in% c(icd10str)), by = "individual_id")

  # disease events
  tmp$timetoevent <- tmp$date - as.Date(tmp$visit_date)

  # experienced event post first visit
  tmp$event <- ifelse(!is.na(tmp$timetoevent), 1, 0)

  # consider death events
  deathevents <- EVENTS %>%
    filter(icd10 == "Death") %>%
    select(individual_id, censor_date = date)

  tmp <- tmp %>% left_join(deathevents)
  tmp$censor_date[is.na(tmp$censor_date)] <- censor_date

  # censor events
  tmp$timetoevent <- ifelse(!is.na(tmp$timetoevent),
    tmp$timetoevent,
    tmp$censor_date - as.Date(tmp$first_visit_date)
  )

  tmp$timetoevent <- as.double(tmp$timetoevent / 365.25)

  return(tmp)
}

getAnyTimeEvents <- function(icd10) {
  tmpdf <- getAllDiseaseEvents(icd10) %>%
    group_by(individual_id, instance_id) %>%
    summarise(timetoevent = min(timetoevent), event = max(event))
  
  anyevents <- X %>%
    select(individual_id, instance_id) %>%
    left_join(tmpdf) %>%
    select(individual_id, timetoevent, event)
  return(anyevents)
}

getPriorDiseaseEvents <- function(icd10, maxtimetoevent = 2) {
  anytimeevents <- getAnyTimeEvents(icd10)
  diseaseRows <- anytimeevents$event & anytimeevents$timetoevent < maxtimetoevent
  return(diseaseRows)
}

getFutureEvents <- function(icd10, mintimetoevent = 2) {
  anytimeevents <- getAnyTimeEvents(icd10)
  diseaseRows <- anytimeevents$event & anytimeevents$timetoevent >= mintimetoevent
  return(diseaseRows)
}

getTimeToEvent <- function(icd10) {
  anytimeevents <- getAnyTimeEvents(icd10)
  return(anytimeevents)
}

# fast = TRUE uses testDiffSplineFast; fast = FALSE runs the original testDiffSpline.
# Biomarkers run in parallel on `cores` forked workers. In fast mode each biomarker's never-diseased
# control group is fitted once and shared by the prior and future comparisons.
modelPastFutureEvents <- function(icd10, fast = TRUE, cores = 6) {
  anytimeevents <- getAnyTimeEvents(icd10) # one event lookup for both case definitions
  diseaseRowsPast <- anytimeevents$event & anytimeevents$timetoevent < 2     # as getPriorDiseaseEvents(icd10, 2)
  diseaseRowsFuture <- anytimeevents$event & anytimeevents$timetoevent >= 2  # as getFutureEvents(icd10, 2)
  diseaseRowsNever <- !diseaseRowsPast & !diseaseRowsFuture

  fitBiomarker <- function(bm) {
    if (!fast) {
      return(list(past = testDiffSpline(bm, icd10, rowsA = diseaseRowsPast, rowsB = diseaseRowsNever),
                  future = testDiffSpline(bm, icd10, rowsA = diseaseRowsFuture, rowsB = diseaseRowsNever)))
    }
    fitB <- getGamFit(bm, rowsubset = diseaseRowsNever & !is.na(X[, bm]), adjustforSex = TRUE, returnFit = TRUE)
    list(past = testDiffSplineFast(bm, icd10, rowsA = diseaseRowsPast, rowsB = diseaseRowsNever, fitB = fitB),
         future = testDiffSplineFast(bm, icd10, rowsA = diseaseRowsFuture, rowsB = diseaseRowsNever, fitB = fitB))
  }

  print(icd10)
  RESULTS <- parallel::mclapply(BIOMARKER_NAMES$biomarker, fitBiomarker, mc.cores = cores)
  failed <- vapply(RESULTS, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    warning(icd10, ": ", sum(failed), " biomarkers failed: ", paste(BIOMARKER_NAMES$biomarker[failed], collapse = ", "))
  }
  RESULTS <- RESULTS[!failed]

  PASTFITS <- bind_rows(lapply(RESULTS, `[[`, "past"))
  PASTFITS$icd10 <- icd10
  PASTFITS$model <- "prior_disease"
  saveRDS(PASTFITS, file = paste0("outputs/", icd10, "_splinefit.Rds"))

  FITS <- bind_rows(lapply(RESULTS, `[[`, "future"))
  FITS$icd10 <- icd10
  FITS$model <- "future_disease"
  saveRDS(FITS, file = paste0("outputs/", icd10, "_future_splinefit.Rds"))
}

outlierTruncation = function(df, colname, lower_quantile = 0.01, upper_quantile = 0.99) {
  lower_bound <- quantile(df[[colname]], probs = lower_quantile, na.rm = TRUE)
  upper_bound <- quantile(df[[colname]], probs = upper_quantile, na.rm = TRUE)
  
  df[[colname]] <- pmin(pmax(df[[colname]], lower_bound), upper_bound)
  return(df)
}



testDiseaseInteraction <- function(eventdf, biomarker, icd10,verbose=F) {
  hasBiomarkerData <- !is.na(X[, biomarker])
  if (sum(hasBiomarkerData) < 10) return(data.frame())
  
  tmpdf = X %>% select("individual_id", all_of(c(SURVIVAL_COVARIATES, biomarker)));
  tmpdf =  tmpdf %>% left_join(eventdf, by = "individual_id")
  tmpdf = outlierTruncation(tmpdf, biomarker) # Truncate outliers to reduce influence on Cox fit
  formula_str = sprintf("Surv(timetoevent,event) ~ pspline(%s*hourminute) + %s", biomarker, covariateTerms(SURVIVAL_COVARIATES));
  
  tryCatch({
    bmfit=coxph(as.formula(formula_str), data=tmpdf)
    if(verbose) { print(summary(bmfit)); } 
    #print(bmfit)
    out=broom::tidy(bmfit) %>% mutate(biomarker=biomarker, icd10=icd10);
    #print(out);
    return(out)
  }, error = function(e) {
    print(e)
    return(data.frame())
  })
}

# Forest plot of model coefficients (estimate +/- 2 SE), keeping the model's term order
forestPlotGlmBroomUnOrder <- function(glm_model, sig_cut = 0.01, termgrep = NA) {
  glm_coef <- broom::tidy(glm_model) %>% filter(!term == "(Intercept)")
  if (!is.na(termgrep)) { glm_coef <- glm_coef %>% filter(grepl(termgrep, term)) }
  glm_coef$color <- ifelse(glm_coef$p.value < sig_cut, "black", "grey")
  glm_coef <- glm_coef %>% filter(p.value < sig_cut)

  p0 <- ggplot(glm_coef, aes(estimate, term,
                             xmin = estimate - std.error * 2,
                             xmax = estimate + std.error * 2)) +
    geom_point(shape = 18, fill = "black", color = "black", size = 5) +
    geom_errorbarh(height = .2) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) +
    ylab("") + theme_bw(base_size = 14)
  return(p0)
}
