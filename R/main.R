loadData <- function(dataset) {
  if (!requireNamespace("Stat2Data", quietly = TRUE)) {
    stop("Package 'Stat2Data' is not installed. Please install it first.")
  }

  available_data <- data(package = "Stat2Data")$results[, "Item"]
  
  exact_match <- available_data[tolower(available_data) == tolower(dataset)]
  
  if (length(exact_match) == 1) {
    data(list = exact_match, package = "Stat2Data", envir = parent.frame())
    return(get(exact_match, envir = parent.frame()))
  }
  
  distances <- adist(tolower(dataset), tolower(available_data))
  best_match_index <- which.min(distances)
  best_match <- available_data[best_match_index]
  best_distance <- distances[best_match_index]
  
  
  data(list = best_match, package = "Stat2Data", envir = parent.frame())
}

residualplot <- function(formula, data, point.col = "blue", line.col = "grey",
                             line.lty = 2, ...) {
  library("lattice")
  exp.var <- all.vars(formula)[2]

  data.aov <- aov(formula, data = data)

  xyplot(data.aov$residuals ~ data[[exp.var]], data = data,
         panel = function(x, y, col, ...) {
           panel.xyplot(x, y, col = point.col, ...)
           panel.abline(h = 0, col = line.col, lty = line.lty)
         },
         xlab = exp.var,
         ylab = "Residuals",
  )
}

xyplot <- function(formula, data, point.col = "blue", line.col = "grey",
                       line.lty = 2, ...) {
    library("lattice")

  response.var <- all.vars(formula)[1]

  grand.mean <- mean(data[[response.var]], na.rm = TRUE)

  xyplot(formula, data = data,
         panel = function(x, y, col, ...) {
           panel.xyplot(x, y, col = point.col, ...)
           panel.abline(h = grand.mean, col = line.col, lty = line.lty)
         },
         ...
  )
}

residualsnormal <- function(formula, data){
  data.aov <- aov(formula, data = data)
  qqnorm(data.aov$residuals)
  qqline(data.aov$residuals)
}

residualsd <- function(formula, data) {
  data.aov <- aov(formula, data = data)
  residual.df <- data.frame(
    Group = model.frame(data.aov)[[all.vars(formula)[2]]],
    Residual = data.aov$residuals
  )
  result <- aggregate(Residual ~ Group, data = residual.df, FUN = sd)
  colnames(result) <- c("Group", "SD")

  return(result)
}

df <- function(n.obs, n.groups) {
  data.frame(
    Type = c("Groups", "Error"),
    df = c(n.groups - 1, n.obs - n.groups)
  )
}

aov.summary.fillTable <- function(df.groups = NA, df.error = NA, df.total = NA,
                                  ss.groups = NA, ss.error = NA, ss.total = NA,
                                  ms.groups = NA, ms.error = NA,
                                  f.stat = NA, n.obs = NA, n.groups = NA) {

  if (!is.na(n.obs) && !is.na(n.groups)) {
    if (is.na(df.groups)) df.groups <- n.groups - 1
    if (is.na(df.error)) df.error <- n.obs - n.groups
    if (is.na(df.total)) df.total <- n.obs - 1
  } else {
    if (sum(!is.na(c(df.groups, df.error, df.total))) >= 2) {
      if (is.na(df.groups) && !is.na(df.error) && !is.na(df.total)) {
        df.groups <- df.total - df.error
      }
      if (is.na(df.error) && !is.na(df.groups) && !is.na(df.total)) {
        df.error <- df.total - df.groups
      }
      if (is.na(df.total) && !is.na(df.groups) && !is.na(df.error)) {
        df.total <- df.groups + df.error
      }
    }
  }

  if (sum(!is.na(c(ss.groups, ss.error, ss.total))) >= 2) {
    if (is.na(ss.groups) && !is.na(ss.error) && !is.na(ss.total)) {
      ss.groups <- ss.total - ss.error
    }
    if (is.na(ss.error) && !is.na(ss.groups) && !is.na(ss.total)) {
      ss.error <- ss.total - ss.groups
    }
    if (is.na(ss.total) && !is.na(ss.groups) && !is.na(ss.error)) {
      ss.total <- ss.groups + ss.error
    }
  }

  if (!is.na(ss.groups) && !is.na(df.groups) && is.na(ms.groups)) {
    ms.groups <- ss.groups / df.groups
  }
  if (!is.na(ss.error) && !is.na(df.error) && is.na(ms.error)) {
    ms.error <- ss.error / df.error
  }

  if (is.na(ss.groups) && !is.na(ms.groups) && !is.na(df.groups)) {
    ss.groups <- ms.groups * df.groups
  }
  if (is.na(ss.error) && !is.na(ms.error) && !is.na(df.error)) {
    ss.error <- ms.error * df.error
  }

  if (!is.na(ms.groups) && !is.na(ms.error) && is.na(f.stat)) {
    f.stat <- ms.groups / ms.error
  }

  anova.table <- data.frame(
    Source = c("Groups", "Error", "Total"),
    DF = c(df.groups, df.error, df.total),
    SS = c(ss.groups, ss.error, ss.total),
    MS = c(ms.groups, ms.error, NA),
    F = c(f.stat, NA, NA)
  )

  return(anova.table)
}

dataset.summary <- function(formula, data){
  exp = all.vars(formula)[2]
  res = all.vars(formula)[1]

  means = tapply(data[[res]], data[[exp]], mean)
  sds = tapply(data[[res]], data[[exp]], sd)
  lengths = tapply(data[[res]], data[[exp]], length)

  return(cbind(means, sds, lengths))
}

aov.residualplot <- function(aov.obj, point.col = "blue", line.col = "grey",
                                  line.lty = 2, ...) {
  library("lattice")

  formula <- formula(aov.obj)
  data <- model.frame(aov.obj)
  exp.var <- all.vars(formula)[2]

  xyplot(aov.obj$residuals ~ data[[exp.var]], data = data,
         panel = function(x, y, col, ...) {
           panel.xyplot(x, y, col = point.col, ...)
           panel.abline(h = 0, col = line.col, lty = line.lty)
         },
         xlab = exp.var,
         ylab = "Residuals",
  )
}

aov.xyplot <- function(aov.obj, point.col = "blue", line.col = "grey",
                            line.lty = 2, ...) {
  library("lattice")

  formula <- formula(aov.obj)
  data <- model.frame(aov.obj)
  response.var <- all.vars(formula)[1]

  grand.mean <- mean(data[[response.var]], na.rm = TRUE)

  xyplot(formula, data = data,
         panel = function(x, y, col, ...) {
           panel.xyplot(x, y, col = point.col, ...)
           panel.abline(h = grand.mean, col = line.col, lty = line.lty)
         },
         ...
  )
}

aov.residualsnormal <- function(aov.obj){
  qqnorm(aov.obj$residuals)
  qqline(aov.obj$residuals)
}

aov.residualsd <- function(aov.obj) {
  formula <- formula(aov.obj)
  residual.df <- data.frame(
    Group = model.frame(aov.obj)[[all.vars(formula)[2]]],
    Residual = aov.obj$residuals
  )
  result <- aggregate(Residual ~ Group, data = residual.df, FUN = sd)
  colnames(result) <- c("Group", "SD")

  return(result)
}

aov.dataset.summary <- function(aov.obj){
  formula <- formula(aov.obj)
  data <- model.frame(aov.obj)

  exp = all.vars(formula)[2]
  res = all.vars(formula)[1]

  means = tapply(data[[res]], data[[exp]], mean)
  sds = tapply(data[[res]], data[[exp]], sd)
  lengths = tapply(data[[res]], data[[exp]], length)

  return(cbind(means, sds, lengths))
}

TukeyHSD <- function(aov.obj, which = NULL, conf.level = 0.95, ordered = FALSE) {
  if (!inherits(aov.obj, "aov"))
    stop("This function requires an object of class 'aov'.")

  mf <- model.frame(aov.obj)
  response <- mf[[1]]

  if (is.null(which))
    which <- names(mf)[-1]

  results <- list()

  for (term in which) {
    group <- mf[[term]]
    if (!is.factor(group))
      stop(paste("Term", term, "is not a factor."))

    means <- tapply(response, group, mean)
    sizes <- tapply(response, group, length)

    mse <- deviance(aov.obj) / df.residual(aov.obj)
    df <- df.residual(aov.obj)
    nlev <- length(levels(group))

    combs <- combn(levels(group), 2)
    ncomp <- ncol(combs)
    comp.names <- apply(combs, 2, function(x) paste(x, collapse = "-"))

    diffs <- numeric(ncomp)
    lwr <- numeric(ncomp)
    upr <- numeric(ncomp)
    pvals <- numeric(ncomp)

    qcrit <- qtukey(conf.level, nlev, df) / sqrt(2)

    for (i in 1:ncomp) {
      g1 <- combs[1, i]
      g2 <- combs[2, i]
      diffs[i] <- means[g1] - means[g2]
      se <- sqrt(mse / 2 * (1 / sizes[g1] + 1 / sizes[g2]))
      lwr[i] <- diffs[i] - qcrit * se
      upr[i] <- diffs[i] + qcrit * se
      qval <- abs(diffs[i]) / se
      pvals[i] <- 1 - ptukey(qval * sqrt(2), nlev, df)
    }

    res <- data.frame(
      diff = diffs,
      lwr = lwr,
      upr = upr,
      p.adj = pvals,
      row.names = comp.names
    )

    results[[term]] <- res
  }

  # Assign attributes that the built-in plot() method expects
  class(results) <- "TukeyHSD"
  attr(results, "conf.level") <- conf.level
  attr(results, "ordered") <- ordered

  return(results)
}

aov2.df = function(nA, nB, nObs) {
  if (nA * nB == nObs) return (c(nA -1, nB - 1, nA * nB))
  return (c(nA - 1, nB - 1, nA * nB * (nObs -1)))
}

anscombe.plot = function(formula, data, xlim = 1000, ylim = 1000,
                         main = "", xlab = "", ylab = "", legendstr = "topleft", verbose = FALSE, roundSlope = 3, flip = FALSE) {
  v = all.vars(formula)
  
  if (main == "") main = paste(v[1], "Blocked By", v[2])
  
  plot(0, 0, type = 'n',
       xlim = c(0, xlim), ylim = c(0, ylim),
       xlab = xlab, ylab = ylab, main = main)
  abline(0, 1, col = 'grey')
  
  c = combn(u <- sort(unique(data[[v[2]]])), 2)
  colors = rainbow(ncol(c))
  by = !flip
  
  for (i in 1:ncol(c)) {
    abline(coe <- coefficients(lm(
      data[[v[1]]][data[[v[2]]] == c[ifelse(by, 1, 2), i]] ~
        data[[v[1]]][data[[v[2]]] == c[ifelse(by, 2, 1), i]]
    )),
    col = colors[i])
    points(data[[v[1]]][data[[v[2]]] == c[ifelse(by, 1, 2), i]] ~
           data[[v[1]]][data[[v[2]]] == c[ifelse(by, 2, 1), i]],
           col = colors[i], pch = i)
    if (verbose) print(paste("Slope of", paste(u[c[1, i]], u[c[2, i]], sep = "-"), "=", round(coe[2], roundSlope)))
  }
  
  legend(legendstr,
         legend = apply(c, 2, function(x) paste(u[x[1]], u[x[2]], sep = "-")),
         col = colors, pch = 1:ncol(c), lty = 1)
}
