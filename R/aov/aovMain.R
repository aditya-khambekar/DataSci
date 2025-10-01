aov.residualsnormal <- function(formula, data){
  data.aov <- aov(formula, data = data)
  qqnorm(data.aov$residuals)
  qqline(data.aov$residuals)
}

aov.residualsd <- function(formula, data) {
  data.aov <- aov(formula, data = data)
  residual.df <- data.frame(
    Group = model.frame(data.aov)[[all.vars(formula)[2]]],
    Residual = data.aov$residuals
  )
  result <- aggregate(Residual ~ Group, data = residual.df, FUN = sd)
  colnames(result) <- c("Group", "SD")

  return(result)
}

aov.df <- function(n.obs, n.groups) {
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

aov.dataset.summary <- function(formula, data){
  exp = all.vars(formula)[2]
  res = all.vars(formula)[1]

  means = tapply(data[[res]], data[[exp]], mean)
  sds = tapply(data[[res]], data[[exp]], sd)
  lengths = tapply(data[[res]], data[[exp]], length)

  return(cbind(means, sds, lengths))
}

ak.aov <- function(formula, data) {
    return(ak.aov$new(aov(formula, data)
}
