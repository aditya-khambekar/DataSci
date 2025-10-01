# Standalone Functions (Minor Fix in aov.residualsd)

loadData <- function(dataset) {
  library("Stat2Data")
  data(list = dataset)
}

aov.residualplot <- function(formula, data, point.col = "blue", line.col = "grey",
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

aov.xyplot <- function(formula, data, point.col = "blue", line.col = "grey",
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

aov.residualsnormal <- function(formula, data){
  data.aov <- aov(formula, data = data)
  qqnorm(data.aov$residuals)
  qqline(data.aov$residuals)
}

# FIX: Changed 'residualsd' to 'aov.residualsd' (consistency)
aov.residualsd <- function(formula, data) {
  data.aov <- aov(formula, data = data)
  # NOTE: The formula variables should be extracted from the model frame, not the raw data
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

# FIX: The R6 (or setRefClass) definition must be defined *before* the generator function.
# Also, the generator function must return the NEW object.
ak.aov = setRefClass(
    "ak.aov", # Recommended to name the class for clarity
    fields = list(aov = "aov"), # Declare field type
    methods = list(
        residualsnormal = function(){
          qqnorm(aov$residuals)
          qqline(aov$residuals)
        },
        residualsd = function() {
          # Use the stored aov object's call for formula and model.frame for data
          residual.df <- data.frame(
            Group = model.frame(aov)[[all.vars(aov$call$formula)[2]]],
            Residual = aov$residuals
          )
          result <- aggregate(Residual ~ Group, data = residual.df, FUN = sd)
          colnames(result) <- c("Group", "SD")

          return(result)
        },
        dataset.summary = function(){
          # Correctly call aov.dataset.summary with formula and data extracted from the aov object
          return(aov.dataset.summary(aov$call$formula, model.frame(aov)))
        },
        # FIX: Missing comma after the previous method definition
        # FIX: The function needs to use the formula and data *stored in the aov object*
        xyplot = function(point.col = "blue", line.col = "grey",
                             line.lty = 2, ...){
          library("lattice")

          formula <- aov$call$formula # Get formula from the aov object
          data <- model.frame(aov)    # Get data from the aov object

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
    )
)

# FIX: The function ak.aov must be defined after the class 'ak.aov' has been defined.
# FIX: It must return a new instance of the Reference Class.
create.aov <- function(formula, data) {
    # The generator object is the one created by setRefClass, which is also named 'ak.aov'
    # We call $new on the generator object, passing the aov object to the 'aov' field.
    return(ak.aov$new(aov = aov(formula, data = data)))
}