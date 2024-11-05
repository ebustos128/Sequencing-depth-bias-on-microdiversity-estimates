#Input files
ins <- read.csv("nucleotide.diversity.table.csv", sep = "\t")
tad <- read.csv("coverage.table.csv", sep = "\t")

ins$Genomospecies <- tolower(ins$Genomospecies)
ins[ins == 0] <- NA
tad[tad == 0] <- NA

gsp <- unique(tad$Genomospecies)
mg <- colnames(tad)[-c(1:2)]

#You can adjust the sequencing depth coverage threshold
depth_limit <- 5

x.breaks <- 10^seq(-6, 3, by = 0.25)
y.val <- vector("list", length(x.breaks))
y.errors <- vector("list", length(x.breaks))
y.n <- rep(0, length(x.breaks))

# Select all data and estimate errors
for (j in gsp) {
  for (i in mg) {
    x <- tad[tad$Genomospecies == j, i]
    y <- ins[ins$Genomospecies == j, i]
    s <- !is.na(x) & !is.na(y)
    if (!any(s)) next
    x <- x[s]
    y <- y[s] / tail(y, n = 1, na.rm = TRUE)
    
    for (k in seq_along(x)) {
      x.k <- head(which(x[k] < x.breaks), n = 1)
      if (length(x.k) > 0 && !is.na(x.k) && x.k <= length(y.val)) {
        y.val[[x.k]] <- c(y.val[[x.k]], y[k])
        
        y.expected <- mean(y.val[[x.k]], na.rm = TRUE)
        y.error <- abs(y[k] - y.expected) / y.expected * 100
        y.errors[[x.k]] <- c(y.errors[[x.k]], y.error)
        
        y.n[x.k] <- y.n[x.k] + 1
      }
    }
}

# Estimate the Average Estimated Error and the Standard Error for the filtered data
error_mean <- sapply(y.errors, function(x) if (length(x) > 0) mean(x, na.rm = TRUE) else NA)
error_sd <- sapply(y.errors, function(x) if (length(x) > 0) sd(x, na.rm = TRUE) else NA)
n <- sapply(y.errors, length)
error_se <- error_sd / sqrt(n)  # Error estándar

# Filter data for coverages >= 5X or the desired threshold
y.val_filtered <- vector("list", length(x.breaks))
y.errors_filtered <- vector("list", length(x.breaks))
y.n_filtered <- rep(0, length(x.breaks))

for (i in seq_along(x.breaks)) {
  if (x.breaks[i] >= depth_limit) {
    y.val_filtered[[i]] <- y.val[[i]]
    y.errors_filtered[[i]] <- y.errors[[i]]
    y.n_filtered[i] <- y.n[i]
  } else {
    y.val_filtered[[i]] <- NULL
    y.errors_filtered[[i]] <- NULL
    y.n_filtered[i] <- 0
  }
}

# Estimate the Average Estimated Error and the Standard Error
error_mean_filtered <- sapply(y.errors_filtered, function(x) if (length(x) > 0) mean(x, na.rm = TRUE) else NA)
error_sd_filtered <- sapply(y.errors_filtered, function(x) if (length(x) > 0) sd(x, na.rm = TRUE) else NA)
n_filtered <- sapply(y.errors_filtered, length)
error_se_filtered <- error_sd_filtered / sqrt(n_filtered)  # Error estándar

# Plot the Average Estimated Error and the Standard Error
plot_with_error_bars <- function(x, y_mean, y_se, title, xlab, ylab) {
  # Remove NA for x and y_mean
  valid <- !is.na(y_mean) & !is.infinite(y_mean)
  x <- x[valid]
  y_mean <- y_mean[valid]
  y_se <- y_se[valid]
  
  if (length(x) > 0) {
    plot(
      x, y_mean, log = "x", type = "o", las = 1,
      ylab = ylab, xlab = xlab,
      ylim = range(c(y_mean - 1.96 * y_se, y_mean + 1.96 * y_se), na.rm = TRUE)
    )
    arrows(
      x, y_mean - 1.96 * y_se, 
      x, y_mean + 1.96 * y_se, 
      angle = 90, code = 3, length = 0.1, col = "blue"
    )
    abline(v = depth_limit, col = "darkred")
    title(title)
  } else {
    plot.new()
    title(title)
    mtext("No data available", side = 3)
  }
}

# Plot of the the Average Estimated Error and the Standard Error for coverages >= 5X or the desired threshold
plot_with_error_bars(x.breaks, error_mean_filtered, error_se_filtered, "Nucleotide diversity", "Sequencing Depth", "Average Absolute Error (%)")
