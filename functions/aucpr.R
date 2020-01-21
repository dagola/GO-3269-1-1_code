
# Define a function that calculates the Precision-Recall AUC
measureAUCPR = function(task, model, pred, feats, extra.args) {
  if (checkmate::anyMissing(pred$data$response) || length(unique(pred$data$truth)) == 1L) {
    return(NA_real_)
  } 
  
  scores_class0 <- getPredictionProbabilities(pred)
  scores_class1 <- scores_class0
  
  weights_class0 <- as.integer(pred$data$truth == pred$task.desc$positive)
  weights_class1 <- 1 - weights_class0
  
  o0 <- order(scores_class0)
  sorted_scores_class0 <- scores_class0[o0]
  if (!is.null(weights_class0)) {
    weights_class0 <- weights_class0[o0]
  }
  o1 <- order(scores_class1)
  sorted_scores_class1 <- scores_class1[o1]
  if (!is.null(weights_class1)) {
    weights_class1 <- weights_class1[o1]
  }
  
  if (
    !is.null(sorted_scores_class1) & 
    (length(sorted_scores_class0) != length(sorted_scores_class1) | 
     suppressWarnings(sum(sorted_scores_class0 != sorted_scores_class1) > 0)) & 
    is.null(weights_class0) & 
    is.null(weights_class1)
  ) {
    weights_class0 <- c(rep(1, length(sorted_scores_class0)), rep(0, length(sorted_scores_class1)))
    sorted_scores_class0 <- c(sorted_scores_class0, sorted_scores_class1)
    o0 <- order(sorted_scores_class0)
    sorted_scores_class0 <- sorted_scores_class0[o0]
    weights_class0 <- weights_class0[o0]
    weights_class1 <- 1 - weights_class0
    sorted_scores_class1 <- sorted_scores_class0
    all_scores <- sorted_scores_class0
    all_weights_pos <- weights_class0
    all_weights_neg <- weights_class1
  } else {
    if (is.null(weights_class0)) {
      weights_class0 <- rep(1, length(sorted_scores_class0))
    }
    if (is.null(weights_class1)) {
      weights_class1 <- rep(1, length(sorted_scores_class1))
    }
    all_scores <- c(sorted_scores_class0, sorted_scores_class1)
    all_weights_pos <- c(weights_class0, rep(0, length(sorted_scores_class1)))
    all_weights_neg <- c(rep(0, length(sorted_scores_class0)), 
                         weights_class1)
  }
  
  davis_and_goadrich <- (length(sorted_scores_class0) == length(sorted_scores_class1) & 
                           suppressWarnings(sum(sorted_scores_class0 != sorted_scores_class1) == 0) & 
                           length(weights_class0) == length(weights_class1) & 
                           suppressWarnings(sum(weights_class0 != (1 - weights_class1)) == 0) &
                           sum(weights_class0 != 0 & weights_class0 != 1) == 0)
  
  auc_dg <- NA_real_
  if (davis_and_goadrich) {
    o <- order(all_scores, decreasing = T)
    all_scores <- all_scores[o]
    all_weights_pos <- all_weights_pos[o]
    all_weights_neg <- all_weights_neg[o]
    cum_weights_pos <- cumsum(all_weights_pos)
    cum_weights_neg <- cumsum(all_weights_neg)
    cum_use <- c(all_scores[-length(all_scores)] != all_scores[-1], TRUE)
    all_scores <- all_scores[cum_use]
    cum_weights_pos <- cum_weights_pos[cum_use]
    cum_weights_neg <- cum_weights_neg[cum_use]
    r_fg <- sum(all_weights_pos)
    tp <- cum_weights_pos
    fp <- cum_weights_neg
    tp_prev <- c(0, cum_weights_pos[-length(cum_weights_pos)])
    fp_prev <- c(0, cum_weights_neg[-length(cum_weights_neg)])
    h <- (fp - fp_prev)/(tp - tp_prev)
    a <- 1 + h
    b <- (fp_prev - h * tp_prev)/r_fg
    h[tp == tp_prev] <- 1
    a[tp == tp_prev] <- 1
    b[tp == tp_prev] <- 0
    min_mat <- cbind(tp_prev, tp, fp_prev, fp)
    idxs <- which(tp - tp_prev > 1 & tp/(tp + fp) != tp_prev/(tp_prev + fp_prev))
    if (length(idxs) > 0) {
      m <- matrix(min_mat[-idxs, ], ncol = ncol(min_mat))
    }
    else {
      m <- min_mat
    }
    auc_dg <- (m[, 2] - m[, 1])/r_fg * (m[, 1]/(m[, 1] + m[, 3]) + m[, 2]/(m[, 2] + m[, 4]))/2
    if (is.nan(auc_dg[1])) {
      auc_dg[1] <- (m[1, 2] - m[1, 1])/r_fg * (m[1, 2]/(m[1, 2] + m[1, 4]))
    }
    diff <- tp - tp_prev
    h2 <- (fp - fp_prev)/diff
    if (length(idxs) > 0) {
      temp_seq <- sapply(1:max(diff[idxs]), function(i) seq(0, i))
      m <- sapply(idxs, function(i) {
        x <- temp_seq[[tp[i] - tp_prev[i]]]
        prcs <- (tp_prev[i] + x)/(tp_prev[i] + x + fp_prev[i] + h2[i] * x)
        sum(1/r_fg * (prcs[-1] + prcs[-length(prcs)])/2)
      })
      auc_dg <- sum(c(auc_dg, m))
    }
    else {
      auc_dg <- sum(auc_dg)
    }
  }
  
  return(auc_dg)
  
}

# Generate the Measure object
aucpr = makeMeasure(
  id = "aucpr", name = "Area Under the Precision-Recall Curve",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = FALSE, best = 1, worst = 0,
  fun = measureAUCPR
)

ci.aucpr <- function(aucpr, n, ci.level) {
  
  eta <- log(aucpr/(1-aucpr))
  tau <- (n*aucpr*(1-aucpr))^(-0.5)
  
  q <- qnorm(1-ci.level/2)
  
  l <- eta - q*tau
  u <- eta + q*tau
  
  return(c(aucpr.lower = unname(exp(l)/(1+exp(l))), aucpr.upper = unname(exp(u)/(1+exp(u)))))
  
}
