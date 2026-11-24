# Retain evaluated terms and actual contrast matrices, including the constants
# that R records for predicting polynomial, spline, and scaled terms.
.mjoint_design <- function(frame, matrix, drop.intercept = FALSE) {
  coding <- attr(matrix, "contrasts")
  for (name in names(coding)) {
    value <- frame[[name]]
    if (is.character(value)) value <- factor(value)
    contrasts(value) <- coding[[name]]
    coding[[name]] <- contrasts(value)
  }
  list(terms = attr(frame, "terms"),
       xlevels = .getXlevels(attr(frame, "terms"), frame),
       contrasts = coding, columns = colnames(matrix),
       drop.intercept = drop.intercept)
}

.mjoint_new_design <- function(design, data) {
  frame <- model.frame(design$terms, data, xlev = design$xlevels,
                       na.action = na.fail)
  matrix <- model.matrix(design$terms, frame, contrasts.arg = design$contrasts)
  if (design$drop.intercept) {
    matrix <- matrix[, attr(matrix, "assign") != 0L, drop = FALSE]
  }
  if (!identical(colnames(matrix), design$columns)) {
    stop("The prediction design columns do not match the fitted model; refit with 'mjoint()'.")
  }
  list(frame = frame, matrix = matrix)
}

.mjoint_prediction_designs <- function(object) {
  if (!is.null(object$designs)) return(object$designs)

  # Older fits can be used only if reconstruction agrees with the matrices
  # actually used for estimation, including those at the event times.
  if (is.null(object$dmats$l$Zik.list) || is.null(object$dmats$z$Zi.fail)) {
    stop("The fitted random-effects basis cannot be verified; refit with 'mjoint()' before predicting.")
  }
  designs <- list(fixed = list(), random = list())
  same.values <- function(x, y) {
    isTRUE(all.equal(unname(x), unname(y), check.attributes = FALSE,
                     tolerance = 1e-10)) && identical(dim(x), dim(y))
  }
  for (k in seq_len(object$dims$K)) {
    frame <- model.frame(object, k = k)
    fixed.matrix <- model.matrix(object, k = k)
    designs$fixed[[k]] <- .mjoint_design(frame, fixed.matrix)
    if (!same.values(fixed.matrix, do.call(rbind, object$dmats$l$Xik.list[[k]]))) {
      stop("The reconstructed fixed-effects design differs from the fitted model; refit with 'mjoint()' before predicting.")
    }
    formula <- nlme::splitFormula(object$formLongRandom[[k]], "|")[[1]]
    frame <- model.frame(formula, object$data[[k]], na.action = na.fail)
    coding <- object$lfit[[k]]$contrasts
    coding <- coding[intersect(names(coding), names(frame))]
    matrix <- model.matrix(attr(frame, "terms"), frame, contrasts.arg = coding)
    designs$random[[k]] <- .mjoint_design(frame, matrix)
    fitted.matrix <- do.call(rbind, object$dmats$l$Zik.list[[k]])
    if (!same.values(matrix, fitted.matrix)) {
      stop("The reconstructed random-effects basis differs from the fitted model; refit with 'mjoint()' before predicting.")
    }
  }
  for (i in seq_len(object$dims$n)) {
    times <- object$dmats$t$tj[seq_len(max(object$dmats$t$survdat2$tj.ind[i], 1L))]
    matrices <- lapply(seq_len(object$dims$K), function(k) {
      data <- setNames(data.frame(times), object$timeVar[k])
      .mjoint_new_design(designs$random[[k]], data)$matrix
    })
    matrix <- as.matrix(Matrix::bdiag(matrices))
    if (!same.values(matrix, object$dmats$z$Zi.fail[[i]])) {
      stop("This model was fitted with inconsistent random-effects bases at measurement and event times; refit with 'mjoint()' before predicting.")
    }
  }
  frame <- model.frame(object$sfit$terms, object$survData,
                       xlev = object$sfit$xlevels, na.action = na.fail)
  designs$event <- .mjoint_design(frame, object$sfit$x, drop.intercept = TRUE)
  if (!same.values(.mjoint_new_design(designs$event, object$survData)$matrix,
                   object$sfit$x)) {
    stop("The reconstructed event design differs from the fitted model; refit with 'mjoint()' before predicting.")
  }
  designs
}
