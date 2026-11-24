# Independent Gaussian and numerical references for conditional simulation.
conditional_case <- function(m, V, A, h) {
  L <- t(chol(V))
  list(m = m, L = L, B = A %*% L,
       logh = as.vector(log(h) + A %*% m))
}

conditional_draws <- function(base, n) {
  envelope <- joineRML:::.mjoint_rejection_envelope(base)
  vapply(seq_len(n), function(i) {
    joineRML:::.mjoint_rejection_draw(envelope, 10000)$b
  }, numeric(length(base$m)))
}

conditional_data <- function() {
  list(r = 1, q = 0, K = 1, nk = 3, y = c(0.5, 1.3, 0.9),
       X = cbind(1, 0:2), Z = matrix(1, 3, 1), v = NULL,
       Z.fail = matrix(1, 2, 1), IW.fail = diag(2), tj.ind = 2, tobs = 2)
}

conditional_theta <- function() {
  list(D = matrix(0.8), beta = c(0.2, 0.1), sigma2 = 0.4,
       gamma = 0.7, haz = c(0.05, 0.12))
}

test_that("longitudinal conditioning preserves correlated random effects", {
  d <- conditional_data()
  d$K <- 2
  d$r <- c(1, 2)
  d$nk <- c(2, 3)
  d$X <- as.matrix(Matrix::bdiag(matrix(1, 2, 1), cbind(1, 0:2)))
  d$Z <- d$X
  d$y <- c(1, 2, -1, 0, 0.5)
  d$tj.ind <- 0
  theta <- list(D = matrix(c(1,.2,.1,.2,2,.5,.1,.5,1), 3),
                beta = c(.1,.2,.3), sigma2 = c(.4,.7))
  base <- joineRML:::.mjoint_conditional_base(theta, d)
  Rinv <- diag(c(2.5, 2.5, rep(1/.7, 3)))
  V <- solve(solve(theta$D) + t(d$Z) %*% Rinv %*% d$Z)
  m <- V %*% t(d$Z) %*% Rinv %*% (d$y - d$X %*% theta$beta)
  expect_equal(base$m, as.vector(m), tolerance = 1e-12)
  expect_equal(tcrossprod(base$L), V, tolerance = 1e-12)
  set.seed(72)
  b <- t(conditional_draws(base, 5000))
  expect_true(all(abs(colMeans(b) - m) < 6 * sqrt(diag(V)/nrow(b))))
  se <- sqrt((V^2 + outer(diag(V), diag(V))) / (nrow(b)-1))
  expect_true(all(abs(cov(b) - V) < 6 * se))
})

test_that("the conditional density agrees with the joint model", {
  d <- conditional_data()
  theta <- conditional_theta()
  base <- joineRML:::.mjoint_conditional_base(theta, d)
  z <- c(-2, -.5, 0, 1, 3)
  b <- as.vector(base$m + base$L[1,1] * z)
  original <- vapply(b, joineRML:::logpb, numeric(1), theta = theta, data = d)
  transformed <- -z^2 / 2 - vapply(z, function(x) {
    sum(exp(base$logh + base$B * x))
  }, numeric(1))
  expect_equal(original - original[1], transformed - transformed[1], tolerance = 1e-12)
  for (j in c(0, 1, 2)) {
    d$tj.ind <- j
    d$Z.fail <- matrix(1, max(1,j), 1)
    d$IW.fail <- diag(max(1,j))
    theta$gamma <- 0
    base <- joineRML:::.mjoint_conditional_base(theta, d)
    expect_equal(nrow(base$B), 0L)
    expect_equal(base$m, (sum(d$y - d$X %*% theta$beta)/.4)/(1/.8+3/.4))
    draw <- joineRML:::.mjoint_rejection_draw(joineRML:::.mjoint_rejection_envelope(base), 1)
    expect_equal(draw$proposals, 1L)
  }
})

test_that("survival-conditioned moments agree with numerical integration", {
  cases <- list(
    moderate = list(m=.3, V=.18, a=.7, h=.2),
    late = list(m=.3, V=.18, a=.7, h=50),
    strong = list(m=.3, V=.18, a=5, h=2),
    opposing = list(m=0, V=1, a=c(-3,3), h=c(10,10)))
  for (name in names(cases)) {
    cc <- cases[[name]]
    target <- function(b) dnorm(b, cc$m, sqrt(cc$V), log=TRUE) - sum(cc$h * exp(cc$a*b))
    mode <- optimize(target, c(-10,10), maximum=TRUE)$maximum
    density <- function(b) vapply(b, function(x) exp(target(x)-target(mode)), numeric(1))
    # Split the domain so integration cannot miss a narrow peak whose
    # second-moment integrand is zero at its centre.
    integral <- function(f, upper=10) {
      breaks <- seq(-10,upper,length.out=max(2,ceiling((upper+10)*2)+1))
      sum(vapply(seq_len(length(breaks)-1L),function(i) {
        integrate(f,breaks[i],breaks[i+1L],rel.tol=1e-10,abs.tol=1e-12)$value
      },numeric(1)))
    }
    mass <- integral(density)
    moments <- vapply(1:2, function(k) {
      integral(function(b) b^k*density(b))/mass
    }, numeric(1))
    base <- conditional_case(cc$m, matrix(cc$V), matrix(cc$a, ncol=1), cc$h)
    set.seed(2400 + match(name, names(cases)))
    draws <- conditional_draws(base, 5000)
    for (k in 1:2) {
      expect_lt(abs(mean(draws^k)-moments[k]), 6*sd(draws^k)/sqrt(length(draws)), label=name)
    }
    for (p in c(.025,.5,.975)) {
      q <- uniroot(function(x) integral(density, x)/mass-p, c(-9,9))$root
      expect_lt(abs(mean(draws <= q)-p), 6*sqrt(p*(1-p)/length(draws)), label=name)
    }
  }
})

test_that("a correlated two-dimensional conditional has the right moments", {
  m <- c(.2,-.3)
  V <- matrix(c(1,.75,.75,2),2)
  a <- c(.7,-1.2)
  h <- .8
  # Integrate the single survival-sensitive combination a'b. Gaussian
  # conditioning supplies the remaining direction independently.
  mu <- sum(a*m)
  variance <- as.numeric(t(a) %*% V %*% a)
  density <- function(x) dnorm(x, mu, sqrt(variance))*exp(-h*exp(x))
  mass <- integrate(density,-20,20)$value
  mean.x <- integrate(function(x) x*density(x),-20,20)$value/mass
  var.x <- integrate(function(x) x^2*density(x),-20,20)$value/mass-mean.x^2
  direction <- as.vector(V %*% a)/variance
  expected.mean <- m + direction*(mean.x-mu)
  expected.V <- V + tcrossprod(direction)*(var.x-variance)
  set.seed(502)
  draws <- t(conditional_draws(conditional_case(m,V,matrix(a,1),h), 8000))
  expect_true(all(abs(colMeans(draws)-expected.mean) < 6*sqrt(diag(expected.V)/nrow(draws))))
  centered <- sweep(draws,2,expected.mean)
  for (i in 1:2) for (j in 1:2) {
    products <- centered[,i]*centered[,j]
    expect_lt(abs(mean(products)-expected.V[i,j]), 6*sd(products)/sqrt(nrow(draws)))
  }
})

test_that("the tangent bound remains valid away from the conditional mode", {
  base <- conditional_case(.3,matrix(.18),matrix(.7),50)
  env <- joineRML:::.mjoint_rejection_envelope(base)
  env$anchor <- env$anchor+.5
  env$logw <- as.vector(base$logh + base$B %*% env$anchor)
  env$center <- -as.vector(crossprod(base$B,exp(env$logw)))
  for (z in seq(-8,8,length.out=31)) {
    log.target <- -z^2/2-sum(exp(base$logh+base$B*z))
    log.proposal <- dnorm(z,env$center,1,log=TRUE)
    log.accept <- joineRML:::.mjoint_rejection_log_accept(z,env)
    # The difference is a constant, independently derived at the anchor.
    constant <- -env$anchor^2/2-sum(exp(env$logw))-
      dnorm(env$anchor,env$center,1,log=TRUE)
    expect_equal(log.proposal+log.accept+constant, log.target, tolerance=1e-12)
    expect_lte(log.accept,0)
  }
})

test_that("acceptance calculations handle cancellation and overflow", {
  gap <- joineRML:::.mjoint_log_tangent_gap
  expect_equal(exp(gap(c(-1e-12,0,1e-12))), c(5e-25,0,5e-25), tolerance=1e-12)
  expect_equal(exp(-1000+gap(1000)),1)
  expect_equal(gap(c(-Inf,Inf)),c(Inf,Inf))
  expect_equal(exp(gap(-1000)),999)
  env <- list(B=matrix(1),anchor=0,logw=-1000)
  expect_equal(joineRML:::.mjoint_rejection_log_accept(1000,env),-1)
  expect_equal(joineRML:::.mjoint_rejection_log_accept(Inf,env),-Inf)
  expect_error(gap(NA_real_),"undefined")
})

test_that("invalid parameters and exhausted budgets stop the calculation", {
  d <- conditional_data()
  theta <- conditional_theta()
  for (bad in c(-1, Inf, NA_real_)) {
    invalid <- theta
    invalid$haz[1] <- bad
    expect_error(joineRML:::.mjoint_conditional_base(invalid,d),"finite and nonnegative")
  }
  invalid <- theta
  invalid$D <- matrix(-1)
  expect_error(joineRML:::.mjoint_conditional_base(invalid,d),"positive definite")
  invalid <- theta
  invalid$sigma2 <- 0
  expect_error(joineRML:::.mjoint_conditional_base(invalid,d),"strictly positive")
  base <- conditional_case(0,matrix(1),matrix(1),1)
  base$logh <- 1000
  expect_error(joineRML:::.mjoint_rejection_envelope(base),"finite.*envelope")
  env <- joineRML:::.mjoint_rejection_envelope(conditional_case(0,matrix(1),matrix(1),1))
  # A valid, deliberately poor tangent makes all these seeded proposals fail.
  env$anchor <- 10
  env$logw <- 10
  env$center <- -exp(10)
  set.seed(33)
  expect_error(joineRML:::.mjoint_rejection_draw(env,3),"exhausted max.proposals = 3")
  calls <- 0L
  draw <- joineRML:::.mjoint_prediction_draw
  local <- new.env(parent=environment(draw))
  local$thetaDraw <- function(object) { calls <<- calls+1L; theta }
  local$.mjoint_rejection_envelope <- function(base) env
  environment(draw) <- local
  expect_error(draw(NULL,d,3,2,10),"Simulation 2 of 10 at landmark 2.*no parameter draw was discarded")
  expect_equal(calls,1L)
})

test_that("conditional survival sums only the requested hazard increments", {
  d <- conditional_data()
  theta <- conditional_theta()
  theta$gamma <- 0
  theta$haz <- c(1000,1)
  expect_equal(joineRML:::S(0,theta,d),0)
  expect_equal(joineRML:::S(0,theta,d,start=1),exp(-1))
  theta$haz[1] <- 1e300
  expect_equal(joineRML:::S(0,theta,d,start=1),exp(-1))
  expect_equal(joineRML:::S(0,theta,d,start=2),1)
  theta$gamma <- 1
  theta$haz <- c(0,exp(-700))
  expect_equal(joineRML:::S(700,theta,d),exp(-1))
  expect_equal(joineRML:::S(1000,theta,d),0)
  expect_error(joineRML:::S(0,theta,d,start=3),"precedes")
})

# Supply a known parameter distribution while retaining both prediction
# functions' actual conditional sampler and summarisation calculations.
prediction_with_parameters <- function(prediction, theta.draw) {
  draw <- joineRML:::.mjoint_prediction_draw
  draw.env <- new.env(parent=environment(draw))
  draw.env$thetaDraw <- theta.draw
  environment(draw) <- draw.env
  prediction.env <- new.env(parent=environment(prediction))
  prediction.env$.mjoint_prediction_draw <- draw
  environment(prediction) <- prediction.env
  prediction
}

test_that("complete longitudinal predictions pair each parameter and random-effects draw", {
  fit <- fit_small_joint(formLongFixed=list(y1=y1~1))
  nd <- subset(fit$data[[1]],id==1)
  nd$y1 <- 0
  theta <- fit$coefficients
  theta$D <- matrix(1)
  theta$sigma2 <- .1
  theta$gamma[] <- 0
  calls <- 0L
  predict <- prediction_with_parameters(dynLong,function(object) {
    calls <<- calls+1L
    theta$beta <- rnorm(1)
    theta
  })
  set.seed(1203)
  n <- 5000
  result <- predict(fit,nd,u=3,type="simulated",M=n,progress=FALSE)
  # Three zero measurements give b | beta,y ~ N(-30 beta/31, 1/31).
  variance <- 1/31 + 1/31^2
  expect_lt(abs(result$pred[[1]]$mean),6*sqrt(variance/n))
  expect_lt(abs(result$pred[[1]]$se^2-variance),6*variance*sqrt(2/(n-1)))
  expect_equal(calls,n)
  expect_identical(result$sampler,"conditional-rejection")
  expect_equal(result$proposals,rep(1L,n))
  expect_equal(result$accept,1)
  expect_output(print(result),"Rejection-sampling acceptance rate")
  set.seed(1203)
  expect_identical(predict(fit,nd,u=3,type="simulated",M=n,progress=FALSE),result)
})

test_that("complete survival predictions agree with a parameter-mixture reference", {
  fit <- fit_small_joint(formLongFixed=list(y1=y1~1))
  nd <- subset(fit$data[[1]],id==1)
  nd$y1 <- 0
  tj <- fit$dmats$t$tj
  landmark <- median(tj)
  nd$time <- seq(0,landmark,length.out=3)
  j <- sum(tj<=landmark)
  theta <- fit$coefficients
  theta$D <- matrix(1)
  theta$sigma2 <- .1
  theta$gamma[] <- .8
  theta$haz <- c(rep(3/j,j),rep(1/(length(tj)-j),length(tj)-j))
  calls <- 0L
  predict <- prediction_with_parameters(dynSurv,function(object) {
    calls <<- calls+1L
    theta$beta <- c(-.7,.4)[(calls-1L)%%2L+1L]
    theta
  })
  # Conditional survival over the future interval is exp(-exp(.8*b)).
  moments <- vapply(c(-.7,.4),function(beta) {
    log.target <- function(b) dnorm(b,0,1,log=TRUE) +
      3*dnorm(0,beta+b,sqrt(.1),log=TRUE)-3*exp(.8*b)
    mode <- optimize(log.target,c(-5,5),maximum=TRUE)$maximum
    density <- function(b) exp(log.target(b)-log.target(mode))
    mass <- integrate(density,-8,8)$value
    vapply(1:2,function(k) integrate(function(b) {
      exp(-k*exp(.8*b))*density(b)
    },-8,8)$value/mass,numeric(1))
  },numeric(2))
  reference <- rowMeans(moments)
  set.seed(3302)
  n <- 5000
  result <- predict(fit,nd,u=c(landmark+.01,max(tj)),type="simulated",M=n,progress=FALSE)
  expect_lt(abs(tail(result$pred$mean,1)-reference[1]),
            6*sqrt((reference[2]-reference[1]^2)/n))
  expect_equal(calls,n)
  expect_true(all(diff(result$pred$mean)<=0))
  expect_true(all(result$pred$lower<=result$pred$median & result$pred$median<=result$pred$upper))
  expect_equal(result$accept,n/sum(result$proposals))
  expect_null(result$b.hat)
  expect_output(print(result),"Rejection-sampling acceptance rate")
})

test_that("both prediction functions validate limits and deprecate scale", {
  fit <- fit_small_joint()
  nd <- subset(fit$data[[1]],id==1)
  for (method in list(dynLong,dynSurv)) {
    predict <- prediction_with_parameters(method,function(object) object$coefficients)
    for (bad in list(0,-1,NA_real_,Inf,1.5,c(1,2),"2")) {
      expect_error(predict(fit,nd,u=3,type="simulated",M=bad,progress=FALSE),
                   "'M' must be a finite, positive integer")
      expect_error(predict(fit,nd,u=3,type="simulated",M=2,max.proposals=bad,progress=FALSE),
                   "'max.proposals' must be a finite, positive integer")
    }
    for (bad in list(NA_real_,Inf,c(.5,.9),numeric(),".95")) {
      expect_error(predict(fit,nd,u=3,type="simulated",M=2,ci=bad,progress=FALSE),"ci argument")
    }
    set.seed(921)
    result <- predict(fit,nd,u=3,type="simulated",M=10,progress=FALSE)
    set.seed(921)
    expect_warning(scaled <- suppressMessages(predict(fit,nd,u=3,type="simulated",M=10,
                                                      progress=FALSE,scale=10)),"scale.*deprecated and ignored")
    expect_equal(scaled$pred,result$pred)
    expect_equal(scaled$proposals,result$proposals)
    # Previously saved prediction objects retain their original print label.
    result$sampler <- NULL
    expect_output(print(result),"M-H acceptance rate")
  }
})

test_that("multivariate survival terms retain outcome-specific associations", {
  d <- list(r=c(1,2),q=1,K=2,nk=c(2,3),y=c(1,2,-1,0,.5),v=1.2,
            X=as.matrix(Matrix::bdiag(matrix(1,2,1),cbind(1,0:2))),
            Z=as.matrix(Matrix::bdiag(matrix(1,2,1),cbind(1,0:2))),
            Z.fail=as.matrix(Matrix::bdiag(matrix(1,2,1),cbind(1,c(1,3)))),
            IW.fail=cbind(diag(2),diag(2)),tj.ind=2)
  theta <- list(D=matrix(c(1,.2,.1,.2,2,.5,.1,.5,1),3),
                beta=c(.1,.2,.3),sigma2=c(.4,.7),gamma=c(.2,.3,-.7),haz=c(.2,.4))
  base <- joineRML:::.mjoint_conditional_base(theta,d)
  reference <- transformed <- numeric(4)
  z <- rbind(c(0,0,0),c(1,-1,.5),c(-2,1,0),c(.2,.5,-1))
  for (i in seq_len(nrow(z))) {
    b <- as.vector(base$m+base$L %*% z[i,])
    log.hazard <- .2*1.2+.3*b[1]-.7*(b[2]+c(1,3)*b[3])
    reference[i] <- mvtnorm::dmvnorm(b,sigma=theta$D,log=TRUE) +
      sum(dnorm(d$y,as.vector(d$X %*% theta$beta+d$Z %*% b),
                sqrt(c(.4,.4,.7,.7,.7)),log=TRUE))-sum(theta$haz*exp(log.hazard))
    transformed[i] <- -sum(z[i,]^2)/2-sum(exp(base$logh+base$B %*% z[i,]))
    expect_equal(joineRML:::S(b,theta,d,start=1),exp(-.4*exp(log.hazard[2])))
  }
  expect_equal(reference-reference[1],transformed-transformed[1],tolerance=1e-12)
})

test_that("complete survival predictions remain defined after survival underflows", {
  fit <- fit_small_joint(formLongFixed=list(y1=y1~1))
  tj <- fit$dmats$t$tj
  landmark <- median(tj)
  nd <- subset(fit$data[[1]],id==1)
  nd$time <- seq(0,landmark,length.out=3)
  nd$y1 <- 0
  fit$coefficients$beta[] <- 0
  fit$coefficients$gamma[] <- 0
  fit$coefficients$haz[] <- 0
  fit$coefficients$haz[c(1,length(tj))] <- c(1000,1)
  for (type in c("first-order","simulated")) {
    predict <- prediction_with_parameters(dynSurv,function(object) object$coefficients)
    result <- predict(fit,nd,type=type,M=10,progress=FALSE)
    expect_equal(result$pred$u[1],landmark)
    values <- if (type=="first-order") result$pred$surv else result$pred$mean
    expect_equal(values[1],1)
    expect_equal(tail(values,1),exp(-1))
    expect_true(all(is.finite(as.matrix(result$pred))))
  }
})
