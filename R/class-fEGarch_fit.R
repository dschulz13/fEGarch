
setClass("fEGarch_fit",
  slots = c(
    pars = "numeric",
    se = "numeric",
    vcov_mat = "matrix",
    rt = "ANY",
    cmeans = "ANY",
    sigt = "ANY",
    etat = "ANY",
    llhood = "numeric",
    inf_criteria = "numeric",
    cond_dist = "character",
    orders = "numeric",
    long_memo = "logical"
  ),
  prototype = list(
    pars = numeric(1),
    se = numeric(1),
    vcov_mat = matrix(numeric(1)),
    rt = numeric(1),
    cmeans = numeric(1),
    sigt = numeric(1),
    etat = numeric(1),
    llhood = numeric(1),
    inf_criteria = numeric(1),
    cond_dist = character(1),
    orders = numeric(1),
    long_memo = logical(1)
  )
)

setClass("fEGarch_fit_egarch",
  contains = "fEGarch_fit",
  slots = c(
    powers = "numeric",
    modulus = "logical"
  ),
  prototype = list(
    orders = numeric(1),
    modulus = logical(1)
  )
)

setClass("fEGarch_fit_loggarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_fiaparch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_lmgarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_rugarch_wrapper",
  contains = "fEGarch_fit",
  slots = c(
    rugarch_model = "ANY",
    model = "character"
  ),
  prototype = list(
    rugarch_model = numeric(1),
    model = character(1)
  )
)

fEGarch_fit_egarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    cmeans,
    sigt,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    powers,
    modulus,
    long_memo
) {

  new("fEGarch_fit_egarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    cmeans = cmeans,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    powers = powers,
    modulus = modulus,
    long_memo = long_memo
  )

}

fEGarch_fit_loggarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo
) {

  new("fEGarch_fit_loggarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo
  )

}


fEGarch_fit_fiaparch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo
) {

  new("fEGarch_fit_fiaparch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo
  )

}

fEGarch_fit_lmgarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo
) {

  new("fEGarch_fit_lmgarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo
  )

}

fEGarch_fit_rugarch_wrapper <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    rugarch_model,
    model
) {

  new("fEGarch_fit_rugarch_wrapper",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    rugarch_model = rugarch_model,
    model = model
  )

}





