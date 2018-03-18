
##----------------------------------------------------------------------------
##                GENERAL FUNCTIONS
##----------------------------------------------------------------------------

make.global <- function(alist, dest.env='.GlobalEnv') {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest.env );
    }
}

expit <- function(x) {
    ex <- exp(x);
    ex/(1+ex);
}

get.covmat <- function(StDevCovar, corrCovar) {
    n.x      <- length(StDevCovar);
    Vars     <- StDevCovar*StDevCovar;
    CovarMat <- matrix(NA, n.x, n.x);
    for (i in 1:n.x) {
        CovarMat[i,i] <- Vars[i];
        for (j in i:n.x) {
            if (j == i) {
                CovarMat[i,i] <- Vars[i];
                next;
            }
            CovarMat[i, j] <- corrCovar*StDevCovar[i]*StDevCovar[j];
            CovarMat[j, i] <- CovarMat[i, j];
        }
    }

    CovarMat
}

get.ps <- function(dta, ps.fml = NULL, ps.cov = NULL, grp="group", delta = 0,
                   type = c("randomforest", "logistic"), ntree = 5000, ...) {

    type <- match.arg(type);

    ## generate formula
    if (is.null(ps.fml))
        ps.fml <- as.formula(paste(grp, "~", paste(ps.cov, collapse="+"), sep=""));


    ## fit model
    switch(type,
           logistic = {glm.fit <- glm(ps.fml, family=binomial, data=dta, ...);
                       est.ps <- glm.fit$fitted;},
           randomforest = {dta[[grp]] <- as.factor(dta[[grp]]);
                           rf.fit     <- randomForest(ps.fml, data = dta, ntree = ntree, ...);
                           est.ps     <- predict(rf.fit, type = "prob")[,2];
                          });

    ##exponential tilting for z=1 only
    if (0 != delta) {
        e.delta        <- exp(delta);
        est.ps.delta   <- sapply(est.ps, function(x) {x*e.delta/(x*e.delta+1-x)});
        inx.z1         <- which(1 == dta[, grp]);
        est.ps[inx.z1] <- est.ps.delta[inx.z1];
    }

    est.ps
}

get.xbeta <- function(covX, regCoeff) {
    apply(covX, 1, function(x) {sum(x * regCoeff)});
}

get.rwe.class <- function(c.str = c("DWITHPS", "PSKL")) {
    c.str <- match.arg(c.str);
    switch(c.str,
           DWITHPS = "RWE_DWITHPS",
           PSKL    = "RWE_PSKL");
}
