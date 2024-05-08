###     -*- Coding: utf-8 -*-          ###
### Analyste: Charles-Édouard Giguère  ###
###                              .~    ###
###  _\\\\\_                    ~.~    ###
### |  ~ ~  |                 .~~.     ###
### #--O-O--#          ==||  ~~.||     ###
### |   L   |        //  ||_____||     ###
### |  \_/  |        \\  ||     ||     ###
###  \_____/           ==\\_____//     ###
##########################################

require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
theme_set(theme_bw() + theme(legend.position = "bottom"))
require(CUFF, quietly = TRUE, warn.conflicts = FALSE)


### Modèle distal
distal <- new.env()
load("./Distal_outcome.Rdata", envir = distal)

distal$test.par
(ft1 <- round(ftable(xtabs(-resMean_bias ~ sep_level + n + method, distal$test.par)), 2))
(ft2 <- round(ftable(xtabs(rmse ~ sep_level + n + method, distal$test.par)), 2))

xtable::xtableFtable(ft1)
xtable::xtableFtable(ft2)

### Modèle covariate.
covar <- new.env()
load("./Covariate.Rdata", envir = covar)

covar$test.par$resMean_bias <- apply((covar$test - 1), 1, mean)
covar$test.par$rmse
covar$test.par$method <- sprintf("%d (%s)",
                                 covar$test.par$n_steps,
                                 covar$test.par$correction) %>%
  factor(levels = c("1 (NA)", "2 (NA)", "3 (NA)", "3 (BCH)", "3 (ML)"))


(ft3 <- round(ftable(xtabs(-resMean_bias ~ sep_level + n + method, covar$test.par)), 2))
(ft4 <- round(ftable(xtabs(rmse ~ sep_level + n + method, covar$test.par)), 2))

xtable::xtableFtable(ft3)
xtable::xtableFtable(ft4)


### Modèle Complete.
complete <- new.env()
load("./complete.Rdata", envir = complete)
complete$test.par$method <- sprintf("%d (%s)",
                                    complete$test.par$n_steps,
                                    complete$test.par$correction) %>%
  factor(levels = c("1 (NA)",
                    "2 (NA)",
                    "3 (NA)",
                    "3 (BCH)",
                    "3 (ML)"))



(ft5 <- round(ftable(xtabs(-Dist_resMean_bias ~ nan_ratio + n + method,
                           complete$test.par)), 2))

(ft6 <- round(ftable(xtabs(Dist_rmse ~ nan_ratio + n + method,
                           complete$test.par)), 2))


xtable::xtableFtable(ft5)
xtable::xtableFtable(ft6)



