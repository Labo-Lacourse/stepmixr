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

### Example of probability.
p <- c(Group1 = 0.71, Group2 = 0.60, Group3 = 0.80)

### Function to reparameterized probability into logit.
p2logit <- function(p, ref=1){
    if(ref > length(p))
        stop("reference group is higher than the number of group")
    ## 1. Change p to odds.
    odds <- p/(1-p)
    ## 2. Set the reference group.
    odds_ref <- odds[ref]
    ## 3. set the coefficients and intercept.
    intercept <- log(odds_ref)
    coef <- log(odds/odds_ref)[-ref]
    res <- c(intercept, coef)
    names(res) <- c(sprintf("Int(grp=%d)", ref),
                    sprintf("Beta(grp=%d)", (1:length(p))[-ref]))
    res
}


p

### Change to logit using default 1st group as reference.
p2logit(p)

### Choosing last group instead as a reference.
p2logit(p, 3)

