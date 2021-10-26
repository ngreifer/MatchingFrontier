getCutpoint <- function(dataset, base.form, cov, cutpoint.method){

    if (cutpoint.method == "median") {
        cutpoint <- median(dataset[[cov]])
    }
    else if (cutpoint.method == "mean") {
        cutpoint <- mean(dataset[[cov]])
    }
    else if (cutpoint.method == "segmented") {

        if (!requireNamespace("segmented", quietly = TRUE)) {
            customStop("package 'segmented' is required for this function to work. Please install it or use a different cutpoint method.")
        }

        mod.form <- update(base.form, as.formula(paste(". ~", cov)))

        if (length(unique(model.response(model.frame(base.form, dataset)))) == 2) {
            base.mod <- glm(mod.form, data = dataset, family = 'binomial')
        }
        else{
            base.mod <- lm(mod.form, data = dataset)
        }

        seg.reg <- segmented::segmented(base.mod, control = segmented::seg.control(it.max = 10000))
        cutpoint <- seg.reg$psi[2]
    }

    return(cutpoint)
}
