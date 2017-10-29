myggsurv <- function (object, ..., alpha = 0.5, shape = 3, xlab = "Time", 
          ylab = "Survival", title = "Marks show times with censoring", 
          legendLabs = NULL, CI = FALSE, bands = FALSE, pval = FALSE, 
          plotTable = FALSE, divideTime = 1, returnTable = FALSE) 
{
  stopifnot(inherits(object, "survfit"))
  if (!is.null(legendLabs)) 
    stopifnot(length(legendLabs) == length(object$strata))
  st1 <- stNames <- surv <- n.risk <- n.censor <- upper <- lower <- NULL
  if (is.null(legendLabs)) {
    stNames <- names(object$strata)
  }
  else {
    stNames <- legendLabs
  }
  st1 <- unlist(sapply(1:length(object$strata), function(i) rep(stNames[i], object$strata[i])))
  if (is.null(object$strata)) 
    st1 <- as.factor(rep(1, length(object$time)))
  dt1 <- data.table(time = object$time, n.risk = object$n.risk, 
                    n.event = object$n.event, n.censor = object$n.censor, 
                    surv = object$surv, upper = object$upper, lower = object$lower, 
                    strata = factor(st1))
  for (i in 1:length(unique(st1))) {
    dt2 <- dt1[dt1$strata == unique(st1)[i]]
    dt3 <- data.table(time = c(0, dt2$time[1]),
                      n.risk = rep(dt2$n.risk[1], 2),
                      n.event = c(0, 0),
                      n.censor = c(0, 0),
                      surv = c(1, 1),
                      upper = c(1, 1), lower = c(1, 1),
                      rep(dt2$strata[i], 2))
    dt1 <- rbindlist(list(dt1, dt3))
  }
  dt1 <- dt1[order(dt1$strata, (1 - dt1$surv)), ]
  if (!CI && !bands) {
    g1 <- ggplot(data = dt1, aes(group = strata, colour = strata)) + 
      geom_step(aes(x = time, y = surv), direction = "hv") + 
      geom_point(data = subset(dt1, n.censor >= 1),
                 aes(x = time, y = surv), shape = shape) + 
      scale_colour_brewer(type = "qual", palette = "Dark2", guide = guide_legend()) + 
      scale_x_continuous(xlab) + 
      scale_y_continuous(ylab) + 
      ggtitle(title)
  }
  if (CI && !bands) {
    g1 <- ggplot(data = dt1, aes(colour = strata, group = strata)) + 
      geom_step(aes(x = time, y = surv), direction = "hv", 
                fill = strata) + geom_step(aes(x = time, y = upper), 
                                           direction = "hv", linetype = 10, alpha = alpha) + 
      geom_step(aes(x = time, y = lower), direction = "hv", 
                linetype = 10, alpha = alpha) + geom_point(data = subset(dt1, 
                                                                         n.censor >= 1), aes(x = time, y = surv), shape = shape) + 
      scale_colour_brewer(type = "qual", palette = "Dark2", 
                          guide = guide_legend(keywidth = 3, keyheight = 3)) + 
      scale_x_continuous(xlab) + scale_y_continuous(ylab) + 
      ggtitle(title) + theme(legend.text = element_text(size = 15), 
                             legend.title = element_text(size = 15))
  }
  if (bands) {
    g1 <- ggplot(data = dt1, aes(colour = strata, group = strata)) + 
      geom_step(aes(x = time, y = surv), direction = "hv") + 
      geom_ribbon(aes(x = time, ymax = upper, ymin = lower, 
                      fill = strata), alpha = alpha) + geom_point(data = subset(dt1, 
                                                                                n.censor >= 1), aes(x = time, y = surv), shape = shape) + 
      scale_fill_brewer(type = "qual", palette = "Set2") + 
      scale_colour_brewer(type = "qual", palette = "Dark2", 
                          guide = guide_legend(keywidth = 3, keyheight = 3)) + 
      scale_x_continuous(xlab) + scale_y_continuous(ylab) + 
      ggtitle(title) + theme(legend.text = element_text(size = 15), 
                             legend.title = element_text(size = 15))
  }
  if (pval) {
    sd1 <- survival::survdiff(eval(object$call$formula), 
                              data = eval(object$call$data))
    p1 <- stats::pchisq(sd1$chisq, length(sd1$n) - 1, lower.tail = FALSE)
    p1txt <- ifelse(p1 < 1e-04, "p < 0.0001", paste("Log-rank test \n p =", 
                                                    signif(p1, 3)))
    g1 <- g1 + annotate("text", x = 0.1 * max(dt1$time), 
                        y = 0.2, label = p1txt)
  }
  times1 <- seq(0, max(object$time), by = divideTime)
  if (is.null(object$strata)) {
    df1 <- data.frame(strata = as.factor(rep(1, length(times1))), 
                      time = summary(object, times = times1, extend = TRUE)$time, 
                      n.risk = summary(object, times = times1, extend = TRUE)$n.risk)
  }
  else {
    df1 <- data.frame(strata = summary(object, times = times1, 
                                       extend = TRUE)$strata, time = summary(object, times = times1, 
                                                                             extend = TRUE)$time, n.risk = summary(object, times = times1, 
                                                                                                                   extend = TRUE)$n.risk)
    if (!is.null(legendLabs)) 
      df1$strata <- factor(df1$strata, labels = legendLabs)
  }
  if (plotTable) {
    tg1 <- ggplot(df1, aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) + 
      geom_text(size = 3.5) +
      scale_y_discrete(breaks = as.character(levels(df1$strata)), labels = levels(df1$strata)) + 
      scale_x_continuous(limits = c(0, max(object$time)), breaks = times1) + ggtitle("Number at risk by time") + 
      theme_grey() + 
      theme(plot.background = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_blank(), plot.title = element_text(size = rel(0.75)), 
            axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  if (!plotTable) 
    print(g1)
  if (plotTable) {
    g2 <- ggplot_gtable(ggplot_build(g1))
    leg <- which(sapply(g2$grobs, function(x) x$name) == 
                   "guide-box")
    legend <- g2$grobs[[leg]]
    grid.arrange(arrangeGrob(g1 + theme(legend.position = "none"), 
                             legend, tg1 + theme(legend.position = "none"), ncol = 2, 
                             widths = c(0.85, 0.15), heights = c(2, 0.75)))
  }
  # if (returnTable) 
  #  return(df1)
  return(g1)
}
