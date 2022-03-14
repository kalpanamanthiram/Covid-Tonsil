plotModelRes <- function(mod.list, 
                         base_size = 11,
                         pvalue.text.size = 2, 
                         scales = "free", 
                         scale_y_discrete_expand = c(0.1, 0.1),
                         labeller = "label_value"){
  
  # coefs:      data.frame with coefficients for each covariate
  # error.bars: data.frame with error (e.g., Std. Error) for each covariate
  # pvals:      data.frame with p-values for each covariate
  
  coefs <- lapply(X = mod.list,
                  FUN = function(mod){
                    if(all(class(mod) == "lm") |
                       all(class(mod) == c("negbin", "glm", "lm")) |
                       all(class(mod) == "glmerMod")){
                      mod %>%
                        summary() %>%
                        .$coefficients %>%
                        .[-1,1]
                    } else if(all(class(mod) == "lme")){
                      mod %>%
                        summary() %>%
                        .$tTable %>%
                        .[-1,1]
                    }
                  }) %>%
    do.call(rbind,
            .) %>%
    as.data.frame()
  
  if(ncol(coefs) == 1){
    colnames(coefs) <- formula(mod.list[[1]]) %>%
      terms() %>%
      attr("term.labels")
  }
  
  ci.lo <- lapply(X = mod.list,
                  FUN = function(mod){
                    if(all(class(mod) == "lm") |
                       all(class(mod) == c("negbin", "glm", "lm"))){
                      mod %>%
                        confint(object = .,
                                level = 0.95) %>%
                        .[-1,1]  %>%
                        suppressMessages()
                    } else if(all(class(mod) == "lme")){
                      mod %>%
                        nlme::intervals(object = .,
                                        level = 0.95, 
                                        which = "fixed") %>%
                        .$fixed %>%
                        .[-1,1]
                    } else if(all(class(mod) == "glmerMod")){
                      (mod %>%
                        summary() %>%
                        .$coefficients %>%
                        .[-1,1]) - (mod %>%
                        summary() %>%
                        .$coefficients %>%
                        .[-1,2] %>%
                        "*"(2))
                    }
                  }) %>%
    do.call(rbind,
            .) %>%
    as.data.frame()
  
  ci.hi <- lapply(X = mod.list,
                  FUN = function(mod){
                    if(all(class(mod) == "lm") |
                       all(class(mod) == c("negbin", "glm", "lm"))){
                      mod %>%
                        confint(object = .,
                                level = 0.95) %>%
                        .[-1,2]   %>%
                        suppressMessages()
                    } else if(all(class(mod) == "lme")){
                      mod %>%
                        nlme::intervals(object = .,
                                        level = 0.95, 
                                        which = "fixed") %>%
                        .$fixed %>%
                        .[-1,3]
                    } else if(all(class(mod) == "glmerMod")){
                      (mod %>%
                         summary() %>%
                         .$coefficients %>%
                         .[-1,1]) + (mod %>%
                                       summary() %>%
                                       .$coefficients %>%
                                       .[-1,2] %>%
                                       "*"(2))
                    }
                  }) %>%
    do.call(rbind,
            .) %>%
    as.data.frame()
  
  p.vals <- lapply(X = mod.list,
                   FUN = function(mod){
                     if(all(class(mod) == "lm") |
                        all(class(mod) == c("negbin", "glm", "lm")) |
                        all(class(mod) == "glmerMod")){
                       mod %>%
                         summary() %>%
                         .$coefficients %>%
                         .[-1,4]  
                     } else if(all(class(mod) == "lme")){
                       mod %>%
                         summary() %>%
                         .$tTable %>%
                         .[-1,5]
                     }
                   }) %>%
    do.call(rbind,
            .) %>%
    as.data.frame()
  
  coefs$cluster <- rownames(coefs)
  coefs$cluster <- factor(x = coefs$cluster,
                          levels = rev(coefs$cluster))
  coefs.melt <- reshape2::melt(coefs,
                               id.vars = "cluster")
  colnames(coefs.melt)[3] <- "coefficient"
  
  
  ci.lo$cluster <- rownames(ci.lo)
  coefs.melt$confint.lo <- reshape2::melt(data = ci.lo,
                                     id.vars = "cluster")$value
  
  ci.hi$cluster <- rownames(ci.hi)
  coefs.melt$confint.hi <- reshape2::melt(data = ci.hi,
                                     id.vars = "cluster")$value
  
  p.vals$cluster <- rownames(p.vals)
  coefs.melt$p.value <- reshape2::melt(data = p.vals,
                                       id.vars = "cluster")$value %>%
    formatC(x = .,
            digits = 0,
            format = "e")
  
  position_nudge_y <- 
    ifelse(test = length(mod.list) == 1,
           yes = -0.1,
           no = -0.4)
  
  ggplot(data = coefs.melt,
         mapping = aes(x = coefficient,
                       y = cluster,
                       label = paste0("(",
                                      p.value,
                                      ")"))) +
    theme_bw(base_size = base_size) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank()) +
    xlab("estimate \u00B1 conf. int.\n(p-values)") +
    ylab("") +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               color = "grey50") +
    geom_point(mapping = aes(color = as.numeric(p.value) <= 0.05)) +
    geom_errorbar(mapping = aes(xmin = confint.lo,
                                xmax = confint.hi,
                                color = as.numeric(p.value) <= 0.05),
                  width = 0.2) +
    geom_text(mapping = aes(color = as.numeric(p.value) <= 0.05),
              size = I(pvalue.text.size),
              hjust = 0.5,
              position = position_nudge(y = position_nudge_y), 
              fontface = "bold") +
    scale_color_manual(values = c("black",
                                  "red")) + 
    #scale_x_continuous(expand = c(0.2, 0.2)) +
    scale_y_discrete(expand = scale_y_discrete_expand) +
    facet_grid(~ variable,
               scales = scales,
               labeller = labeller)
  
}
