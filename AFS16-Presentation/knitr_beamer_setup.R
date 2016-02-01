library(knitr)
## knitr custom hooks
# set default plot options
knit_hooks$set(par1 = function(before, options, envir) {
  if (before) par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)
})
## knitr options -- output look
opts_chunk$set(size = 'footnotesize',prompt=TRUE,comment='')
## knitr options -- figure handling
opts_chunk$set(fig.path='Figs/', fig.width=5, fig.height=5, out.width='.5\\linewidth', 
               fig.show='hold', fig.align='center', par1=TRUE)
## r options
options(show.signif.stars=FALSE,continue=" ")

