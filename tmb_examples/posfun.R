## Crash on 32 bit ?
library(TMB)
compile("posfun.cpp","-O0 -g")
dyn.load(dynlib("posfun"))
x <- rnorm(10)
mod = MakeADFun(data=list(x=x, eps=1e-3), parameters=list(p=0,Dummy=0), random=NULL)
o <- nlminb(mod$par, mod$fn, mod$gr)
##but doesn't work when "p" (the parameter governing the variable that enters posfun) is random:
mod = MakeADFun(data=list(x=x, eps=1e-3), parameters=list(p=0,Dummy=0), random="p")
o <- nlminb(mod$par, mod$fn, mod$gr)
