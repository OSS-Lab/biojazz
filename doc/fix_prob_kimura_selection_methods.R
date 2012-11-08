#!/usr/bin/R
# this is a script to generation plots of kimura selections

# the fixing probability is calculated by f = (1-exp(-2*s))/(1-exp(-4*Neff*s))
# where s = (w' - w) / w

fixation <- function(alpha, Neff, s) {
  poss = s[s>0];
  negs = s[s<0];
  fpos = (1 - exp(-2 * poss)) / (1 - exp(-4 * Neff * poss));
  fneg = (1 - exp(-2 * negs)) / (1 - exp(-4 * Neff * negs));
  if (0 %in% s) {
    zeros = s[s==0];
    fzero = array(1 / (2 * Neff), dim = dim(zeros));
    f = c(fneg, fzero, fpos);
  } else {
    f = c(fneg, fpos);    
  }
    
  p = alpha * f;
  return(p);
}

e = c(0.001:10);
pos = log(e);
Neff = 1e3;
alpha = 10;

p = fixation(alpha, Neff, s);

plot(s, p, type = "l", log = "y")#, xlim = c(-1e-4, 1e-5), ylim = c(1e-50, 100));