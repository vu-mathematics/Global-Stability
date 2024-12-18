function [alphaopt,muopt,aaopt,mopt,ropt] = optimise_ppGpp(s);

for i=1:length(s)
  i
  PPGPP = @(alpha) -ppGpp_alpha(s(i),alpha,0);
  alphaopt(i) = fminsearch(PPGPP,0.5);
  [muopt(i),aaopt(i),mopt(i),ropt(i)] = ppGpp_alpha(s(i),alphaopt(i),0);
end  