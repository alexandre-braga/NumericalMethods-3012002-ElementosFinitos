function [y] = dfuncaoExata(x,E,K)
  y = 1/K * ( K*exp(K*x/E) / (E*(1-exp(K/E))) ) + 1 ;
endfunction
