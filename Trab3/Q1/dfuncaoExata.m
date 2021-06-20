function [y] = dfuncaoExata(x,E,K)
  y = ( E*(1-exp(K/E)) + K*(1-exp(K*x/E)) ) / E*K*(1-exp(K/E)) ;
endfunction
