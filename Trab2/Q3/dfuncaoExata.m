function [y] = dfuncaoExata(x,E)
  %constantes da solucão exata
  c2 = (exp(-1/sqrt(E)) - 1) / ( exp(1/sqrt(E)) - exp(-1/sqrt(E)) );
  c1 = -1 - c2;
  y = ( c2*exp(x/sqrt(E)) - c1*exp(-x/sqrt(E)) ) / sqrt(E) ;
endfunction
