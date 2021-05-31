function [y] = funcaoExata(x,E)
  %constantes da solucão exata
  c2 = (exp(-1/sqrt(E)) - 1) / ( exp(1/sqrt(E)) - exp(-1/sqrt(E)) );
  c1 = -1 - c2;
  y = c1*exp(-x/sqrt(E)) + c2*exp(x/sqrt(E)) + 1;
endfunction
