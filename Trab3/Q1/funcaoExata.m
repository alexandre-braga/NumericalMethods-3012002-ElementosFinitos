function [y] = funcaoExata(x,E,K)
    y = 1/K * (x - ( (1 - exp(K*x/E)) / (1 - exp(K/E)) ) );
endfunction
