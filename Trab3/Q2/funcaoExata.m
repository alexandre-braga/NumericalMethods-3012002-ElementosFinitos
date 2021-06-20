function [y] = funcaoExata(x)
    y = -1.90476*exp(-5*x) + 3.9604*exp(-x) + 5.99498*(1.e-44)*exp(100*x) - 2.05563;
endfunction
