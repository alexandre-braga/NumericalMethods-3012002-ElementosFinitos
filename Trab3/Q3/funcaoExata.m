function [y] = funcaoExata(x,t,E,K)
    y = 1/(sqrt(4*t + 1)) * exp(- ((x-K*t-0.5)^2) / (E*(4*t+1)) );
endfunction
