function [shg, w] = shgGeraPetrovGalerkin(nen,nint,E,h)
  %define pontos de integração e pesos
  switch (nint)
    %grau 1, 2 pts de integração
    case 2
      pt(1) = -sqrt(3)/3;
      pt(2) = sqrt(3)/3;
      w(1) = 1;
      w(2) = 1;
  endswitch
  %monta as funções de base
  for l = 1:nint
    t = pt(l);
    switch (nen)
      %grau 1
      case 2
        shg(1,1,l) = -( exp((t - 1)/sqrt(E)) - exp((1 - t)/sqrt(E)) )/( exp(h/sqrt(E)) - exp(-h/sqrt(E)) );
        shg(1,2,l) = ( exp((t - 1)/sqrt(E)) - exp((1 - t)/sqrt(E)) )/( exp(h/sqrt(E)) - exp(-h/sqrt(E)) );
        shg(2,1,l) = -( (1*sqrt(E))*exp((t - 1)/1*sqrt(E)) - (1*sqrt(E))*exp((1 - t)/1*sqrt(E)) ) / ( exp(h/sqrt(E)) - exp(-h/sqrt(E)) );
        shg(2,2,l) = ( (1*sqrt(E))*exp((t - 1)/1*sqrt(E)) - (1*sqrt(E))*exp((1 - t)/1*sqrt(E)) ) / ( exp(h/sqrt(E)) - exp(-h/sqrt(E)) );
    endswitch
  endfor 
endfunction
