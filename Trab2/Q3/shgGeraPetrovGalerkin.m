function [shg, w] = shgGeraPetrovGalerkin(nen,nint,E,h)
  %define pontos de integração e pesos
  switch (nint)
    case 2
      pt(1) = -sqrt(3)/3;
      pt(2) = sqrt(3)/3;
      w(1) = 1;
      w(2) = 1;
    %grau 2, 3 pts de integração
    case 3
      pt(1) = -sqrt(3/5);
      pt(2) = 0.0;
      pt(3) = sqrt(3/5);
      w(1) = 5/9;
      w(2) = 8/9;
      w(3) = 5/9;
    %grau 3, 4 pts de integração
    case 4
      pt(1) = -sqrt((3 + 2*sqrt(6/5))/7);
      pt(2) = -sqrt((3 - 2*sqrt(6/5))/7);
      pt(3) = sqrt((3 - 2*sqrt(6/5))/7);
      pt(4) = sqrt((3 + 2*sqrt(6/5))/7);
      w(1) = (18 - sqrt(30))/36;
      w(2) = (18 + sqrt(30))/36;
      w(3) = (18 + sqrt(30))/36;
      w(4) = (18 - sqrt(30))/36;
    case 5
      pt(1) = -1/3*sqrt(5 + 2*sqrt(10/7));
      pt(2) = -1/3*sqrt(5 - 2*sqrt(10/7));
      pt(3) = 0.0;
      pt(4) = 1/3*sqrt(5 - 2*sqrt(10/7));
      pt(5) = 1/3*sqrt(5 + 2*sqrt(10/7));
      w(1) = (322 - 13*sqrt(70))/900;
      w(2) = (322 + 13*sqrt(70))/900;
      w(3) = 128/225;
      w(4) = (322 + 13*sqrt(70))/900;
      w(5) = (322 - 13*sqrt(70))/900;
  endswitch
  %monta as funções de base
  for l = 1:nint
    t = pt(l);
    shg(1,1,l) = -( exp((t - 1)/sqrt(E)) - exp((1 - t)/sqrt(E)) )/( exp(h/sqrt(E)) - exp(-h/sqrt(E)) );
    shg(1,2,l) = ( exp((t - 1)/sqrt(E)) - exp((1 - t)/sqrt(E)) )/( exp(h/sqrt(E)) - exp(-h/sqrt(E)) );
    shg(2,1,l) = -( (1*sqrt(E))*exp((t - 1)/1*sqrt(E)) - (1*sqrt(E))*exp((1 - t)/1*sqrt(E)) ) / ( exp(h/sqrt(E)) - exp(-h/sqrt(E)) );
    shg(2,2,l) = ( (1*sqrt(E))*exp((t - 1)/1*sqrt(E)) - (1*sqrt(E))*exp((1 - t)/1*sqrt(E)) ) / ( exp(h/sqrt(E)) - exp(-h/sqrt(E)) );
  endfor 
endfunction
