function [shg, w] = shgGera(nen,nint)
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
        shg(1,1,l) = 0.5*(1.0 - t);
        shg(1,2,l) = 0.5*(1.0 + t);
        shg(2,1,l) = -1.0/2.0;
        shg(2,2,l) = 1.0/2.0;
    endswitch
  endfor 
endfunction
