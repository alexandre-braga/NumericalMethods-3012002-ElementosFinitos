clc
clear all
close all

format long;
%dominio
a = 0.00;
b = 1.00;

delta1 = -1/2;
delta2 = 1/2;

erro = zeros(4,1);
hh = zeros(4,1);

for grau = 1:4
  for cont = 1:4
    nel = 4^cont;
    h = (b-a)/nel;
    k = grau;
    nen = k+1;
    np =  k*nel+1;
    nint = nen;
    
    A = zeros(np,np);
    B = zeros(np,np);
    C = zeros(np,np);
    
    F = zeros(np,1);

    %montagem do xl
    xl = zeros(np,1);
    xl(1) = a;
    for i = 2:np
      xl(i) = xl(i-1) + h/k;
    endfor

    %gera shg e pega as funções peso
    [shg, w]= shgGera(nen,nint);

    %Montagem da Fonte e da Matriz
    for n = 1:nel
      Ge = zeros(np,1);
      Fe = zeros(np,1);
      Ae = zeros(np,np);
      Be = zeros(np,np);
      Ce = zeros(np,np);
      for l = 1:nint
        xx = 0.;
        for i = 1:nen
          xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
        endfor
        for j = 1:nen
          Ge(k*(n-1)+j) = Ge(k*(n-1)+j) + delta2*funcao(xx)*shg(2,j,l)*2/h*w(l)*h/2;
          Fe(k*(n-1)+j) = Fe(k*(n-1)+j) + shg(1,j,l)*w(l)*h/2;
          for i = 1:nen
            Ae((k*(n-1)+i),(k*(n-1)+j)) = Ae((k*(n-1)+i),(k*(n-1)+j)) + ((1+delta1)*shg(1,i,l)*shg(1,j,l) + delta2*shg(2,i,l)*2/h*shg(2,j,l)*2/h)*w(l)*h/2;
            Be((k*(n-1)+i),(k*(n-1)+j)) = Be((k*(n-1)+i),(k*(n-1)+j)) + (-shg(1,i,l)*shg(2,j,l)*2/h + delta1*shg(2,i,l)*2/h*shg(2,j,l)*2/h)*w(l)*h/2;
            Ce((k*(n-1)+j),(k*(n-1)+i)) = Ce((k*(n-1)+j),(k*(n-1)+i)) + delta1*shg(2,i,l)*2/h*shg(2,j,l)*2/h)*w(l)*h/2;
          endfor
        endfor
      endfor
      F(n) = Fe(n);
      F(n + nel) = Ge(n);
      A = A + Ae;
      B = B + Be;
      C = C + Ce;
    endfor
    BT = B';
        
    %Matriz global
    M = zeros(np,np);
    
    %Caso o tamanho difira, usar outros loops
    %E modificar o tamanho da M, o mesmo pra F
    for i = 1:np
      for j = 1:np
        M(i,j) = A(i,j);
        M(i,j + nel) = B(i,j);
        M(i + nel,j) = BT(i,j);
        M(i + nel,j + nel) = C(i,j);
      endfor
    endfor
    
    %função exata
    x = a;
    t = T0;
    j = 1;
    espacoT = ceil((T-T0)/deltaT + 1)
    exata = zeros(np,espacoT);
    while t <= T
      for i = 1:np
        exata(i,j) = funcaoExata(x,t,E,Kappa);
        x += h/k;
      endfor
      x = a;
      j++;
      t += deltaT;
    endwhile
    x = a:h/k:b;
    t = T0:deltaT:T;
 
    %salva a resolucao
    nome = sprintf("log/PesosEPontosIntegracaoPeh%dGrau%d.txt", Peh, grau);
    save(nome, 'xl', 'h', 'deltaT', 'U', 'x', 't', 'exata');
      
  endfor 
  
endfor 
