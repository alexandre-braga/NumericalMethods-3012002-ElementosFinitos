clc
clear all
close all

format long;
%dominio
a = 0.00;
b = 1.00;

grau = 1
beta = 0

erro = zeros(4,1);
hh = zeros(4,1);

for alfa = -1:1
  for cont = 1:4
    nel = 4^cont;
    h = (b-a)/nel;
    k = grau;
    nen = k+1;
    np =  k*nel+1
    nint = nen;
    
    A = zeros(np,np);

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
      for l = 1:nint
        xx = 0.;
        for j = 1:nen
          xx = xx + shg(1,j,l)*xl(k*(n-1)+j);
        endfor
        for i = 1:nen
          F(k*(n-1)+i) = F(k*(n-1)+i) + funcao(xx)*shg(1,i,l)*w(l)*h/2;
          for j = 1:nen
            A((k*(n-1)+i),(k*(n-1)+j)) = A((k*(n-1)+i),(k*(n-1)+j)) + shg(2,i,l)*2/h*shg(2,j,l)*2/h**w(l)*h/2;
          endfor
        endfor
      endfor
    endfor
    
    %Contorno Nitsche
    A(1,1)     += -shg(2,1,1)*-1  + alfa*shg(2,1,1)*-1  + beta*funcaoExata(a)*-1;
    A(nel,nel) += -shg(2,nen,2)*1 + alfa*shg(2,nen,2)*1 + beta*funcaoExata(b)*1;
    
    F(1)   += alfa*shg(2,1,1)*funcaoExata(a) + beta*funcaoExata(a)*-1;
    F(nel) += alfa*shg(2,nen,2)*funcaoExata(b) + beta*funcaoExata(b)*1;
    
    u = A\F;
    
    %função exata
    x = a;
    exata = zeros(np,1);
    for i = 1:np
      exata(i) = funcaoExata(x);
      x += h/k;
    endfor
    x = a:h/k:b;
    
    %cálculo do erro L2
    erul2 = 0;
    for n = 1:nel
      eru = 0;
      for l = 1:nint
        uh = 0;
        xx = 0;
        for i = 1:nen
          uh = uh + shg(1,i,l)*u(k*(n-1)+i);
          xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
        endfor
        eru = eru + ((funcaoExata(xx) - uh)**2) * w(l) * h/2;
      endfor
      erul2 = erul2 + eru;
    endfor
    erul2 = sqrt(erul2);
    erro(cont) = erul2;
    hh(cont) = h;
    
    %salva a resolucao
    nome = sprintf("log/PesosEPontosIntegracao%dalfa%d.txt", cont, alfa);
    save(nome, 'alfa', 'beta', 'h', 'xl', 'u', 'x', 'exata');
    
  endfor 
  
  %salva os erros
  nome = sprintf("log/Erros%d.txt", alfa);
  save(nome, 'erro', 'hh' );
  
  alfa++;
  
endfor 
