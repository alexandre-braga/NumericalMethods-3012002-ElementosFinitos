clc
clear all
close all

format long;
%dominio
a = 0.00;
b = 1.00;

##  F(1) += + alfa*(shg(2,nen,2)*2/h*funcaoExata(xl(n+1)) - shg(2,1,1)*2/h*funcaoExata(xl(n))) + beta*(shg(1,nen,2)*funcaoExata(xl(n+1)) - shg(1,1,1)*funcaoExata(xl(n)));
##  F(np) += + alfa*(shg(2,nen,2)*2/h*funcaoExata(xl(n+1)) - shg(2,1,1)*2/h*funcaoExata(xl(n))) + beta*(shg(1,nen,2)*funcaoExata(xl(n+1)) - shg(1,1,1)*funcaoExata(xl(n)));
##  A(1,1) += -(shg(2,nen,2)*2/h*shg(1,1,2) - shg(2,nen,1)*2/h*shg(1,1,1)) + alfa*(shg(2,nen,2)*2/h*shg(1,1,2) - shg(2,nen,1)*2/h*shg(1,1,1)) + beta*(shg(1,nen,2)*shg(1,nen,2) - shg(1,1,1)*shg(1,1,1));
##  A(np,np) += -(shg(2,nen,2)*2/h*shg(1,1,2) - shg(2,nen,1)*2/h*shg(1,1,1)) + alfa*(shg(2,nen,2)*2/h*shg(1,1,2) - shg(2,nen,1)*2/h*shg(1,1,1)) + beta*(shg(1,nen,2)*shg(1,nen,2) - shg(1,1,1)*shg(1,1,1));
## 

grau = 1;
alfa = input('Insira o valor de alfa: ');
beta = input('Insira o valor de beta: ');

erro = zeros(4,1);
hh = zeros(4,1);

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
  [shge] = shgeGera(nen,nint);

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

  %Contorno Nitsche-  
  for i = 1:nen
    if i == 1
      F(i) += alfa*(shge(2,i,2)*2/h*funcaoExata(xl(np)) - shge(2,i,1)*2/h*funcaoExata(xl(1)) ) + beta*(funcaoExata(xl(np))*shg(1,i,2)-funcaoExata(xl(1))*shg(1,i,1));
    endif
    if i == nen
      F(np) += alfa*(shge(2,i,2)*2/h*funcaoExata(xl(np)) - shge(2,i,1)*2/h*funcaoExata(xl(1)) )  + beta*(funcaoExata(xl(np))*shg(1,i,2)-funcaoExata(xl(1))*shg(1,i,1));
    endif
    
    for j = 1:nen
      if (i == 1 && j == 1)
        A(i,j)  += -(shge(2,j,2)*2/h*shge(1,i,2) - shge(2,j,1)*2/h*shge(1,i,1)) + alfa*(shge(2,i,2)*2/h*shge(1,j,2) - shge(2,i,1)*2/h*shge(1,j,1)) + beta*(shg(1,j,2)*shg(1,i,2) - shg(1,j,1)*shg(1,i,1));
      endif
      if (i == nen && j == nen)
        A(np,np) += -(shge(2,j,2)*2/h*shge(1,i,2)  - shge(2,j,1)*2/h*shg(1,i,1)) + alfa*(shge(2,i,2)*2/h*shge(1,j,2) - shge(2,i,1)*2/h*shge(1,j,1)) + beta*(shg(1,j,2)*shg(1,i,2) - shg(1,j,1)*shg(1,i,1));
      endif
    endfor
  endfor

  u = A\F;
  
  %função exata
  x = a;
  exata = zeros(np,1);
  for i = 1:np
    exata(i) = funcaoExata(x);
    x += h/k;
  endfor
  x = a:h/k:b;
  
  %salva a resolucao
  nome = sprintf("log/PesosEPontosIntegracao%dalfa%dbeta%d.txt", cont, alfa, beta);
  save(nome, 'alfa', 'beta', 'h', 'xl', 'u', 'x', 'exata');
  
endfor 

