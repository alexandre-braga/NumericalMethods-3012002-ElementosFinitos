clc
clear all
close all
%dominio
a = 0.0;
b = 1.5;
erro = zeros(4,1);
hh = zeros(4,1);

for cont = 2:5
  %n de elementos
  nel = 2^(cont);
  %tamanho do elemento
  h = (b-a)/nel
  %grau do polinomio de interpolação
  k = 1;
  %n de nós do elemento
  nen = k+1
  %n de nós global
  np =  k*nel+1;
  %n de pontos de integração
  nint = nen;
  %matriz global zerada
  K = zeros(np,np);
  %vetor fonte zerado
  F = zeros(np,1);

  %montagem do xl !!!LINEAR!!!
  xl(1) = a;
  for i = 2:np
    xl(i) = xl(i-1) + h/k;
  endfor

  %gera shg e pega as funções peso
  [shg, w]= shgGera(nen,nint)

  %montagem global
  for n = 1:nel
    %matriz do elemento %Ke = zeros(nen,nen); %vetor fonte do elemento %Fe = zeros(nen,1);
    for l = 1:nint
      xx = 0;
      for i = 1:nen
        %xl !!!LINEAR!!!
        xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
      endfor
      for j = 1:nen
        %construção direto global (aula 6)
        F(k*(n-1)+j) = F(k*(n-1)+j) + funcao(xx)*shg(1,j,l)*w(l)*h/2;
        for i = 1:nen
          K((k*(n-1)+i),(k*(n-1)+j)) = K((k*(n-1)+i),(k*(n-1)+j)) + shg(2,i,l)*2/h*shg(2,j,l)*2/h*w(l)*h/2;
        endfor
      endfor
    endfor
  endfor

  %Condição de Dirichlet
  K(nel+1,nel+1) += 10e4;
  K(nel+1,nel) = 0;
  F(nel+1) = 10e4*(-1);
  F(nel) = F(nel) + (1*K(nel,nel+1));
  K(nel,nel+1) = 0;

  %Condição de Neumann
  F(1) = F(1) - pi

  %função exata
  x = a:h:b;
  exata = funcaoExata(x);
  u = zeros(np);
  u = K\F;

  %cálculo do erro
    erdul2 = 0;
    for n = 1:nel
      erdu = 0;
      for l = 1:nint
        duh = 0;
        xx = 0;
        for i = 1:nen
          duh = duh + shg(2,i,l)*2/h*u(k*(n-1)+i);
          xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
        endfor
        erdu = erdu + ((dfuncaoExata(xx) - duh)**2) * w(l) * h/2;
      endfor
      erdul2 = erdul2 + erdu;
    endfor
    erdul2 = sqrt(erdul2);
    erro(cont) = erdul2;
    hh(cont) = h;
    
    %salva a resolucao
    nome = sprintf("log/PesosEPontosIntegrecao%d.txt", cont-1);
    save(nome, 'xl', 'h', 'u', 'x', 'exata');
    
endfor 

%salva os erros
nome = sprintf("log/Erros.txt");
save(nome, 'erro', 'hh');