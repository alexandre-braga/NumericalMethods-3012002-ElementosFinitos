clc
clear all
close all

format long;
%dominio
a = 0.;
b = 1.;
erro = zeros(4,1);
hh = zeros(4,1);

for grau = 1:4
  for cont = 2:5
    %n de elementos
    nel = 4^(cont);
    %tamanho do elemento
    h = (b-a)/nel
    %grau do polinomio de interpolação
    k = grau;
    %n de nós do elemento
    nen = k+1;
    %n de nós global
    np =  k*nel+1;
    %n de pontos de integração
    nint = nen;
    %matriz global de massa zerada
    K = zeros(np,np);
    %matriz global de rigidez zerada
    M = zeros(np,np);
    %matriz global de difusao-reação zerada
    KM = zeros(np,np);
    %vetor fonte zerado
    F = zeros(np,1);

    %montagem do xl
    xl = zeros(np,1);
    xl(1) = a;
    for i = 2:np
      xl(i) = xl(i-1) + h/k;
    endfor
    %gera shg e pega as funções peso
    [shg, w]= shgGera(nen,nint);

    %montagem global
    for n = 1:nel
      for l = 1:nint
        xx = 0.;
        for i = 1:nen
          xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
        endfor
        for j = 1:nen
          F(k*(n-1)+j) = F(k*(n-1)+j) + funcao(xx)*shg(1,j,l)*w(l)*h/2;
          for i = 1:nen
            K((k*(n-1)+i),(k*(n-1)+j)) = K((k*(n-1)+i),(k*(n-1)+j)) + funcaok(xx)*shg(2,i,l)*2/h*shg(2,j,l)*2/h*w(l)*h/2;
            M((k*(n-1)+i),(k*(n-1)+j)) = M((k*(n-1)+i),(k*(n-1)+j)) + funcaoY(xx)*shg(1,i,l)*shg(1,j,l)*w(l)*h/2;
          endfor
        endfor
      endfor
    endfor

    KM = K + M;
    
    %Condição de Dirichlet entrada
    KM(1,1) = 1;
    F(1) = 1;
    for i = 2:k+1
      F(i) = F(i) - (1*KM(i,1));
      KM(1,i) = 0;
      KM(i,1) = 0;
    endfor
    
    %Condição de Robin saida
    
    KM(np,np) += 10e6*0.5;
    F(np) += 10e6*0.5+5/2;
    %função exata
    x = a;
    exata = zeros(np,1);
    for i = 1:np
      exata(i) = funcaoExata(x);
      x += h/k;
    endfor  
    x = a:h/k:b;
    u = zeros(np);
    u = KM\F;

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
      nome = sprintf("log/PesosEPontosIntegrecao%dGrau%d.txt", cont-1, grau);
      save(nome, 'xl', 'h', 'u', 'x', 'exata');
      
  endfor 
  
  %salva os erros
  nome = sprintf("log/Erros%d.txt", grau);
  save(nome, 'erro', 'hh');

endfor 