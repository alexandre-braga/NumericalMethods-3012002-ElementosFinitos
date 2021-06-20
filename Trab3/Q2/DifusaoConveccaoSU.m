clc
clear all
close all

format long;
%dominio
a = 0.0;
b = 1.0;
erroSU = zeros(1001,1);
hhSU = zeros(1001,1);
betasSU = zeros(1001,1);

for grau = 1:1
  Beta = 0.0;
  for cont = 1:1001
    Kappa = 1.;
    E = 1.e-2;
    %n de elementos
    nel = 10;
    %tamanho do elemento
    hSU = (b-a)/nel;
    %grau do polinomio de interpolação
    k = grau;
    %n de nós do elemento
    nen = k+1;
    %n de nós global
    np =  k*nel+1;
    %n de pontos de integração
    nint = nen;
    %matriz global de rigidez zerada
    K = zeros(np,np);
    %matriz global de convecção zerada
    C = zeros(np,np);
    %matriz global de difusao-convecção zerada
    KC = zeros(np,np);
    %vetor fonte zerado
    F = zeros(np,1);

    %montagem do xl
    xlSU = zeros(np,1);
    xlSU(1) = a;
    for i = 2:np
      xlSU(i) = xlSU(i-1) + hSU/k;
    endfor

    %gera shg e pega as funções peso
    [shg, w]= shgGera(nen,nint);
    [shgSU, wSU]= shgGeraSU(nen,nint,Beta);
        
    %montagem global
    for n = 1:nel
      for l = 1:nint
        xx = 0.;
        for i = 1:nen
          xx = xx + shg(1,i,l)*xlSU(k*(n-1)+i);
        endfor
        for j = 1:nen
          F(k*(n-1)+j) = F(k*(n-1)+j) + funcao(xx)*shgSU(1,j,l)*wSU(l)*hSU/2;
          for i = 1:nen
            K((k*(n-1)+i),(k*(n-1)+j)) = K((k*(n-1)+i),(k*(n-1)+j)) + funcaok(xx,E)*shg(2,j,l)*2/hSU*shgSU(2,i,l)*2/hSU*w(l)*hSU/2;
            C((k*(n-1)+i),(k*(n-1)+j)) = C((k*(n-1)+i),(k*(n-1)+j)) + funcaoKappa(xx,Kappa)*shg(2,j,l)*2/hSU*shgSU(1,i,l)*w(l)*hSU/2;
          endfor
        endfor
      endfor
    endfor
    KC = K + C;
    
    %Condição de Dirichlet entrada
    KC(1,1) = 1;
    F(1) = 0;
    for i = 2:k+1
      F(i) = F(i) - (F(1)*KC(i,1));
      KC(1,i) = 0;
      KC(i,1) = 0;
    endfor
    
    %Condição de Dirichlet saida
    for i = np-(k+1):np
      F(i) = F(i) - (F(np)*KC(i,np));
      KC(np,i) = 0.;
      KC(i,np) = 0.;
    endfor
    KC(np,np) = 1;
    F(np) = 1;
    
    %função exata
    x = a;
    exata = zeros(np,1);
    for i = 1:np
      exata(i) = funcaoExata(x);
      x += hSU/k;
    endfor  
    x = a:hSU/k:b;
    uSU = zeros(np);
    uSU = KC\F;

    %cálculo do erro
    erdul2 = 0;
    for n = 1:nel
      erdu = 0;
      for l = 1:nint
        duh = 0;
        xx = 0;
        for i = 1:nen
          duh = duh + shgSU(2,i,l)*2/hSU*uSU(k*(n-1)+i);
          xx = xx + shgSU(1,i,l)*xlSU(k*(n-1)+i);
        endfor
        erdu = erdu + ((dfuncaoExata(xx) - duh)**2) * w(l) * hSU/2;
       endfor
       erdul2 = erdul2 + erdu;
     endfor
    erdul2 = sqrt(erdul2);
    erroSU(cont) = erdul2;
    hhSU(cont) = hSU;
    betasSU(cont) = Beta;
    Beta += 0.001;
      
  endfor 
  
  %salva os erros
  nome = sprintf("log/ErrosBetaSU.txt");
  save(nome, 'betasSU', 'erroSU', 'hhSU');

endfor 
