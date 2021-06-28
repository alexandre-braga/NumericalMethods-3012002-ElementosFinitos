clc
clear all
close all
%salva os erros
format long;
%dominio
a = 0.0;
b = 1.0;
erro = zeros(4,1);
hh = zeros(4,1);

for grau = 1:4
  for cont = 1:4
    if cont < 3
      E = 1.e-3;
    else
      E = 1.e-4;
    endif
    %n de elementos
    nel = 4*(cont);
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

    %montagem do xl !!!LINEAR!!!
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
            K((k*(n-1)+i),(k*(n-1)+j)) = K((k*(n-1)+i),(k*(n-1)+j)) + funcaok(xx,E)*shg(2,i,l)*2/h*shg(2,j,l)*2/h*w(l)*h/2;
            M((k*(n-1)+i),(k*(n-1)+j)) = M((k*(n-1)+i),(k*(n-1)+j)) + funcaoY(xx)*shg(1,i,l)*shg(1,j,l)*w(l)*h/2;
          endfor
        endfor
      endfor
    endfor
    
    KM = K + M;
    
    %Condição de Dirichlet entrada
    KM(1,1) = 1;
    F(1) = 0;
    for i = 2:k+1
      F(i) = F(i) - (F(1)*KM(i,1));
      KM(1,i) = 0;
      KM(i,1) = 0;
    endfor
    
    %Condição de Dirichlet saida
    for i = np-(k+1):np
      F(i) = F(i) - (F(np)*KM(i,np));
      KM(np,i) = 0.;
      KM(i,np) = 0.;
    endfor
    KM(np,np) = 1;
    F(np) = 0;
    
    %função exata
    x = a;
    exata = zeros(np,1);
    for i = 1:np
      exata(i) = funcaoExata(x,E);
      x += h/k;
    endfor  
    x = a:h/k:b;
    u = zeros(np);
    u = KM\F;

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
    nome = sprintf("log/PesosEPontosIntegrecao%dGrau%d.txt", cont, grau);
    save(nome, 'xl', 'h', 'u', 'x', 'exata');
    
  
  endfor 
  
  %verifica estabilidade
  for i = 1:4
    if hh(i) < sqrt(6*E)
      estabilidade(i) = true;
    else
      estabilidade(i) = false;
      difestabilidade(i) = hh(i) - sqrt(6*E);
    endif
  endfor
  %salva os erros
  nome = sprintf("log/Erros%d.txt", grau);
  save(nome, 'erro', 'hh', 'estabilidade', 'difestabilidade');

endfor 
