clc
clear all
close all
%dominio
a = -2.0;
b = 2.0;
erro = zeros(4,1);
hh = zeros(4,1);

for grau = 1:4
  for cont = 2:5
    %n de elementos
    nel = 4^(cont);
    %tamanho do elemento
    h = (b-a)/nel
    %grau do polinomio de interpolação
    k = grau
    %n de nós do elemento
    nen = k+1;
    %n de nós global
    np =  k*nel+1;
    %n de pontos de integração
    nint = nen;
    
    %montagem do xl !!!LINEAR!!!
    xl = zeros(np,1);
    xl(1) = a;
    for i = 2:np
      xl(i) = xl(i-1) + h/k;
    endfor
    
    %gera shg e pega as funções peso
    [shg, w]= shgGera(nen,nint);

    %sem montagem global
    uu = zeros(nel,2);
    for n = 1:nel
      Ke = zeros(nen,nen);
      Fe = zeros(nen,1);
      u = zeros(nen);
      for l = 1:nint
        xx = 0;
        for i = 1:nen
          xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
        endfor
        for j = 1:nen
          Fe(j) = Fe(j) + funcao(xx)*shg(1,j,l)*w(l)*h/2;
          for i = 1:nen
            Ke(i,j) = Ke(i,j) + shg(1,i,l)*shg(1,j,l)*w(l)*h/2;
          endfor
        endfor
      endfor
      u = Ke\Fe;
      uu(n,1) = u(1);
      uu(n,2) = u(2);    
    endfor
    
    %função exata
    x = a;
    exata = zeros(np,1);
    for i = 1:np
      exata(i) = funcao(x);
      x += h/k;
    endfor  
    x = a:h/k:b;

    %cálculo do erro
    %erul2 = 0;
    %for n = 1:nel
    %  eru = 0;
    %  for l = 1:nint
    %    uh = 0;
    %    xx = 0;
    %    for i = 1:nen
    %      %uh !!!LINEAR!!!
    %      uh = uh + shg(1,i,l)*u(k*(n-1)+i);
    %      %xl !!!LINEAR!!!
    %      xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
    %    endfor
    %    eru = eru + ((funcao(xx) - uh)**2) * w(l) * h/2;
    %  endfor
    %  erul2 = erul2 + eru;
    %endfor
    %erul2 = sqrt(erul2);
    %erro(cont) = erul2;
    %hh(cont) = h;
    
    %salva a resolucao
    nome = sprintf("log/PesosEPontosIntegrecaoDescontinua%dGrau%d.txt", cont-1, grau);
    save(nome, 'uu', 'x', 'h', 'exata', 'xl');
    
  endfor
  
  %salva os erros
  %nome = sprintf("log/ErrosDescontinua%d.txt", grau);
  %save(nome, 'erro', 'hh');

endfor 