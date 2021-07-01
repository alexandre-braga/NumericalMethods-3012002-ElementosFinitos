clc
clear all
close all

format long;
%dominio
a = 0.00;
b = 2.00;

c = 0.00;
d = 1.25;

erro = zeros(4,1);
hh = zeros(4,1);

for grau = 1:1
  for cont = 1:1
    Kappa = 1.;
    E = 1.e-2;
    Peh = 5;
    %if cont >=2
    %  Peh = 5;
    %  if cont >=3
    %    Peh = 10;
    %  endif
    %endif
    h = (Peh*2*E)/abs(Kappa);
    deltaT = h^2;
    nel = (b-a)/h;
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
    xl = zeros(np,1);
    xl(1) = a;
    for i = 2:np
      xl(i) = xl(i-1) + h/k;
    endfor

    %gera shg e pega as funções peso
    [shg, w]= shgGera(nen,nint);

##    %montagem global
##    for n = 1:nel
##      for l = 1:nint
##        xx = 0.;
##        for i = 1:nen
##          xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
##        endfor
##        for j = 1:nen
##          F(k*(n-1)+j) = F(k*(n-1)+j) + funcao(xx)*shg(1,j,l)*w(l)*h/2;
##          for i = 1:nen
##            K((k*(n-1)+i),(k*(n-1)+j)) = K((k*(n-1)+i),(k*(n-1)+j)) + funcaok(xx,E)*shg(2,i,l)*2/h*shg(2,j,l)*2/h*w(l)*h/2;
##            C((k*(n-1)+j),(k*(n-1)+i)) = C((k*(n-1)+j),(k*(n-1)+i)) + funcaoKappa(xx,Kappa)*shg(2,i,l)*2/h*shg(1,j,l)*w(l)*h/2;
##          endfor
##        endfor
##      endfor
##    endfor
##    KC = K + C;
##    %Condição de Dirichlet entrada
##    KC(1,1) = 1;
##    F(1) = funcaoExata(0,t);
##    for i = 2:k+1
##      F(i) = F(i) - (F(1)*KC(i,1));
##      KC(1,i) = 0;
##      KC(i,1) = 0;
##    endfor
##    %Condição de Dirichlet saida
##    for i = np-(k+1):np
##      F(i) = F(i) - (F(np)*KC(i,np));
##      KC(np,i) = 0.;
##      KC(i,np) = 0.;
##    endfor
##    KC(np,np) = 1;
##    F(np) = funcaoExata(np,t);
    
    %função exata
    x = a;
    t = c;
    j = 1;
    espacoT = ceil((d-c)/deltaT + 1)
    exata = zeros(np,espacoT);
    while t <= d
      for i = 1:np
        exata(i,j) = funcaoExata(x,t,E,Kappa);
        x += h/k;
      endfor
      x = a;
      j++;
      t += deltaT;
    endwhile
    j
    x = a:h/k:b;
    t = c:deltaT:d;
    u = zeros(np);
    u = KC\F;
    
    %salva a resolucao
    nome = sprintf("log/PesosEPontosIntegracaoPeh%dGrau%d.txt", Peh, grau);
    save(nome, 'xl', 'h', 'deltaT', 'u', 'x', 't', 'exata');
      
  endfor 
  
endfor 
