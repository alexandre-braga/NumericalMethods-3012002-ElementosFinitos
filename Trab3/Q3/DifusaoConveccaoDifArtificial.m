clc
clear all
close all

format long;
%dominio
a = 0.00;
b = 2.00;

T0 = 0.00;
T = 1.25;

erro = zeros(4,1);
hh = zeros(4,1);

for grau = 1:1
  for cont = 1:3
    Kappa = 1.;
    E = 1.e-2;
    Peh = 1;
    if cont >=2
      Peh = 5;
      if cont >=3
        Peh = 10;
      endif
    endif
    h = (Peh*2*E)/abs(Kappa);
    Beta = 1+(1 - 1/Peh);
    Gamma = Beta*E*Peh;
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
    %matriz global de massa zerada
    M = zeros(np,np);
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
            M((k*(n-1)+i),(k*(n-1)+j)) = M((k*(n-1)+i),(k*(n-1)+j)) + shg(1,i,l)*shg(1,j,l)*w(l)*h/2*1/deltaT;
            K((k*(n-1)+i),(k*(n-1)+j)) = K((k*(n-1)+i),(k*(n-1)+j)) + (funcaok(xx,E) + Gamma)*shg(2,i,l)*2/h*shg(2,j,l)*2/h*w(l)*h/2;
            C((k*(n-1)+j),(k*(n-1)+i)) = C((k*(n-1)+j),(k*(n-1)+i)) + funcaoKappa(xx,Kappa)*shg(2,i,l)*2/h*shg(1,j,l)*w(l)*h/2;
          endfor
        endfor
      endfor
    endfor
    
    %U zero  = phiX(x,0)
    unext = zeros(np,1);
    for i = 1:np
      unext(i) = funcaoExata(xl(i),0,E,Kappa);
    endfor
    
    n = 0;
    t = T0;
    %tam aprox pra dim t
    espacoT = ceil((T-T0)/deltaT + 1)
    U = zeros(np,espacoT);
    while (t < T)
        t += deltaT;
        n++;
        if n+1 > espacoT + 1
          unext
          break;
        endif
        %remontagem da matriz pra refazer  o contorno
        A = zeros(np,np);
        A = M + K + C;
        Fonte = zeros(np,1);
        Fonte = M*unext + F;
        
        %Condição de Dirichlet entrada
        A(1,1) = 1.;
        Fonte(1) = funcaoExata(a,t,E,Kappa);
        for i = 2:k+1
          Fonte(i) = Fonte(i) - (Fonte(1)*A(i,1));
          A(1,i) = 0.;
          A(i,1) = 0.;
        endfor
        
       %Condição de Dirichlet saida
       Fonte(np) = funcaoExata(b,t,E,Kappa);
       for i = np-(k+1):np
         Fonte(i) = Fonte(i) - (Fonte(np)*A(i,np));
         A(np,i) = 0.;
         A(i,np) = 0.;
       endfor
       A(np,np) = 1.;
       Fonte(np) = funcaoExata(b,t,E,Kappa);
        
       unext =  A\Fonte;
       
       %armazena novo unext na matriz solucao
       for i = 1:np
         U(i,n+1) = unext(i);
       endfor
       
    endwhile

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
    j
    x = a:h/k:b;
    t = T0:deltaT:T;
    
    %ajusta tam U pra plotagem
    if length(t) < columns(U)
      U(:,length(t)+1) = [];
    endif
    
    %ajusta tam exata pra plotagem
    if length(t) < columns(exata)
      exata(:,length(t)+1) = [];
    endif
    
    
    %salva a resolucao
    nome = sprintf("log/PesosEPontosIntegracaoDAPeh%dGrau%d.txt", Peh, grau);
    save(nome, 'Gamma', 'xl', 'h', 'deltaT', 'U', 'x', 't', 'exata');
      
  endfor 
  
endfor 
