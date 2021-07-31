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
    nel = 16*cont;
    h = (b-a)/nel;
    k = grau;
    nen = k+1;
    np =  k*nel+1
    nint = nen;
    
    A = zeros(np,np);
    B = zeros(np,np);
    BT = zeros(np,np);
    C = zeros(np,np);
    
    G = zeros(np,1);
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
          Ge(k*(n-1)+j) = Ge(k*(n-1)+j) - 0. + delta2*funcao(xx)*shg(2,j,l)*2/h*w(l)*h/2;
          Fe(k*(n-1)+j) = Fe(k*(n-1)+j) - funcao(xx)*shg(1,j,l)*w(l)*h/2;
          for i = 1:nen
            Ae((k*(n-1)+i),(k*(n-1)+j)) = Ae((k*(n-1)+i),(k*(n-1)+j)) + (shg(1,i,l)*shg(1,j,l) + delta1*shg(1,i,l)*shg(1,j,l) + delta2*shg(2,i,l)*2/h*shg(2,j,l)*2/h)*w(l)*h/2;
            Be((k*(n-1)+i),(k*(n-1)+j)) = Be((k*(n-1)+i),(k*(n-1)+j)) + (-shg(1,i,l)*shg(2,j,l)*2/h + delta1*shg(2,i,l)*2/h*shg(1,j,l))*w(l)*h/2;
            Ce((k*(n-1)+i),(k*(n-1)+j)) = Ce((k*(n-1)+i),(k*(n-1)+j)) + (delta1*shg(2,i,l)*2/h*shg(2,j,l)*2/h)*w(l)*h/2;
          endfor
        endfor
      endfor
      A = A + Ae;
      B = B + Be;
      BT = B';
      C = C + Ce;
      
      G = G + Ge;
      F = F + Fe;
    endfor
    %Matriz global
    M = zeros(2*np,2*np);
    Fonte = zeros(2*np,1);
    
    #Condição de Dirichlet entrada forçando p
    C(1,1) = 1.;
    BT(1,1) = 0.;
    B(1,1) = 0.;
    F(1) = funcaoExata(a);
    for i = 2:k+1
      F(i) = F(i) - (F(1)*C(i,1));
      F(i) = F(i) - (F(1)*BT(i,1));
      G(i) = G(i) - (G(1)*B(i,1));
      C(1,i) = 0.;
      BT(1,i) = 0.;
      B(1,i) = 0.;
      C(i,1) = 0.;
      BT(i,1) = 0.;
      B(i,1) = 0.;
    endfor
        
    #Condição de Dirichlet saida forçando p
    F(np) = funcaoExata(b);
    for i = np -(k+1):np
      F(i) = F(i) - (F(np)*C(i,np));
      F(i) = F(i) - (F(np)*BT(i,np));
      G(i) = G(i) - (G(np)*B(i,np));
      C(np,i) = 0.;
      BT(np,i) = 0.;
      B(np,i) = 0.;
      C(i,np) = 0.;
      BT(i,np) = 0.;
      B(i,np) = 0.;
    endfor
    C(np,np) = 1.;
    BT(np,np) = 0.;
    B(np,np) = 0.;
    F(np) = funcaoExata(b);
    
    
    %Caso o tamanho difira, usar outros loops
    %E modificar o tamanho da M, o mesmo pra F
    for i = 1:np
      for j = 1:np
        M(i,j) = A(i,j);
        M(i,j + np) = B(i,j);
        M(i + np,j) = BT(i,j);
        M(i + np,j + np) = C(i,j);
      endfor
      Fonte(i) = G(i);
      Fonte(i + np) = F(i);
    endfor

    U = M\Fonte;

    u = zeros(np,1);
    p = zeros(np,1);
    for i = 1:np
      u(i) = U(i);
      p(i) = U(i + np);
    endfor
    
    %função exata
    x = a;
    exata = zeros(np,1);
    for i = 1:np
      exata(i) = funcaoExata(x);
      x += h/k;
    endfor
    x = a:h/k:b;
    
    %salva a resolucao
    nome = sprintf("log/PesosEPontosIntegracao%dGrau%d.txt", cont, grau);
    save(nome, 'xl', 'h', 'u', 'p', 'x', 'exata');
      
  endfor 
  
endfor 
