clc
clear all
close all

format long;
%dominio
a = 0.00;
b = 1.00;

grau = 1;
beta0 = 10;

erro = zeros(4,1);
hh = zeros(4,1);

for grau = 1:1
  for alfa = -1:1
    for cont = 1:4
      nel = 4^cont;
      h = (b-a)/nel;
      k = grau;
      
      beta = k*k*beta0/h
      
      nen = k+1;
      np =  k*nel+1
      nint = nen;
      
      Lambda = zeros(np);
      lambdak = zeros(nen);
      U = zeros(nel,nen);
      u = zeros(nen);

      elementosK = zeros(np,np);
      FkGLobal = zeros(np,1);
      
      %montagem do xl
      xl = zeros(np,1);
      xl(1) = a;
      for i = 2:np
        xl(i) = xl(i-1) + h/k;
      endfor
      
      %gera shg e pega as funções peso
      [shg, w]= shgGera(nen,nint);
      [shge]= shgeGera(nen,nint);
      
      %Problema Global
      for n = 1:nel      
        Ak = zeros(nen,nen);
        Bk = zeros(nen,2);
        BTk = zeros(2,nen);
        Ck = zeros(2,2);
        Fk = zeros(nen);
        
        #parte integral de A e F no elemento
        for l = 1:nint
          xx = 0.;
          for j = 1:nen
            xx = xx + shg(1,j,l)*xl(k*(n-1)+j);
          endfor
          for i = 1:nen
            Fk(i) = Fk(i) + funcao(xx)*shg(1,i,l)*w(l)*h/2;
            for j = 1:nen
              Ak(i,j) = Ak(i,j) + shg(2,j,l)*2/h*shg(2,i,l)*2/h**w(l)*h/2;
            endfor
          endfor
        endfor
        
        #Calcula B, A , BT, C pra cada elemento
        for i = 1:nen
          Bk(i,1) += alfa*shge(2,i,1)*2/h + beta*shge(1,i,1);
          Bk(i,2) += -alfa*shge(2,i,2)*2/h - beta*shge(1,i,2);
          for j = 1:nen
            Ak(i,j) += -(shge(2,j,2)*2/h*shge(1,i,2) - shge(2,j,1)*2/h*shge(1,i,1)) + alfa*(shge(2,i,2)*2/h*shge(1,j,2) - shge(2,i,1)*2/h*shge(1,j,1)) + beta*(shge(1,j,2)*shge(1,i,2) - shge(1,j,1)*shge(1,i,1));
          endfor
        endfor
        for j = 1:nen
          BTk(1,j) += alfa*shge(2,j,1)*2/h + beta*shge(1,j,1);
          BTk(2,j) += -alfa*shge(2,j,2)*2/h - beta*shge(1,j,2);
        endfor
        Ck(1,1) = -beta;
        Ck(2,2) = beta;
        
        #calcula K e F elemento
        elementoK = zeros(2,2);
        elementoK = Ck - BTk*inverse(Ak)*Bk;
        Fkk = zeros(2);
        Fkk = -BTk*inverse(Ak)*Fk;
        
        #armazena como global o elemento
        for i = 1:nen
          FkGLobal(k*(n-1)+i) += Fkk(i);
          for j = 1:nen
            elementosK(k*(n-1)+i,k*(n-1)+j) += elementoK(1,j);
          endfor
        endfor
        
      endfor
      
      Lambda = elementosK\FkGLobal
      
      %Problema Local
      for n = 1:nel
        Ae = zeros(nen,nen);
        Fe = zeros(nen,1);
        for l = 1:nint
          xx = 0.;
          for j = 1:nen
            xx = xx + shg(1,j,l)*xl(k*(n-1)+j);
          endfor
          for i = 1:nen
            Fe(i) = Fe(i) + funcao(xx)*shg(1,i,l)*w(l)*h/2;
            for j = 1:nen
              Ae(i,j) = Ae(i,j) + shg(2,j,l)*2/h*shg(2,i,l)*2/h**w(l)*h/2;
            endfor
          endfor
        endfor
        
        for i = 1:nen
          Fe(i) += alfa*(shge(2,i,2)*2/h*Lambda(n+1) - shge(2,i,1)*2/h*Lambda(n) ) + beta*(shge(1,i,2)*Lambda(n+1) - shge(1,i,1)*Lambda(n));
          for j = 1:nen
            Ae(i,j) += -(shge(2,j,2)*2/h*shge(1,i,2) - shge(2,j,1)*2/h*shge(1,i,1)) + alfa*(shge(2,i,2)*2/h*shge(1,j,2) - shge(2,i,1)*2/h*shge(1,j,1)) + beta*(shge(1,j,2)*shge(1,i,2) - shge(1,j,1)*shge(1,i,1));
          endfor
        endfor
        
        u = Ae\Fe;
        #pos = nen*(n-1);
        U(n,1) = u(1);
        U(n,2) = u(2);
        
      endfor
      
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
            uh = uh + shg(1,i,l)*u(i);
            xx = xx + shg(1,i,l)*xl(k*(n-1)+i);
          endfor
          eru = eru + ((funcaoExata(xx) - uh)**2) * w(l) * h/2;
        endfor
        erul2 = erul2 + eru;
      endfor
      erul2 = sqrt(erul2);
      erro(cont) = erul2;
      hh(cont) = h;
      
      xlU1 = xl(1:end-1);
      xlU2 = xl(2:end);
      
      %salva a resolucao
      nome = sprintf("log/PesosEPontosIntegracao%dalfa%dgrau%d.txt", cont, alfa, grau);
      save(nome, 'alfa', 'beta', 'h', 'xlU1', 'xlU2', 'U', 'x', 'exata');
      
    endfor 
    
    %salva os erros
    nome = sprintf("log/Erros%d.txt", alfa);
    save(nome, 'erro', 'hh' );
    
    alfa++;
    
  endfor 
endfor