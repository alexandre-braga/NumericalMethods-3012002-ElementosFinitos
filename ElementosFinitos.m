1
%dominio
a = -2.0;
b = 2.0;
%n de elementos
nel = 16;
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
  %matriz do elemento
  Ke = zeros(nen,nen);
  %vetor fonte do elemento
  Fe = zeros(nen,1);
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

%função exata
x = a:h:b;
exata = funcao(x);

u = zeros(np);
u = K\F;

save saida.txt xl u x exata;
%figure;
%plot(xl,u,x,exata);