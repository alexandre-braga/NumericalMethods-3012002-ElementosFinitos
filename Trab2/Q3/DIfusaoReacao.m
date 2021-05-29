1
%dominio
a = 0.0;
b = 1.0;
%n de elementos
nel = 16;
%constante de difusão
E = 10e-3;
%constantes da solucão exata
c2 = (exp(-1/sqrt(E)) - 1) / ( exp(1/sqrt(E)) - exp(-1/sqrt(E)) );
c1 = -1 - c2;

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
%matriz de rigidez global zerada
K = zeros(np,np);
%matriz de massa global zerada
M = zeros(np,np);
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
        %matriz de rigidez - derivadas phi i e phi j
        K((k*(n-1)+i),(k*(n-1)+j)) = K((k*(n-1)+i),(k*(n-1)+j)) + funcaok(xx)*shg(2,i,l)*2/h*shg(2,j,l)*2/h*w(l)*h/2;
        %matriz de massa - phi i e phi j
        M((k*(n-1)+i),(k*(n-1)+j)) = M((k*(n-1)+i),(k*(n-1)+j)) + shg(1,i,l)*shg(1,j,l)*w(l)*h/2;
      endfor
    endfor
  endfor
endfor

%solução exata
x = a;
i = 1;
while x != b+h
  exata(i) = c1*exp(-x/sqrt(E)) + c2*exp(x/sqrt(E)) + 1;
  i++;
  x += h;
endwhile
x = a:h:b;
u = zeros(np);
DIFREAC = K+M;

%Condição de Dirichlet
DIFREAC(1,1) = 1;
DIFREAC(1,2) = 0;
F(1) = 0;
F(2) = F(2) - (1*DIFREAC(2,1));
DIFREAC(2,1) = 0;

%Condição de Dirichlet
DIFREAC(nel+1,nel+1) = 1;
DIFREAC(nel+1,nel) = 0;
F(nel+1) = 0;
F(nel) = F(nel) - (1*DIFREAC(nel,nel+1));
DIFREAC(nel,nel+1) = 0;

u = DIFREAC\F

save saida.txt xl u x exata;
%figure;
%plot(xl,u,x,exata);