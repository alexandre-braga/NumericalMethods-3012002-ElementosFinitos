clear;
load log/erros/Erros1.txt;
plot (-log10(hh),log10(erro));
hold on;

load log/erros/Erros2.txt;
plot (-log10(hh),log10(erro));
hold on;

load log/erros/Erros3.txt;
plot (-log10(hh),log10(erro));
hold on;

load log/erros/Erros4.txt;
plot (-log10(hh),log10(erro));
hold on;

xlabel -log10(hh);
ylabel log10(erro); 
legend('grau1', 'grau2', 'grau3', 'grau4');

plot(x,exata,xlU1,U(:,1),xlU2,U(:,2)); legend('exata','U1','U2');
