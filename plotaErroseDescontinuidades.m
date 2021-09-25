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

plot(x,exata); hold on; legend('exata'); hold on; j = 1; i=1; while i<=256 plot(xl(j:j+(nen-1)),U(i,:)); hold on; j+=(nen-1); i++; endwhile

load log/erros/Erros1.txt;  plot (-log10(hh), log10(derro)); hold on; load log/erros/Erros2.txt;  plot (-log10(hh), log10(derro)); hold on; load log/erros/Erros3.txt;  plot (-log10(hh), log10(derro)); hold on; load log/erros/Erros4.txt;  plot (-log10(hh), log10(derro)); hold on; xlabel -log10(hh); ylabel log10(derro); legend('grau1', 'grau2', 'grau3', 'grau4');
