clear all;
close all;
x= logspace(-2, 0, 11);
% y6=[0.0361371        0.0360592      0.0362702    0.0367406   0.0377336  0.0401291];
% y12=[0.005049          0.006228        0.005367    0.005186       0.005080    0.008293];
% y18=[0.000450         0.000458         0.000497     0.000628     0.001127    0.002940];

y6=[0.0361371        0.0360592      0.0362702    0.0367406   0.0377336  0.0401291  0.047299   0.064290    0.105701  0.190095    0.499974];
y12=[0.005049          0.005080        0.005186   0.005367   0.006228   0.008293 0.014733   0.034106    0.081284  0.176480    0.499746];
y18=[0.000450         0.000458         0.000497     0.000628     0.001127    0.002940 0.009699   0.030835    0.082598  0.181904    0.499221];


figure()
box on; hold on;
plot(x,y6,'bo-');
plot(x,y12,'bo-');
plot(x,y18,'bo-');
xlabel('ρ ');
ylabel('BER');
legend('CSI AFDM/MMSE BPSK')
grid on;
set(gca,'Yscale','log');
set(gca,'Xscale','log');
ylim([1e-6 1]);
xlim([1e-2,  1e-0]); 
xlabel('ρ ');
ylabel('BER');
legend('CSI AFDM/MMSE BPSK')
grid on;


% figure;
% semilogx(x, y18, 'LineWidth', 2); %  
% grid on;
% xlim([1e-2, 1e-1]); %
% set(gca,'Yscale','log');
% ylim([1e-6 1]);
% ylabel('BER');
% xticks([1e-2, 1e-1]);
% xticklabels({'10^{-2}', '10^{-1}'});
% xlabel('ρ ');


% x= logspace(-2, -1, 6);
% 
% figure;
% semilogx(x, y18, 'LineWidth', 2); %  
% grid on;

% x_ticks = [10^(-2), 10^(-1)];
% xticks(x_ticks);
% xticklabels({'10^{-2}', '10^{-1}'});
% y_ticks = [10^(-6), 10^(-5),10^(-4),10^(-3),10^(-2),10^(-1), 10^(0)];
% xticks(y_ticks);
% xticklabels({'10^{-6}', '10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'});