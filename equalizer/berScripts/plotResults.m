clear;
clc;
close all;




addpath(['.' filesep 'results']);



load testSM_PAPA_LinEqBer.mat;
figure
for i = 1:size(ber,2)
    semilogy(SNR,ber(:,i))
    hold on
end
xlabel('SNR [dB]','interpreter','latex');
ylabel('BER','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([min(SNR) max(SNR)]);
ylim([1e-8 1]);




load testSM_PAPA_VolterraEqBer.mat;
figure
for i = 1:size(ber,2)
    semilogy(SNR,ber(:,i))
    hold on
end
xlabel('SNR [dB]','interpreter','latex');
ylabel('BER','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([min(SNR) max(SNR)]);
ylim([1e-8 1]);




load testSM_PAPA_DFELinBer.mat;
figure
for i = 1:size(ber,2)
    semilogy(SNR,ber(:,i))
    hold on
end
xlabel('SNR [dB]','interpreter','latex');
ylabel('BER','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([min(SNR) max(SNR)]);
ylim([1e-8 1]);



load testSM_PAPA_DFEVolterraBer.mat;
figure
for i = 1:size(ber,2)
    semilogy(SNR,ber(:,i))
    hold on
end
xlabel('SNR [dB]','interpreter','latex');
ylabel('BER','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([min(SNR) max(SNR)]);
ylim([1e-8 1]);







rmpath(['.' filesep 'results']);
