clear;
clc;
close all;



addpath(['.' filesep 'results']);
addpath(['..' filesep '..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);



linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

fileVector = {'resultsTest2'};


for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    figure;
    for i = 1:size(e3,1)
        
        
        x = e3{i,1};
        aux = find(x,1);
        xAux = 10*log10(x(aux:end));
        plot(xAux)
        
        hold on
        updatesLin = mean(meanCountLin{i,1}(aux:999))*100;
        updatesNonLin = mean(meanCountNonLin{i,1}(aux:999))*100;
        
        updatesAux(i,l,:) = [updatesLin;updatesNonLin];
        
    end
    
    H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
        '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
        '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
    set(H,'interpreter','latex')
    ylim([-35 10]);
    
    ylabel('MSE [dB]','interpreter','latex');
    xlabel('Iterations [$k$]','interpreter','latex');
    xlim([0 3000])
end