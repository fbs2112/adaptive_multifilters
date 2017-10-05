clear;
clc;
close all;


addpath(['.' filesep 'results']);
addpath(['..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);


fileVector = [1 5 13 9];



for i = 4:5
    
    figure
    
    for j = 1:length(fileVector)

        if fileVector(j) < 10
            load (['resultsBER0' num2str(fileVector(j)) '.mat']);
        else
            load (['resultsBER' num2str(fileVector(j)) '.mat']);
        end

        semilogy(SNR,squeeze(ber(i,:)))
        hold on;
    end
    H = legend('V-PNLMS','VSM-PNLMS','V-RLS','VM-BEACON');
    set(H,'location','SouthWest');
    set(H,'interpreter','latex')
    
    xlabel('SNR [dB]','interpreter','latex');
    ylim([1e-6 1e0])
    xlim([0 30]);
    
    set(gca,'xtick',SNR);
  
    if i == 4
        ylabel('BER','interpreter','latex');
    else
        set(gca,'ytick',[]);
    end
        
    
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'berFF'  num2str(i)],'en' , figProp );
    
end

close all;
fileVector = [3 7 15 11];

figure

for i = 5:5
    for k = 1:1
    
        for j = 1:length(fileVector)
            
            if fileVector(j) > 9
                load (['resultsBER' num2str(fileVector(j)) '.mat']);
            else
                load (['resultsBER0' num2str(fileVector(j)) '.mat']);
            end

            semilogy(SNR,squeeze(ber(i,k,:)))
            hold on;
        end
        H = legend('V-PNLMS','VSM-PNLMS','V-RLS','VM-BEACON');
        set(H,'location','SouthWest');
        set(H,'interpreter','latex')

        xlabel('SNR [dB]','interpreter','latex');
        ylim([1e-6 1e0])
        xlim([0 30]);

        set(gca,'xtick',SNR);

        ylabel('BER','interpreter','latex');
        
    end
    
    formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'berDFE'],'en' , figProp );
    
end





