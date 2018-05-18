clear;
clc;
close all;



addpath(['.' filesep 'results']);
addpath(['..' filesep '..' filesep 'Misc' filesep]); 



linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  

fileVector = {'results_SMDTPU_03'};
N = 2;

M = 6:8;
numberOfIterations = 7000;
upVector = [47.59 40.18];

for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    
    for i = 3:size(e4,1)
        linSamples = M(i) + 1;
        nonLinSamples = nchoosek(M(i)+1,N) + M(i) + 1;
        for j = 1:size(e4{i,1},1)
            figure;
            
            for k = 1:size(e4{i,1},2)
                
                x = e4{i,1}{j,k};
                aux = find(x,1);
                xAux = 10*log10(x(aux:end));
                plot(xAux)
                
                hold on
                
                meanMSE = mean(10*log10((e4{i,1}{j,k}(end - 999:end))));
                stdMse = std(10*log10((e4{i,1}{j,k}(end - 999:end))));
                convergenceSample(i,k) = find(10*log10(e4{i,1}{j,k}(aux:end)) < meanMSE + 2*stdMse*sign(2*stdMse),1,'first');
                
                updatesLin = mean(meanCountLin2{i,1}{j,k}(aux:end))*100;
                updatesNonLin = mean(meanCountNonLin2{i,1}{j,k}(aux:end))*100;
                
                updatesAux(i,j,k,l,:) = [updatesLin;updatesNonLin];
                
                updates(i,j,k,l,:) = [updatesLin - updatesNonLin;updatesNonLin];
                
                
                tueLin = ((updatesLin - updatesNonLin) * linSamples * numberOfIterations)/100;
                tueNonLin = (updatesNonLin * (nonLinSamples + linSamples) * numberOfIterations)/100;

                TUE(k,i) = tueNonLin  + tueLin;

%                 TUEAux(k,i) = 2*(updatesLin - updatesNonLin)/(upVector(i-3)-updatesNonLin) - 4;

                
            end
            
            
            H = legend('$\bar{\gamma}_2 = \bar{\gamma}_1$','$\bar{\gamma}_2 = 1.5\bar{\gamma}_1$','$\bar{\gamma}_2 = 2\bar{\gamma}_1$',...
                '$\bar{\gamma}_2 = 2.5\bar{\gamma}_1$','$\bar{\gamma}_2 = 3\bar{\gamma}_1$','$\bar{\gamma}_2 = 3.5\bar{\gamma}_1$',...
                '$\bar{\gamma}_2 = 4\bar{\gamma}_1$');
            set(H,'interpreter','latex')
            ylim([-15 10]);
            
            ylabel('MSE [dB]','interpreter','latex');
            xlabel('Iterations [$k$]','interpreter','latex');
            xlim([0 numberOfIterations])
%             formatFig( gcf ,['.' filesep 'figs' filesep 'mseEqDT_' num2str(i)],'en' , figProp );
            
            
        end
%         figure;
%         bar(1:0.5:4,TUE(:,i))
%         TUE_s = upVector(i-3)*(nonLinSamples + linSamples)*numberOfIterations/100;
%         ylabel('NUE$_{\mathrm{DT}}$','interpreter','latex');
%         xlabel('$\times\bar{\gamma}_1$','interpreter','latex');
%         xlim([0.8 4.2])
%         hold on;
%         plot([0.8 4.2],[TUE_s TUE_s],'r','linewidth',2)
%         xticks(1:0.5:4)
%         xticklabels({'1','1.5','2','2.5','3','3.5','4'})
        
        
    end
    
%     figure;
%     bar(1:0.5:4,TUE(:,4:5))
%     TUE_s = upVector*(nonLinSamples + linSamples)*numberOfIterations/100;
%     ylabel('NUE$_{\mathrm{DT}}$','interpreter','latex');
%     xlabel('$\times\bar{\gamma}_1$','interpreter','latex');
%     xlim([0.8 4.2])
%     hold on;
%     plot([0.8 4.2],[TUE_s(1) TUE_s(1)],'color',[0    0.4470    0.7410],'linewidth',2)
%     plot([0.8 4.2],[TUE_s(2) TUE_s(2)],'color','r','linewidth',2)
%     xticks(1:0.5:4)
%     xticklabels({'1','1.5','2','2.5','3','3.5','4'})
%     H = legend('$M_{\mathrm{FF}} = 3$','$M_{\mathrm{FF}} = 4$');
%     set(H,'interpreter','latex')
%     set(H,'location','best');
%     formatFig( gcf ,['.' filesep 'figs' filesep 'TUEEqDT'],'en' , figProp );
end
  



fileVector = {'results_SMDTPU_04'};

MFF = 6:8;
MFB = 0;

for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    for i = 3:size(e4,1)
        for j = 1:1%size(e4,1)
            linSamples = MFF(i) + 1 + MFB(j) + 1;
            nonLinSamples = nchoosek(MFF(i)+1,N) + MFF(i) + 1;
            figure
            
            for k = 1:size(e4{i,j},2)
            
                x = e4{i,j}{1,k};
                aux = find(x,1);
                xAux = 10*log10(x(aux:end));
                plot(xAux)
                
                meanMSE = mean(10*log10((e4{i,j}{j,k}(end - 999:end))));
                stdMse = std(10*log10((e4{i,j}{j,k}(end - 999:end))));
                convergenceSample(i,j,k) = find(10*log10(e4{i,j}{j,k}(aux:end)) < meanMSE + 2*stdMse*sign(2*stdMse),1,'first');
                
                hold on
                updatesLin = mean(meanCountLin2{i,j}{1,k}(aux:end))*100;  
                updatesNonLin = mean(meanCountNonLin2{i,j}{1,k}(aux:end))*100; 
                updatesAux2(i,j,k,l,:) = [updatesLin;updatesNonLin];
                updates(i,j,k,l,:) = [updatesLin - updatesNonLin;updatesNonLin];
                
                
                tueLin = ((updatesLin - updatesNonLin) * linSamples * numberOfIterations)/100;
                tueNonLin = (updatesNonLin * (nonLinSamples + linSamples) * numberOfIterations)/100;

                TUE(k,i,j) = tueNonLin  + tueLin;

%                 TUEAux(k,i,j) = 2*(updatesLin - updatesNonLin)/(upVector-updatesNonLin) - 4;

            end
             
            H = legend('$\bar{\gamma}_2 = \bar{\gamma}_1$','$\bar{\gamma}_2 = 1.5\bar{\gamma}_1$','$\bar{\gamma}_2 = 2\bar{\gamma}_1$',...
                '$\bar{\gamma}_2 = 2.5\bar{\gamma}_1$','$\bar{\gamma}_2 = 3\bar{\gamma}_1$','$\bar{\gamma}_2 = 3.5\bar{\gamma}_1$',...
                '$\bar{\gamma}_2 = 4\bar{\gamma}_1$');
            set(H,'interpreter','latex')
            ylim([-15 10]);
            
            ylabel('MSE [dB]','interpreter','latex');
            xlabel('Iterations [$k$]','interpreter','latex');
            xlim([0 numberOfIterations])
%             formatFig( gcf ,['.' filesep 'figs' filesep 'mseEqDFEDT_' num2str(i)],'en' , figProp );
            
        end
        
%         xlim([0 3000])
    end
    
    
%     figure;
%     bar(1:0.5:4,TUE(:,5,1))
%     TUE_s = upVector*(nonLinSamples + linSamples)*numberOfIterations/100;
%     ylabel('NUE$_{\mathrm{DT}}$','interpreter','latex');
%     xlabel('$\times\bar{\gamma}_1$','interpreter','latex');
%     xlim([0.8 4.2])
%     hold on;
%     plot([0.8 4.2],[TUE_s(1) TUE_s(1)],'color',[0.8500 0.3250 0.0980],'linewidth',2)
%     xticks(1:0.5:4)
%     xticklabels({'1','1.5','2','2.5','3','3.5','4'})
%     H = legend('$M_{\mathrm{FF}} = 3$','$M_{\mathrm{FF}} = 4$');
%     set(H,'interpreter','latex')
%     set(H,'location','best');
%     formatFig( gcf ,['.' filesep 'figs' filesep 'TUEEqDFEDT'],'en' , figProp );
  
end




rmpath(['.' filesep 'results']);
rmpath(['..' filesep '..' filesep 'Misc' filesep]); 
