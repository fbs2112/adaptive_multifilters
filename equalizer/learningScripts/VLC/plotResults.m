clear;
close all;
clc;

addpath(['.' filesep 'results']);


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

fileVector = {'resultsTest'};

for l = 1:length(fileVector)
    figure
    load([fileVector{l} '.mat']);
    for m = 1:size(e4,2)
        for i = 1:size(e4,1)
            for j = 1:size(e4{i,m},1)
                
                for k = 1:size(e4{i,m},2)
                    
                    x = e4{i,m}{j,k};
                    aux = find(x,1);
                    xAux = 10*log10(x(aux:end));
                    plot(xAux)
                    
%                     updatesLin = mean(meanCountLin2{i,m}{j,k}(aux:4999))*100;
%                     updatesNonLin = mean(meanCountNonLin2{i,m}{j,k}(aux:4999))*100;
                    
                   
                    
                    hold on
                    
%                     updatesAux(i,j,k,l,m,:) = [updatesLin updatesNonLin];
                    updatesAux(i,j,k,l,m,:) = mean(meanCount2{i,m}{j,k}(aux:4999))*100;
                end
                
                
                H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
                    '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
                    '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
                set(H,'interpreter','latex')
%                 ylim([-15 10]);
                
                ylabel('MSE [dB]','interpreter','latex');
                xlabel('Iterations [$k$]','interpreter','latex');
                xlim([0 3000])
                
                
            end            
        end
       
    end
    
end
    



% clear
% 
% blindIt = [2000 2000 2000];
% 
% fileVector = {'resultsTestDFE_FF_FB' 'resultsTestDFE_FF' 'resultsTestDFE_FB'};
% 
% for l = 1:length(fileVector)
%     
%     load([fileVector{l} '.mat']);
%     figure
%     for i = 5:size(e4,2)
%         for j = 5:size(e4,1)
%             for k = 1:size(e4{i,j},2)
%             
%                 x = e4{i,j}{1,k};
%                 aux = find(x,1);
%                 xAux = 10*log10(x(aux:end));
%                 plot(xAux)
%                 
%                 hold on
%                 updates = mean(meanCount2{i,j}{1,k}(aux:4999 + aux))*100;
%                 updatesSup = mean(meanCount2{i,j}{1,k}(aux:blindIt(l)))*100;
%                 updatesBlind = mean(meanCount2{i,j}{1,k}(blindIt(l):4999 + aux))*100; 
%                 updatesAux2(i,j,k,l,:) = [updates;updatesSup;updatesBlind];
%             end
%                 
%             
%         end
%          H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
%             '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
%             '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
%         set(H,'interpreter','latex')
%         ylim([-15 10]);
%         
%         ylabel('MSE [dB]','interpreter','latex');
%         xlabel('Iterations [$k$]','interpreter','latex');
% %         xlim([0 3000])
%     end
%   
% end
% 
% 
% 

clear
close all;

fileVector = {'resultsTestDFE_FF_FB'};

for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    for i = 1:size(e4,1)
        figure
    for m = 1:size(e4,3)
        
            for j = 1:1%size(e4,1)
                for k = 1:size(e4{i,j,m},2)
                    
                    x = e4{i,j,m}{1,k};
                    aux = find(x,1);
                    xAux = 10*log10(x(aux:end));
                    plot(xAux)
                    hold on
%                     updatesLin = mean(meanCountLin2{i,j,m}{1,k}(aux:4999))*100;
%                     updatesNonLin = mean(meanCountNonLin2{i,j,m}{1,k}(aux:4999))*100;
%                     
%                     updatesAux(i,j,k,l,m,:) = [updatesLin updatesNonLin];
                    updatesAux(i,j,k,l,m,:) =  mean(meanCount2{i,j,m}{1,k}(aux:4999))*100;
                    
                end
                
                
            end
            H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
                '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
                '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
            set(H,'interpreter','latex')
            ylim([-15 10]);
            
            ylabel('MSE [dB]','interpreter','latex');
            xlabel('Iterations [$k$]','interpreter','latex');
            %         xlim([0 3000])
        end
    end
    
end



% clear
% close all;
% 
% fileVector = {'resultsTestDFE_FF_FB'};
% 
% for l = 1:length(fileVector)
%     
%     load([fileVector{l} '.mat']);
%     figure
%     for m = 1:size(e4,3)
%         for i = 1:size(e4,2)
%             for j = 1:size(e4,1)
%                 for k = 1:size(e4{i,j,m},2)
%                     
%                     x = e4{i,j,m}{1,k};
%                     aux = find(x,1);
%                     xAux = 10*log10(x(aux:end));
%                     plot(xAux)
%                     
%                     blindItAux = blindIt(:,k,i,j,m);
%                     blindItAux2 = blindIt2(:,k,i,j,m);
%                     if m == 1
%                         updatesSup = mean(meanCount2{i,j,m}{1,k}(aux:4999))*100;
%                         updatesBlind = mean(meanCount2{i,j,m}{1,k}(aux:4999))*100;
%                         
%                     else
%                         meanBlindIt(l,m,i,j,k) = round(mean(blindItAux(blindItAux~=0)));
%                         meanBlindIt2(l,m,i,j,k) = round(mean(blindItAux2(blindItAux2~=0)));
%                         updatesSup = mean(meanCount2{i,j,m}{1,k}(aux:meanBlindIt(l,m,i,j,k)))*100;
%                         updatesBlind = mean(meanCount2{i,j,m}{1,k}(meanBlindIt(l,m,i,j,k):4999))*100;
%                     end
%                     
%                     hold on
%                     updates = mean(meanCount2{i,j,m}{1,k}(aux:4999))*100;
%                     
%                     updatesAux(i,j,k,l,m,:) = [updates;updatesSup;updatesBlind];
%                 end
%                 
%                 
%             end
%             H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
%                 '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
%                 '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
%             set(H,'interpreter','latex')
%             ylim([-15 10]);
%             
%             ylabel('MSE [dB]','interpreter','latex');
%             xlabel('Iterations [$k$]','interpreter','latex');
%             %         xlim([0 3000])
%         end
%     end
%     
% end





rmpath(['.' filesep 'results']);
