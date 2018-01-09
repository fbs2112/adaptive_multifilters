clear;
clc;
close all;



addpath(['.' filesep 'results']);
addpath(['..' filesep '..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);



linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  

% fileVector = {'teste6' 'teste5' 'teste4' 'teste' 'teste2' 'teste3'};

% fileVector = {'teste6' 'teste5' 'teste4' 'teste' 'teste2'};

fileVector = {'results4' 'results5'};


for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    figure;
    for i = 1:size(e3,2)
        for j = 1:size(e3,1)
            
                x = e3{j,i};
                aux = find(x,1);
                xAux = 10*log10(x(aux:4999 + aux));
                plot(xAux)
                
                hold on
                updatesLin = mean(meanCountLin{j,i}(aux:4999 + aux))*100;
                updatesNonLin = mean(meanCountNonLin{j,i}(aux:4999 + aux))*100;
                updatesTotal = mean(meanCountTotal{j,i}(aux:4999 + aux))*100;
                
                updatesAux(i,j,l,:) = [updatesLin;updatesNonLin];
                
            
        end
        H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
            '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
            '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
        set(H,'interpreter','latex')
        ylim([-15 10]);
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations [$k$]','interpreter','latex');
        xlim([0 3000])
    end
  
end

close all
clear;
fileVector = {'resultsTest2'};


for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    
    for i = 5:size(e4,1)
        for j = 1:size(e4{i,1},1)
            figure;
            
            for k = 1:size(e4{i,1},2)
                
                x = e4{i,1}{j,k};
                aux = find(x,1);
                xAux = 10*log10(x(aux:4999 + aux));
                plot(xAux)
                
                hold on
                updatesLin = mean(meanCountLin2{i,1}{j,k}(aux:4999 + aux))*100;
                updatesNonLin = mean(meanCountNonLin2{i,1}{j,k}(aux:4999 + aux))*100;
                
                updatesAux(i,j,k,l,:) = [updatesLin;updatesNonLin];
            end
            
            
            H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
                '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
                '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
            set(H,'interpreter','latex')
            ylim([-15 10]);
            
            ylabel('MSE [dB]','interpreter','latex');
            xlabel('Iterations [$k$]','interpreter','latex');
            xlim([0 3000])
            
            
        end
        
        
    end
    
end
  




clear
fileVector = {'resultsTestDFE_FF_FB' 'resultsTestDFE_FF'};

for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    figure
    for i = 5:size(e4,2)
        for j = 5:size(e4,1)
            for k = 1:size(e4{i,j},2)
            
                x = e4{i,j}{1,k};
                aux = find(x,1);
                xAux = 10*log10(x(aux:4999 + aux));
                plot(xAux)
                
                hold on
                updatesLin = mean(meanCountLin2{i,j}{1,k}(aux:4999 + aux))*100;  
                updatesNonLin = mean(meanCountNonLin2{i,j}{1,k}(aux:4999 + aux))*100; 
                updatesAux2(i,j,k,l,:) = [updatesLin;updatesNonLin];
            end
                
            
        end
         H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
            '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
            '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
        set(H,'interpreter','latex')
        ylim([-15 10]);
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations [$k$]','interpreter','latex');
        xlim([0 3000])
    end
  
end


fileVector = {'results7' 'results8'};

figure;

for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    for i = 5:size(e4,2)
        for j = 5:size(e4,1)
            for k = 1:size(e4{i,j},2)
            
                x = e4{i,j}{1}(:,k);
                aux = find(x,1);
                xAux = 10*log10(x(aux:4999 + aux));
                plot(xAux)
                
                hold on
                updatesLin = mean(meanCountLin2{i,j}{1}(aux:4999 + aux))*100;  
                updatesNonLin = mean(meanCountLin2{i,j}{1}(aux:4999 + aux))*100; 
                updatesAux2(i,j,k,l,:) = [updatesLin;updatesNonLin];
            end
                
            
        end
        H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$','$L = 4$');
        set(H,'interpreter','latex')
        ylim([-15 10]);
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations [$k$]','interpreter','latex');
        xlim([0 3000])
    end
  
end



clear;

fileVector = {'results42' 'results43'};

for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    for i = 5:size(e4,2)
        for j = 5:size(e4,1)
           
            
            x = e4{i,j}{1,1};
            aux = find(x,1);
            xAux = 10*log10(x(aux:4999 + aux));
            plot(xAux)
            
            hold on
%             updatesLin = mean(meanCountLin2{i,j}{1}(aux:4999 + aux))*100;
%             updatesNonLin = mean(meanCountLin2{i,j}{1}(aux:4999 + aux))*100;
%             updatesAux2(i,j,k,l,:) = [updatesLin;updatesNonLin];
                
            
        end
%         H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}',...
%             '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}',...
%             '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}');
%         set(H,'interpreter','latex')
        ylim([-15 10]);
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations [$k$]','interpreter','latex');
        xlim([0 3000])
    end
  
end

fileVector = {'results9' 'results10'};

figure;

for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    for i = 5:size(e4,2)
        for j = 5:size(e4,1)
            for k = 1:size(e4{i,j},2)
            
                x = e4{i,j}{1}(:,k);
                aux = find(x,1);
                xAux = 10*log10(x(aux:4999 + aux));
                plot(xAux)
                
                hold on
                updatesLin = mean(meanCountLin2{i,j}{1}(aux:4999 + aux))*100;  
                updatesNonLin = mean(meanCountLin2{i,j}{1}(aux:4999 + aux))*100; 
                updatesAux2(i,j,k,l,:) = [updatesLin;updatesNonLin];
            end
                
            
        end
        H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$','$L = 4$');
        set(H,'interpreter','latex')
        ylim([-15 10]);
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations [$k$]','interpreter','latex');
        xlim([0 3000])
    end
  
end



close all;
fileVector = {'results11' 'results12'};

for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    figure
    for i = 5:size(e4,2)
        for j = 5:size(e4,1)
            for k = 1:size(e4{i,j},1)
            
                x = e4{i,j}{k,1};
                aux = find(x,1);
                xAux = 10*log10(x(aux:4999 + aux));
                plot(xAux)
                
                hold on
                updatesLin = mean(meanCountLin2{i,j}{k,1}(aux:4999 + aux))*100;  
                updatesNonLin = mean(meanCountLin2{i,j}{k,1}(aux:4999 + aux))*100; 
                updatesAux2(i,j,k,l,:) = [updatesLin;updatesNonLin];
            end
                
            
        end
         H = legend('$\gamma_{\mathrm{NL}} = \gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 1.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 2\gamma_{\mathrm{L}}$',...
            '$\gamma_{\mathrm{NL}} = 2.5\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3\gamma_{\mathrm{L}}$','$\gamma_{\mathrm{NL}} = 3.5\gamma_{\mathrm{L}}$',...
            '$\gamma_{\mathrm{NL}} = 4\gamma_{\mathrm{L}}$');
        set(H,'interpreter','latex')
        ylim([-15 10]);
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations [$k$]','interpreter','latex');
        xlim([0 3000])
    end
  
end





rmpath(['.' filesep 'results']);
rmpath(['..' filesep '..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);

