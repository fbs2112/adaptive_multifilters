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

fileVector = {'results4', 'results5'};


for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
    for i = 1:size(e3,2)
        figure
        for j = 1:size(e3,1)
            
                x = e3{j,i};
                aux = find(x,1);
                xAux = 10*log10(x(aux:4999 + aux));
                plot(xAux)
                
                hold on
                updatesLin = mean(meanCountLin{j,i}(aux:4999 + aux))*100;
                updatesNonLin = mean(meanCountNonLin{j,i}(aux:4999 + aux))*100;
                updatesTotal = mean(meanCountTotal{j,i}(aux:4999 + aux))*100;
                
                updatesAux(l,i,j,:) = [updatesLin;updatesNonLin];
                
            
        end
        H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$','$L = 4$');
        set(H,'interpreter','latex')
        ylim([-15 10]);
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations [$k$]','interpreter','latex');
        xlim([0 3000])
    end
  
end




rmpath(['.' filesep 'results']);
rmpath(['..' filesep '..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);

