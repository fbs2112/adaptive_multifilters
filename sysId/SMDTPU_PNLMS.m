%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['.' filesep 'simParameters' filesep]);

load param01.mat;

barGammaLin = sqrt(5*noisePower);
barGammaNonLin = (1:0.5:4)*barGammaLin;

gamma = 1e-1;
lambdaUp = 0.5;

CLin = diag([ones(N,1).' zeros(adapFiltLength - N,1).'].');
CNonLin = diag([zeros(N,1).' ones(adapFiltLength - N,1).'].');

e3 = cell(length(barGammaNonLin),1);
w3 = cell(length(barGammaNonLin),1);
meanCountLin = cell(length(barGammaNonLin),1);
meanCountNonLin = cell(length(barGammaNonLin),1);
meanCountTotal = cell(length(barGammaNonLin),1);

for barGammaNonLinIndex = 1:length(barGammaNonLin)
    barGammaNonLinIndex
    
    
    
    globalLength = maxRuns + N - 1;
    countLin = zeros(globalLength,maxIt);
    countNonLin = zeros(globalLength,maxIt);
    countTotal = zeros(globalLength,maxIt);
    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);
    
    for index = 1:maxIt
        index
        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        muLin = zeros(globalLength,1);
        muNonLin = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);
        
        G1 = zeros(adapFiltLength,adapFiltLength,globalLength);
        
        G2 = zeros(adapFiltLength,adapFiltLength,globalLength);
        
        xFiltered = zeros(globalLength,1);
        xLin = zeros(N,globalLength);
        input = randi([0,pamOrder-1],globalLength,1);
        
        modSignal = pammod(input,pamOrder,0,'gray');
        
        modSignal = modSignal.*sqrt(signalPower/var(modSignal));
        
        n = randn(globalLength,1);
        n = n.*sqrt(noisePower/var(n));
        
        w = zeros(adapFiltLength,globalLength) + 1e-6;
        
        w1 = zeros(adapFiltLength,globalLength) + 1e-6;
        w2 = zeros(adapFiltLength,globalLength) + 1e-6;
        
        woIndex = 1;
        
        for k = N:globalLength
            
            if k >= changingIteration
                woIndex = 2;
            end
            
            xLin(:,k) = modSignal(k:-1:k-N+1);
            
            
            xAP = zeros(adapFiltLength,1);
            
            for windowIndex = 1:1
                xNonLin = zeros(length(l1),1);
                
                for lIndex = 1:length(l1)
                    xNonLin(lIndex,1) = xLin(l1(lIndex),k)*(xLin(l2(lIndex),k));
                end
                xAP(:,windowIndex) = [xLin(:,k);xNonLin];
                
            end
            
            xAP = fliplr(xAP);
            
            d(k) = ((wo(:,woIndex)'*xAP))  + n(k);
            
            e(k) = d(k) - w(:,k)'*xAP(:,1);
            absoluteValueError = abs(e(k));
            
%             if pow2db(abs(e(k)^2 - e(k-1)^2)) > 10 && k > 200 && k < changingIteration
%                 e(k)
%             end
            
            if absoluteValueError > barGammaLin
                muLin(k) = 1 - barGammaLin/absoluteValueError;
                G1(:,:,k) = diag(((1 - kappa*muLin(k))/adapFiltLength) + (kappa*muLin(k)*abs(w1(:,k))/(norm(w1(:,k),1))));
                
                w1(:,k+1) = w1(:,k) + muLin(k) * CLin*G1(:,:,k)*xAP*((xAP'*CLin*G1(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                countLin(k,index) = 1;
            else
                w1(:,k+1) = w1(:,k);
            end
            
            
            if absoluteValueError > barGammaNonLin(barGammaNonLinIndex)
                muNonLin(k) = 1 - barGammaNonLin(barGammaNonLinIndex)/absoluteValueError;
                
                G2(:,:,k) = diag(((1 - kappa*muNonLin(k))/adapFiltLength) + (kappa*muNonLin(k)*abs(w2(:,k))/(norm(w2(:,k),1))));
                                
                w2(:,k+1) = w2(:,k) + muNonLin(k)*CNonLin*G2(:,:,k)*xAP*((xAP'*CNonLin*G2(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                countNonLin(k,index) = 1;
            else
                w2(:,k+1) = w2(:,k);
            end
            
            if countLin(k,index) && countNonLin(k,index)
                countTotal(k,index) = 1;
            else
                countTotal(k,index) = 0;
            end
            
            
            w(:,k+1) = (lambdaUp*w1(:,k+1)+ (1-lambdaUp)*w2(:,k+1));
            w1(:,k+1) = w(:,k+1);
            w2(:,k+1) = w(:,k+1);
            
        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
    end
    
    meanCountLin{barGammaNonLinIndex} = mean(countLin,2);
    meanCountNonLin{barGammaNonLinIndex} = mean(countNonLin,2);
    
    w3{barGammaNonLinIndex} = mean(wIndex,3);
    e3{barGammaNonLinIndex} = mean(e2,2);
    
end


save(['.' filesep 'results' filesep 'results_sysId01.mat'],'w3','e3','meanCountLin','meanCountNonLin');

rmpath(['.' filesep 'simParameters' filesep]);
