%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['.' filesep 'simParameters' filesep]);

load param01.mat;

maxIt = 100;

barGammaLin = sqrt(5*noisePower);
barGammaNonLin = (1:0.5:4)*barGammaLin;
barGammaNonLin = (1)*barGammaLin;

lambdaUp = 0.5;

CLin = diag([ones(N,1).' zeros(adapFiltLength - N,1).'].');

CLin = eye(adapFiltLength);
CNonLin = eye(adapFiltLength);
% CNonLin = diag([zeros(N,1).' ones(adapFiltLength - N,1).'].');

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
        delta = zeros(globalLength,1);
        muLin = zeros(globalLength,1);
        muNonLin = zeros(globalLength,1);
        
        G1 = zeros(globalLength,1);
        G2 = zeros(globalLength,1);
        
        P1 = zeros(adapFiltLength,adapFiltLength,globalLength);
        P1(:,:,N) = eye(adapFiltLength)*1e-6;
        
        P2 = zeros(adapFiltLength,adapFiltLength,globalLength);
        P2(:,:,N) = eye(adapFiltLength)*1e-6;
        
        xFiltered = zeros(globalLength,1);
        xLin = zeros(N,globalLength);
        input = randi([0,pamOrder-1],globalLength,1);
        
        modSignal = pammod(input,pamOrder,0,'gray');
        
        modSignal = modSignal.*sqrt(signalPower/var(modSignal));
        
        n = randn(globalLength,1);
        n = n.*sqrt(noisePower/var(n));
        
        theta = zeros(adapFiltLength,globalLength);
        
        theta1 = zeros(adapFiltLength,globalLength);
        theta2 = zeros(adapFiltLength,globalLength);
        
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
            
            delta(k) = d(k) - theta(:,k)'*xAP(:,1);
            absoluteValueError = abs(delta(k));
            
%             if pow2db(abs(e(k)^2 - e(k-1)^2)) > 10 && k > 200 && k < changingIteration
%                 e(k)
%             end
            
           
            if absoluteValueError > barGammaLin
                G1(k) = xAP.'*CLin*P1(:,:,k)*conj(xAP);
                lambda(k) = (1/G1(k))*((abs(delta(k))/barGammaLin) - 1);
                lambda2 = 1/lambda(k);
                
                P1(:,:,k+1) = lambda(k)*(CLin*P1(:,:,k) - (CLin*P1(:,:,k)*conj(xAP)*xAP.'*CLin*P1(:,:,k))/(lambda2+G1(k)));
                
                theta1(:,k+1) = theta1(:,k) + P1(:,:,k+1)*conj(xAP)*delta(k);
                                
                countLin(k,index) = 1;
            else
                theta1(:,k+1) = theta1(:,k);
                P1(:,:,k+1) = P1(:,:,k);
            end
            
            
            if absoluteValueError > barGammaNonLin(barGammaNonLinIndex) 
                G2(k) = xAP.'*CNonLin*P2(:,:,k)*conj(xAP);
                lambda(k) = (1/G2(k))*((abs(delta(k))/barGammaNonLin(barGammaNonLinIndex)) - 1);
                lambda2 = 1/lambda(k);
                
                P2(:,:,k+1) = lambda(k)*(CNonLin*P2(:,:,k) - (CNonLin*P2(:,:,k)*conj(xAP)*xAP.'*CNonLin*P2(:,:,k))/(lambda2+G2(k)));
                
                theta2(:,k+1) = theta2(:,k) + P2(:,:,k+1)*conj(xAP)*delta(k);
                countNonLin(k,index) = 1;
            else
                theta2(:,k+1) = theta2(:,k);
                P2(:,:,k+1) = P2(:,:,k);
            end
            
            if countLin(k,index) && countNonLin(k,index)
                countTotal(k,index) = 1;
            else
                countTotal(k,index) = 0;
            end
            
            
            theta(:,k+1) = 2*(lambdaUp*theta1(:,k+1)+ (1-lambdaUp)*theta2(:,k+1));
            if isnan(theta(:,k+1))
                pause;
            end
            
        end
        wIndex(:,:,index) = conj(theta(:,1:globalLength));
        e2(:,index) = abs(delta).^2;
    end
    
    meanCountLin{barGammaNonLinIndex} = mean(countLin,2);
    meanCountNonLin{barGammaNonLinIndex} = mean(countNonLin,2);
    
    w3{barGammaNonLinIndex} = mean(wIndex,3);
    e3{barGammaNonLinIndex} = mean(e2,2);
    
end


save(['.' filesep 'results' filesep 'resultsTest2.mat'],'w3','e3','meanCountLin','meanCountNonLin');

rmpath(['.' filesep 'simParameters' filesep]);
