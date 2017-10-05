%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);

load paramEq.mat;

numberOfSymbols = 2^numberOfBits;

L = 0;%PAPA NLMS

e3 = cell(length(L));
meanCountLin = cell(length(L));
meanCountNonLin = cell(length(L));
meanCountTotal = cell(length(L));

maxIt = 100;

NIndex = 5;
delayVector = N(NIndex)+1;
delayVector2 = [N(NIndex)+1 N(NIndex)-2];

delay = 1;

barGammaLin = 4*sqrt(5*noisePower);
barGammaNonLin = barGammaLin*2;

muLambda = 0.2;
eta = 0.2;


for LIndex = 1:length(L)
    
    globalLength = maxRuns + N(5) + L(LIndex) - 1;
    countLin = zeros(globalLength,maxIt);
    countNonLin = zeros(globalLength,maxIt);
    countTotal = zeros(globalLength,maxIt);
    
    u = zeros(L(LIndex)+1,1);
    u(1) = 1;
    wIndex = zeros(adapFiltLength(N(5)),globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);
    
    for index = 1:maxIt
        index
        
        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        eLin = zeros(globalLength,1);
        eNonLin = zeros(globalLength,1);
        yLin = zeros(globalLength,1);
        yNonLin = zeros(globalLength,1);
        lambda = zeros(globalLength,1);
        p = zeros(globalLength,1);
        
        muLin = zeros(globalLength,1);
        muNonLin = zeros(globalLength,1);
        G = zeros(adapFiltLength(N(5)),adapFiltLength(N(5)),globalLength);
        
        xFiltered = zeros(globalLength,1);
        x = zeros(N(5),globalLength);
        input = randi([0,numberOfSymbols-1],globalLength,1);
        
        pilot = pammod(input,pamOrder,0,'gray');
        
        pilot2 = pilot.*sqrt(signalPower/var(pilot));
        
        xAux = zeros(length(pilot),size(h,2));
        
        for channelIndex = 1:size(h,2)
            aux2 = zeros(length(l1Pilot),1);
            xAux2 = zeros(length(pilot),1);
            
            for i = memoryChannelLength:length(pilot) %Channel 1
                xPilot = (pilot2(i:-1:i-memoryChannelLength+1));
                for lIndex = 1:length(l1Pilot)
                    aux2(lIndex,1) = xPilot(l1Pilot(lIndex),1)*(xPilot(l2Pilot(lIndex),1));
                end
                xConc = [xPilot;(aux2)];
                xAux2(i,1) = xConc.'*h(:,channelIndex);
            end
                        
            n = randn(globalLength,1);
            powerSignal = xAux2'*xAux2./(globalLength);
            powerNoiseAux = n'*n/(globalLength);
            powerNoise = (powerSignal/SNR);
            n = n.*sqrt(powerNoise/powerNoiseAux);
            
            xAux(:,channelIndex) = xAux2 + n;
            
        end
        
        
        xFlip = xAux(:,1);
        w = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;
        
        channelIndex = 1;
        
        for k = (adapFiltLength(NIndex) + max(delayVector2)):globalLength
            
            x(:,k) = xAux(k:-1:k-N(NIndex)+1,channelIndex);
            
            xTDLAux = zeros(length(l1{NIndex}),1);
            
            for lIndex = 1:length(l1{NIndex})
                xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex),k)*(x(l2{NIndex}(lIndex),k));
            end
            
            xAP = [x(:,k);xTDLAux];
            
            d(k) = (pilot(-delayVector(delay) + k + 1));
            
            
            yLin(k) = w(1:N(5),k)'*xAP(1:N(5),1);
            yNonLin(k) = w(N(5)+1:end,k)'*xAP(N(5)+1:end,1);
            
            eLin(k) = d(k) - yLin(k);
            
            eNonLin(k) = d(k) - yNonLin(k);
            
            
            yFinal =  yNonLin(k) + lambda(k)*(yLin(k) - yNonLin(k));
            
            e(k) = d(k) - yFinal;
            
            lambda(k+1) = lambda(k) + muLambda*e(k)*(yLin(k) - yNonLin(k))/(p(k) + 1e-8);
            p(k+1) = eta*p(k) + (1-eta)*(yLin(k) - yNonLin(k))^2;
           
            absoluteValueErrorLin = abs(eLin(k));
            absoluteValueErrorNonLin = abs(eNonLin(k));
            
            
            if absoluteValueErrorLin > barGammaLin
                muLin(k) = 1 - barGammaLin/absoluteValueErrorLin;
                G(1:N(5),1:N(5),k) = diag(((1 - kappa*muLin(k))/N(5)) + (kappa*muLin(k)*abs(w(1:N(5),k))/norm(w(1:N(5),k),1)));
                w(1:N(5),k+1) = w(1:N(5),k) + muLin(k)*G(1:N(5),1:N(5),k)*xAP(1:N(5),1)*((xAP(1:N(5),1)'*G(1:N(5),1:N(5),k)*xAP(1:N(5),1)+gamma*eye(L+1))\eye(L+1))*conj(eLin(k))*u;
                countLin(k,index) = 1;
            else
                w(1:N(5),k+1) = w(1:N(5),k);
            end
            
            
            if absoluteValueErrorNonLin > barGammaNonLin
                muNonLin(k) = 1 - barGammaNonLin/absoluteValueErrorNonLin;
                G(N(5)+1:end,N(5)+1:end,k) = diag(((1 - kappa*muNonLin(k))/(adapFiltLength(5)-N(5))) + (kappa*muNonLin(k)*abs(w(N(5)+1:end,k))/norm(w(N(5)+1:end,k),1)));
                w(N(5)+1:end,k+1) = w(N(5)+1:end,k) + muNonLin(k)*G(N(5)+1:end,N(5)+1:end,k)*...
                    xAP(N(5)+1:end,1)*((xAP(N(5):end,1)'*G(N(5):end,N(5):end,k)*xAP(N(5):end,1)+gamma*eye(L+1))\eye(L+1))*conj(eNonLin(k))*u;
                countNonLin(k,index) = 1;
            else
                w(N(5)+1:end,k+1) = w(N(5)+1:end,k);
            end
            
            if countLin(k,index) && countNonLin(k,index)
                countTotal(k,index) = 1;
            else
                countTotal(k,index) = 0;
            end
            
            
        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
    end
    
    meanCountLin{LIndex} = mean(countLin,2);
    meanCountNonLin{LIndex} = mean(countNonLin,2);
    meanCountTotal{LIndex} = mean(countTotal,2);
    
    w3 = mean(wIndex,3);
    e3{LIndex} = mean(e2,2);
    
end
save(['.' filesep 'results' filesep 'testeWam.mat'],'w3','e3','meanCountLin','meanCountNonLin','meanCountTotal');


rmpath(['..' filesep 'simParameters' filesep]);