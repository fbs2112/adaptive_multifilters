%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['..' filesep '..' filesep 'simParameters' filesep]);

load paramEq.mat;

numberOfSymbols = 2^numberOfBits;

L = 0;%PAPA NLMS

maxIt = 100;

NIndex = 5;
delayVector = N(NIndex)+1;
delayVector2 = [N(NIndex)+1 N(NIndex)-2];

delay = 1;

barGammaLin = 4*sqrt(5*noisePower);

barGammaNonLin = (1:0.5:4)*barGammaLin;
% barGammaNonLin = 1*barGammaLin;

gamma = 1e-8;

e3 = cell(length(L),length(barGammaNonLin));
w3 = cell(length(L),length(barGammaNonLin));
meanCountLin = cell(length(L),length(barGammaNonLin));
meanCountNonLin = cell(length(L),length(barGammaNonLin));
meanCountTotal = cell(length(L),length(barGammaNonLin));


for barGammaNonLinIndex = 1:length(barGammaNonLin)
    barGammaNonLinIndex
    

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
            muLin = zeros(globalLength,1);
            muNonLin = zeros(globalLength,1);
            G = zeros(adapFiltLength(N(5)),adapFiltLength(N(5)),globalLength);
            
            G1 = zeros(N(5),adapFiltLength(N(5)),globalLength);
            
            G2 = zeros(adapFiltLength(N(5)) - N(5),adapFiltLength(N(5)),globalLength);

            xFiltered = zeros(globalLength,1);
            xLin = zeros(N(5),globalLength);
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

            w = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;
            
            w1 = zeros(N(5),globalLength) + 1e-6;
            w2 = zeros(adapFiltLength(NIndex) - N(5),globalLength) + 1e-6;


            channelIndex = 1;

            for k = (adapFiltLength(NIndex) + max(delayVector2)):globalLength

                xLin(:,k) = xAux(k:-1:k-N(NIndex)+1,channelIndex);


                xAP = zeros(adapFiltLength(NIndex),L(LIndex)+1);

                for windowIndex = 1:L(LIndex)+1
                    xNonLin = zeros(length(l1{NIndex}),1);

                    for lIndex = 1:length(l1{NIndex})
                        xNonLin(lIndex,1) = xLin(l1{NIndex}(lIndex),k)*(xLin(l2{NIndex}(lIndex),k));
                    end
                    xAP(:,windowIndex) = [xLin(:,k);xNonLin];

                end

                xAP = fliplr(xAP);

                d(k) = (pilot(-delayVector(delay) + k + 1));


                e(k) = d(k) - w(:,k)'*xAP(:,1);
                absoluteValueError = abs(e(k));
                
                if absoluteValueError > barGammaLin
                    muLin(k) = 1 - barGammaLin/absoluteValueError;
%                     G(1:N(5),1:N(5),k) = diag(((1 - kappa*muLin(k))/N(5)) + (kappa*muLin(k)*abs(w(1:N(5),k))/norm(w(1:N(5),k),1)));
                    G(1:N(5),1:N(5),k) = diag(((1 - kappa*muLin(k))/N(5)) + (kappa*muLin(k)*abs(w1(:,k))/norm(w1(:,k),1)));

                    w1(:,k+1) = w1(:,k) + muLin(k) * G(1:N(5),1:N(5),k)*xLin(:,k)*((xLin(:,k)'* G(1:N(5),1:N(5),k)*xLin(:,k)+gamma*eye(L(LIndex)+1))\eye(L(LIndex)+1))*conj(e(k))*u;
%                     w(1:N(5),k+1) = w(1:N(5),k) + muLin(k)*G(1:N(5),1:N(5),k)*xAP(1:N(5),:)*((xAP'*G(:,:,k)*xAP+gamma*eye(L(LIndex)+1))\eye(L(LIndex)+1))*conj(e(k))*u;
                    countLin(k,index) = 1;
                else
                    w1(:,k+1) = w1(:,k);
%                     w(1:N(5),k+1) = w(1:N(5),k);
                end


                if absoluteValueError > barGammaNonLin(barGammaNonLinIndex)
                    muNonLin(k) = 1 - barGammaNonLin(barGammaNonLinIndex)/absoluteValueError;

%                     G(N(5)+1:end,N(5)+1:end,k) = diag(((1 - kappa*muNonLin(k))/(adapFiltLength(5)-N(5))) + (kappa*muNonLin(k)*abs(w(N(5)+1:end,k))/norm(w(N(5)+1:end,k),1)));
                    G(N(5)+1:end,N(5)+1:end,k) = diag(((1 - kappa*muNonLin(k))/(adapFiltLength(5)-N(5))) + (kappa*muNonLin(k)*abs(w2(:,k))/norm(w2(:,k),1)));

%                     w(N(5)+1:end,k+1) = w(N(5)+1:end,k) + muNonLin(k)*G(N(5)+1:end,N(5)+1:end,k)*...
%                         xAP(N(5)+1:end,:)*((xAP'*G(:,:,k)*xAP+gamma*eye(L(LIndex)+1))\eye(L(LIndex)+1))*conj(e(k))*u;
                    w2(:,k+1) = w2(:,k) + muNonLin(k)*G(N(5)+1:end,N(5)+1:end,k)*xNonLin*((xNonLin'*G(N(5)+1:end,N(5)+1:end,k)*xNonLin+gamma*eye(L(LIndex)+1))\eye(L(LIndex)+1))*conj(e(k))*u;
                    countNonLin(k,index) = 1;
                else
                    w2(:,k+1) = w2(:,k);
%                     w(N(5)+1:end,k+1) = w(N(5)+1:end,k);
                end

                if countLin(k,index) && countNonLin(k,index)
                    countTotal(k,index) = 1;
                else
                    countTotal(k,index) = 0;
                end

                
                w(:,k+1) = [w1(:,k+1);w2(:,k+1)];

            end
            wIndex(:,:,index) = conj(w(:,1:globalLength));
            e2(:,index) = abs(e).^2;
        end

        meanCountLin{LIndex,barGammaNonLinIndex} = mean(countLin,2);
        meanCountNonLin{LIndex,barGammaNonLinIndex} = mean(countNonLin,2);
        meanCountTotal{LIndex,barGammaNonLinIndex} = mean(countTotal,2);

        w3{LIndex,barGammaNonLinIndex} = mean(wIndex,3);
        e3{LIndex,barGammaNonLinIndex} = mean(e2,2);

    end
end


save(['.' filesep 'results' filesep 'results5.mat'],'w3','e3','meanCountLin','meanCountNonLin','meanCountTotal');

rmpath(['..' filesep '..' filesep 'simParameters' filesep]);
