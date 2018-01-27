%Teste Volterra SM-NLMS

clear;
clc;
close all;


addpath(['..' filesep '..' filesep 'VLC_param' filesep]);

load paramEqVLC.mat;
load VLC_param01.mat;

numberOfSymbols = 2^numberOfBits;

modulationIndexVector = [0.05 0.075 0.1];

e4 = cell(length(N),length(modulationIndexVector));
w4 = cell(length(N),length(modulationIndexVector));
meanCount2 = cell(length(N),length(modulationIndexVector));

for modulationIndex = 1:length(modulationIndexVector)

    maxVoltage = VDC*(1+modulationIndexVector(modulationIndex));
    for NIndex = 1:length(N)
        
        delayVector = N(NIndex) + 1;

        e3 = cell(length(delayVector),1);
        w3 = cell(length(delayVector),1);
        meanCount = cell(length(delayVector));


        for delay = 1:length(delayVector)
            delay
            globalLength = maxRuns + adapFiltLength(NIndex) + delayVector - 1;

            wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
            e2 = zeros(globalLength,maxIt);

            count = zeros(globalLength,maxIt);

            w2 = zeros(adapFiltLength(NIndex),globalLength,maxIt);
            for index = 1:maxIt
                index

                d = zeros(globalLength,1);
                P = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);
                P(:,:,adapFiltLength(NIndex) + delayVector) = eye(adapFiltLength(NIndex))*1e-6;
                sigma = zeros(globalLength,1);
                sigma(adapFiltLength(NIndex) + delayVector(delay)) = 1;
                delta = zeros(globalLength,1);
                lambda = zeros(globalLength,1);
                G = zeros(globalLength,1);

                x = zeros(N(NIndex),globalLength);
                xAP = zeros(adapFiltLength(NIndex),globalLength);

                input = randi([0,numberOfSymbols-1],globalLength,1);

                pilot = pammod(input,pamOrder,0,'gray');

                pilot2 = pilot.*sqrt(signalPower/var(pilot));

                xAux = VLC_channel(pilot2, modulationIndexVector(modulationIndex), maxVoltage, VDC, SNR);

                xAux = [zeros(N(NIndex)-1,1);xAux];

                theta = zeros(adapFiltLength(NIndex),globalLength);

                channelIndex = 1;

                for k = (adapFiltLength(NIndex) + delayVector):globalLength

                    x(:,k) = xAux(k:-1:k-N(NIndex)+1,channelIndex);

                    xTDLAux = zeros(length(l1{NIndex}),1);

                    for lIndex = 1:length(l1{NIndex})
                        xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex),k)*(x(l2{NIndex}(lIndex),k));
                    end

                    xAP(:,k) = [x(:,k);xTDLAux];

                    d(k) = (pilot(-delayVector(delay) + k + 1)); 

                    delta(k) = d(k) - theta(:,k).'*xAP(:,k);

                    if abs(delta(k)) > barGamma 
                        G(k) = xAP(:,k).'*P(:,:,k)*conj(xAP(:,k));
                        lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);
                        lambda2 = 1/lambda(k);

                        P(:,:,k+1) = lambda(k)*(P(:,:,k) - (P(:,:,k)*conj(xAP(:,k))*xAP(:,k).'*P(:,:,k))/(lambda2+G(k)));

                        theta(:,k+1) = theta(:,k) + P(:,:,k+1)*conj(xAP(:,k))*delta(k);

                        count(k,index) = 1; 
                    else
                        lambda(k) = 0;
                        P(:,:,k+1) = P(:,:,k);
                        theta(:,k+1) = theta(:,k);
                    end

               end
               wIndex(:,:,index) = (theta(:,1:globalLength));
               e2(:,index) = abs(delta).^2;
            end

            meanCount{delay} = mean(count,2);

            w3{delay} = mean(wIndex,3);

            e3{delay} = mean(e2,2);
        end
        w4{NIndex,modulationIndex} = w3;
        e4{NIndex,modulationIndex} = e3;
        meanCount2{NIndex,modulationIndex} = meanCount;

    end
end


save(['.' filesep 'results' filesep 'results_BEACON_VLC_01.mat'],'e4','w4','meanCount2');

rmpath(['..' filesep '..' filesep 'VLC_param' filesep]);

