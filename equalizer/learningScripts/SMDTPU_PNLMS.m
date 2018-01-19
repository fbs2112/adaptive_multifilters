%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);

load paramEq_DT.mat;

numberOfSymbols = 2^numberOfBits;

barGammaLin = 4*sqrt(5*noisePower);
barGammaNonLin = (1:0.5:4)*barGammaLin;

e4 = cell(length(N),1);
w4 = cell(length(N),1);
meanCountLin2 = cell(length(N),1);
meanCountNonLin2 = cell(length(N),1);

for NIndex = 3:length(N)
    
    CLin = diag([ones(N(NIndex),1).' zeros(adapFiltLength(N(NIndex)) - N(NIndex),1).'].');
    CNonLin = diag([zeros(N(NIndex),1).' ones(adapFiltLength(N(NIndex)) - N(NIndex),1).'].');
    
    delayVector = N(NIndex)+1;%adapFiltLength + 10;
    delayVector2 = [N(NIndex)+1 N(NIndex)-2];
    
    e3 = cell(length(delayVector),length(barGammaNonLin));
    w3 = cell(length(delayVector),length(barGammaNonLin));
    meanCountLin = cell(length(delayVector),length(barGammaNonLin));
    meanCountNonLin = cell(length(delayVector),length(barGammaNonLin));
    meanCountTotal = cell(length(delayVector),length(barGammaNonLin));

    for barGammaNonLinIndex = 1:length(barGammaNonLin)
        barGammaNonLinIndex


         for delay = 1:length(delayVector)

            globalLength = maxRuns + N(NIndex) - 1;
            countLin = zeros(globalLength,maxIt);
            countNonLin = zeros(globalLength,maxIt);
            countTotal = zeros(globalLength,maxIt);
            wIndex = zeros(adapFiltLength(N(NIndex)),globalLength,maxIt);
            e2 = zeros(globalLength,maxIt);

            for index = 1:maxIt
                index
                d = zeros(globalLength,1);
                e = zeros(globalLength,1);
                muLin = zeros(globalLength,1);
                muNonLin = zeros(globalLength,1);
                G = zeros(adapFiltLength(N(NIndex)),adapFiltLength(N(NIndex)),globalLength);

                G1 = zeros(adapFiltLength(N(NIndex)),adapFiltLength(N(NIndex)),globalLength);

                G2 = zeros(adapFiltLength(N(NIndex)),adapFiltLength(N(NIndex)),globalLength);

                xFiltered = zeros(globalLength,1);
                xLin = zeros(N(NIndex),globalLength);
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

                w1 = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;
                w2 = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;


                channelIndex = 1;

                for k = (adapFiltLength(NIndex) + max(delayVector2)):globalLength
                    
                    if k >= changingIteration
                        if N(NIndex) > 2
                            delayVector = delayVector2(2);
                        else
                            delayVector = delayVector2(2) + 2;
                        end
                        channelIndex = 2;
                    else
                        delayVector = delayVector2(1);
                        channelIndex = 1;
                    end

                    xLin(:,k) = xAux(k:-1:k-N(NIndex)+1,channelIndex);


                    xAP = zeros(adapFiltLength(NIndex),1);

                    for windowIndex = 1:1
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
                        G1(:,:,k) = diag(((1 - kappa*muLin(k))/adapFiltLength(N(NIndex))) + (kappa*muLin(k)*abs(w1(:,k))/norm(w1(:,k),1)));

                        w1(:,k+1) = w1(:,k) + muLin(k) * CLin*G1(:,:,k)*xAP*((xAP'*CLin*G1(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                        countLin(k,index) = 1;
                    else
                        w1(:,k+1) = w1(:,k);
                    end


                    if absoluteValueError > barGammaNonLin(barGammaNonLinIndex)
                        muNonLin(k) = 1 - barGammaNonLin(barGammaNonLinIndex)/absoluteValueError;

                        G2(:,:,k) = diag(((1 - kappa*muNonLin(k))/adapFiltLength(N(NIndex))) + (kappa*muNonLin(k)*abs(w2(:,k))/norm(w2(:,k),1)));

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

            meanCountLin{delay,barGammaNonLinIndex} = mean(countLin,2);
            meanCountNonLin{delay,barGammaNonLinIndex} = mean(countNonLin,2);

            w3{delay,barGammaNonLinIndex} = mean(wIndex,3);
            e3{delay,barGammaNonLinIndex} = mean(e2,2);

        end
    end
    meanCountLin2{NIndex} = meanCountLin;
    meanCountNonLin2{NIndex} = meanCountNonLin;
    w4{NIndex} = w3;
    e4{NIndex} = e3;
end

save(['.' filesep 'results' filesep 'results_SMDTPU_01.mat'],'w4','e4','meanCountLin2','meanCountNonLin2');

rmpath(['..' filesep 'simParameters' filesep]);
