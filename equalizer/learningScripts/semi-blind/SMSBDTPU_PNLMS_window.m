%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['..' filesep '..' filesep 'simParameters' filesep]);

load paramEq.mat;

numberOfSymbols = 2^numberOfBits;

barGammaLin = 4*sqrt(5*noisePower);
barGammaNonLin = (1:0.5:4)*barGammaLin;

windowLength = 100;
eta = 0.1:0.1:0.3;
lambdaUp = 0.5;

e4 = cell(length(N),length(eta));
w4 = cell(length(N),length(eta));
meanCountLin2 = cell(length(N),length(eta));
meanCountNonLin2 = cell(length(N),length(eta));
blindIt = zeros(maxIt,1,length(N),length(eta),length(barGammaNonLin));
blindIt2 = zeros(maxIt,1,length(N),length(eta),length(barGammaNonLin));

for etaIndex = 1:1%length(eta)
    for NIndex = 3:length(N)
        
        CLin = diag([ones(N(NIndex),1).' zeros(adapFiltLength(NIndex) - N(NIndex),1).'].');
        CNonLin = diag([zeros(N(NIndex),1).' ones(adapFiltLength(NIndex) - N(NIndex),1).'].');
        
        delayVector = N(NIndex)+1;%adapFiltLength + 10;
        delayVector2 = [N(NIndex)+1 N(NIndex)-2];

        e3 = cell(length(delayVector),length(barGammaNonLin));
        w3 = cell(length(delayVector),length(barGammaNonLin));
        meanCountLin = cell(length(delayVector),length(barGammaNonLin));
        meanCountNonLin = cell(length(delayVector),length(barGammaNonLin));

        for barGammaNonLinIndex = 1:length(barGammaNonLin)
            barGammaNonLinIndex


             for delay = 1:length(delayVector)

                globalLength = maxRuns + N(NIndex) - 1;
                countLin = zeros(globalLength,maxIt);
                countNonLin = zeros(globalLength,maxIt);
                wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
                e2 = zeros(globalLength,maxIt);

                for index = 1:maxIt
                    index
                    d = zeros(globalLength,1);
                    e = zeros(globalLength,1);
                    muLin = zeros(globalLength,1);
                    muNonLin = zeros(globalLength,1);
                    gammaAux = zeros(globalLength,1);
                    medianAux = zeros(globalLength,1);
                    G1 = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);
                    G2 = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);
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
                    blindFlag = 0;

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
                        y = w(:,k)'*xAP(:,1);

                        if (k > (adapFiltLength(NIndex) + max(delayVector2)) + windowLength && k < changingIteration) || k > changingIteration + windowLength
                            medianAux(k) = median(abs(e(k-1 - windowLength:k-1)));
                            if medianAux(k) <= 2*eta(etaIndex) || blindFlag == 1
                                d(k) = pamHardThreshold(y);

                                if ~blindFlag && k < changingIteration
                                    blindIt(index,delay,NIndex,etaIndex,barGammaNonLinIndex) = k;
                                elseif ~blindFlag && k > changingIteration
                                    blindIt2(index,delay,NIndex,etaIndex,barGammaNonLinIndex) = k;
                                end
                                blindFlag = 1;

                            else
                                d(k) = (pilot(-delayVector(delay) + k + 1));
                            end

                        else
                            d(k) = (pilot(-delayVector(delay) + k + 1));
                        end

                        if k == changingIteration
                            blindFlag = 0;
                        end

                        e(k) = d(k) - y;
                        absoluteValueError = abs(e(k));

                        if absoluteValueError > barGammaLin
                            muLin(k) = 1 - barGammaLin/absoluteValueError;
                            G1(:,:,k) = diag(((1 - kappa*muLin(k))/adapFiltLength(NIndex)) + (kappa*muLin(k)*abs(w1(:,k))/norm(w1(:,k),1)));
                            
                            w1(:,k+1) = w1(:,k) + muLin(k) * CLin*G1(:,:,k)*xAP*((xAP'*CLin*G1(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                            countLin(k,index) = 1;
                        else
                            w1(:,k+1) = w1(:,k);
                        end
                        
                        
                        if absoluteValueError > barGammaNonLin(barGammaNonLinIndex)
                            muNonLin(k) = 1 - barGammaNonLin(barGammaNonLinIndex)/absoluteValueError;
                            
                            G2(:,:,k) = diag(((1 - kappa*muNonLin(k))/adapFiltLength(NIndex)) + (kappa*muNonLin(k)*abs(w2(:,k))/norm(w2(:,k),1)));
                            
                            w2(:,k+1) = w2(:,k) + muNonLin(k)*CNonLin*G2(:,:,k)*xAP*((xAP'*CNonLin*G2(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                            countNonLin(k,index) = 1;
                        else
                            w2(:,k+1) = w2(:,k);
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
        meanCountLin2{NIndex,etaIndex} = meanCountLin;
        meanCountNonLin2{NIndex,etaIndex} = meanCountNonLin;
        w4{NIndex,etaIndex} = w3;
        e4{NIndex,etaIndex} = e3;
    end
end

save(['.' filesep 'results' filesep 'SBSMDTPU_results03.mat'],'w4','e4','meanCountLin2','meanCountNonLin2','blindIt','blindIt2');

rmpath(['..' filesep '..' filesep 'simParameters' filesep]);
