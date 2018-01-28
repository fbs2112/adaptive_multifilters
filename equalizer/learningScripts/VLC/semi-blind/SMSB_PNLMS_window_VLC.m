%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['..' filesep '..' filesep '..' filesep 'VLC_param' filesep]); 
addpath(['..' filesep '..' filesep '..' filesep 'VLC_param' filesep 'Utils' filesep]); 
addpath(['..' filesep '..' filesep  '..' filesep 'VLC_param' filesep 'LED Parameters' filesep]); 

load paramEqVLC.mat;
load VLC_param01.mat;

numberOfSymbols = 2^numberOfBits;

barGammaNonLin = 1;
windowLength = 100;
eta = 0.1:0.1:0.3;

modulationIndexVector = [0.05 0.075 0.1];

e4 = cell(length(N),length(eta),length(modulationIndexVector));
w4 = cell(length(N),length(eta),length(modulationIndexVector));
meanCount2 = cell(length(N),length(eta),length(modulationIndexVector));
blindIt = zeros(maxIt,1,length(N),length(eta),length(modulationIndexVector));
blindIt2 = zeros(maxIt,1,length(N),length(eta),length(modulationIndexVector));

for modulationIndex = 1:length(modulationIndexVector)
    maxVoltage = VDC*(1+modulationIndexVector(modulationIndex));

    for etaIndex = 1:length(eta)
        for NIndex = 1:length(N)

            delayVector = N(NIndex)+1;%adapFiltLength + 10;

            e3 = cell(length(delayVector),length(barGammaNonLin));
            w3 = cell(length(delayVector),length(barGammaNonLin));
            meanCount = cell(length(delayVector),length(barGammaNonLin));

            for barGammaNonLinIndex = 1:length(barGammaNonLin)
                barGammaNonLinIndex


                 for delay = 1:length(delayVector)

                    globalLength = maxRuns + N(NIndex) - 1;
                    count = zeros(globalLength,maxIt);
                    wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
                    e2 = zeros(globalLength,maxIt);

                    for index = 1:maxIt
                        index
                        d = zeros(globalLength,1);
                        e = zeros(globalLength,1);
                        mu = zeros(globalLength,1);
                        gammaAux = zeros(globalLength,1);
                        medianAux = zeros(globalLength,1);
                        G = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);
                        xFiltered = zeros(globalLength,1);
                        xLin = zeros(N(NIndex),globalLength);
                        input = randi([0,numberOfSymbols-1],globalLength,1);

                        pilot = pammod(input,pamOrder,0,'gray');

                        pilot2 = pilot.*sqrt(signalPower/var(pilot));

                        xAux = VLC_channel(pilot2, modulationIndexVector(modulationIndex), maxVoltage, VDC, SNR);

                        xAux = [zeros(N(NIndex)-1,1);xAux];

                        w = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;

                        channelIndex = 1;
                        blindFlag = 0;

                        for k = (adapFiltLength(NIndex) + delayVector):globalLength

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

                            if (k > (adapFiltLength(NIndex) + delayVector) + windowLength)
                                medianAux(k) = median(abs(e(k-1 - windowLength:k-1)));
                                if medianAux(k) <= 2*eta(etaIndex) || blindFlag == 1
                                    d(k) = pamHardThreshold(y);

                                    if ~blindFlag 
                                        blindIt(index,delay,NIndex,etaIndex,modulationIndex) = k;
                                    end
                                    blindFlag = 1;

                                else
                                    d(k) = (pilot(-delayVector(delay) + k + 1));
                                end

                            else
                                d(k) = (pilot(-delayVector(delay) + k + 1));
                            end



                            e(k) = d(k) - y;
                            absoluteValueError = abs(e(k));

        %                     gammaAux(k+1) = alpha*gammaAux(k) + (1-alpha)*sqrt(beta*w(:,k)'*w(:,k)*noisePower);
        %                     barGamma = gammaAux(k+1);

                            if absoluteValueError > barGamma
                                mu(k) = 1 - barGamma/absoluteValueError;
                                G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength(NIndex)) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));

                                w(:,k+1) = w(:,k) + mu(k) * G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                                count(k,index) = 1;
                            else
                                w(:,k+1) = w(:,k);
                            end

                        end
                        wIndex(:,:,index) = conj(w(:,1:globalLength));
                        e2(:,index) = abs(e).^2;
                    end

                    meanCount{delay,barGammaNonLinIndex} = mean(count,2);

                    w3{delay,barGammaNonLinIndex} = mean(wIndex,3);
                    e3{delay,barGammaNonLinIndex} = mean(e2,2);

                end
            end
            meanCount2{NIndex,etaIndex,modulationIndex} = meanCount;
            w4{NIndex,etaIndex,modulationIndex} = w3;
            e4{NIndex,etaIndex,modulationIndex} = e3;
        end
    end
end

save(['.' filesep 'results' filesep 'resultsSB_VLC_01.mat'],'w4','e4','meanCount2','blindIt','blindIt2');

rmpath(['..' filesep '..' filesep '..' filesep 'VLC_param' filesep]); 
rmpath(['..' filesep '..' filesep '..' filesep 'VLC_param' filesep 'Utils' filesep]); 
rmpath(['..' filesep '..' filesep  '..' filesep 'VLC_param' filesep 'LED Parameters' filesep]); 
