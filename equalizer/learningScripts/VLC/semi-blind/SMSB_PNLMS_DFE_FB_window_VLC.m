%Volterra NLMS DFE

clear;
clc;
close all;

addpath(['..' filesep '..' filesep '..' filesep 'VLC_param' filesep]);

load paramDFE_FB_VLC.mat;
load VLC_param01.mat;

numberOfSymbols = 2^numberOfBits;
barGammaNonLin = 1;
windowLength = 100;

modulationIndexVector = [0.05 0.075 0.1];
eta = 0:0.1:0.3;

e4 = cell(length(feedforwardLength),length(feedbackLength),length(eta));
w4 = cell(length(feedforwardLength),length(feedbackLength),length(eta));
meanCount2 = cell(length(feedforwardLength),length(feedbackLength),length(eta));

blindIt = zeros(maxIt,1,length(feedforwardLength),length(feedbackLength),length(eta));
blindIt2 = zeros(maxIt,1,length(feedforwardLength),length(feedbackLength),length(eta));

maxIt = 50;


for modulationIndex = 1:length(modulationIndexVector)
    
    maxVoltage = VDC*(1+modulationIndexVector(modulationIndex));

    for etaIndex = 1:length(eta)

        for FFIndex = 1:length(feedforwardLength)
            FFIndex
            for FBIndex = 5:length(feedbackLength)
                FBIndex
                %         delayVector = 1:feedforwardLength+length(h);%adapFiltLength + 10;

                delayVector = feedforwardLength(FFIndex)+1;

                e3 = cell(length(delayVector),1);
                w3 = cell(length(delayVector),1);
                meanCount = cell(length(delayVector),length(barGammaNonLin));

                for delay = 1:length(delayVector)
                    globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + delayVector - 1;


                    for barGammaNonLinIndex = 1:length(barGammaNonLin)
                        barGammaNonLinIndex
                        wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
                        e2 = zeros(globalLength,maxIt);

                        count = zeros(globalLength,maxIt);

                        for index = 1:maxIt
                            index

                            mu = zeros(globalLength,1);
                            d = zeros(globalLength,1);
                            e = zeros(globalLength,1);
                            medianAux = zeros(globalLength,1);
                            G = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);

                            x = zeros(feedforwardLength(FFIndex),globalLength);
                            yHat = zeros(feedbackLength(FBIndex),1);

                            input = randi([0,numberOfSymbols-1],globalLength,1);

                            pilot = pammod(input,pamOrder,0,'gray');

                            pilot2 = pilot.*sqrt(signalPower/var(pilot));

                            xAux = VLC_channel(pilot2, modulationIndexVector(modulationIndex), maxVoltage, VDC, SNR);

                            xAux = [zeros(feedforwardLength(FFIndex)-1,1);xAux];
                            
                            w = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;

                            channelIndex = 1;
                            blindFlag = 0;

                            for k = (adapFiltLength(FFIndex,FBIndex) + delayVector):globalLength

                                x(:,k) = xAux(k:-1:k-feedforwardLength(FFIndex)+1,channelIndex);

                                yHat(:,k) = (pilot(-delayVector(delay) + k + 1 -1:-1:-delayVector(delay) + k + 1 - feedbackLength(FBIndex) - 1 + 1));

                                if volterraFFFlag

                                    aux = zeros((feedforwardLength(FFIndex)^2+feedforwardLength(FFIndex))/2,1);

                                    for lIndex = 1:length(l1FF{FFIndex})
                                        aux(lIndex,1) = x(l1FF{FFIndex}(lIndex),k)*(x(l2FF{FFIndex}(lIndex),k));
                                    end
                                    xConc = [x(:,k);aux];
                                else
                                    xConc = x(:,k);
                                end


                                if volterraFBFlag
                                    aux = zeros((feedbackLength(FBIndex)^2+feedbackLength(FBIndex))/2,1);
                                    for lIndex = 1:length(l1FB{FBIndex})
                                        aux(lIndex,1) = yHat(l1FB{FBIndex}(lIndex),k)*(yHat(l2FB{FBIndex}(lIndex),k));
                                    end

                                    yHatConc = [yHat(:,k);aux];
                                else
                                    yHatConc = yHat(:,k);
                                end

                                if ~volterraFFFlag && ~volterraFBFlag
                                    xConc = x(:,k);
                                    yHatConc = yHat(:,k);
                                end

                                z = [xConc;yHatConc];
                                y = w(:,k)'*z;

                                if k > adapFiltLength(FFIndex,FBIndex) + delayVector + windowLength 
                                    medianAux(k) = median(abs(e(k-1 - windowLength:k-1)));
                                    if medianAux(k) <= 2*eta(etaIndex) || blindFlag == 1
                                        d(k) = pamHardThreshold(y);

                                        if ~blindFlag 
                                            blindIt(index,delay,FFIndex,FBIndex,etaIndex) = k;
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

                                if absoluteValueError > barGamma
                                    mu(k) = 1 - barGamma/absoluteValueError;
                                    G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength(FFIndex,FBIndex)) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));

                                    w(:,k+1) = w(:,k) + mu(k) * G(:,:,k)*z*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
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
                meanCount2{FFIndex,FBIndex,etaIndex} = meanCount;
                w4{FFIndex,FBIndex,etaIndex} = w3;
                e4{FFIndex,FBIndex,etaIndex} = e3;
            end
        end
    end
end

save(['.' filesep 'results' filesep 'resultsTestDFE_FB.mat'],'w4','e4','meanCount2','blindIt','blindIt2');

rmpath(['..' filesep '..' filesep '..' filesep 'VLC_param' filesep]);

