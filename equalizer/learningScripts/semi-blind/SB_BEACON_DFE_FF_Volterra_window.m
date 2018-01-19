%Volterra NLMS DFE

clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'simParameters' filesep]);

load paramDFE_FF.mat;

numberOfSymbols = 2^numberOfBits;
barGammaNonLin = 1;
windowLength = 100;

eta = 0.1:0.1:0.3;

e4 = cell(length(feedforwardLength),length(feedbackLength),length(eta));
w4 = cell(length(feedforwardLength),length(feedbackLength),length(eta));
meanCount2 = cell(length(feedforwardLength),length(feedbackLength),length(eta));

blindIt = zeros(maxIt,1,length(feedforwardLength),length(feedbackLength),length(eta));
blindIt2 = zeros(maxIt,1,length(feedforwardLength),length(feedbackLength),length(eta));

for etaIndex = 1:length(eta)

    for FFIndex = 1:length(feedforwardLength)
        FFIndex
        for FBIndex = 1:length(feedbackLength)
            FBIndex
            %         delayVector = 1:feedforwardLength+length(h);%adapFiltLength + 10;

            delayVector = feedforwardLength(FFIndex)+1;

            e3 = cell(length(delayVector),1);
            w3 = cell(length(delayVector),1);
            meanCount = cell(length(delayVector),length(barGammaNonLin));

            for delay = 1:length(delayVector)
                delayVector2 = [feedforwardLength(FFIndex)+1 feedforwardLength(FFIndex)-2];
                globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + max(delayVector2) - 1;


                for barGammaNonLinIndex = 1:length(barGammaNonLin)
                    barGammaNonLinIndex
                    wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
                    e2 = zeros(globalLength,maxIt);

                    count = zeros(globalLength,maxIt);

                    for index = 1:maxIt
                        index

                        mu = zeros(globalLength,1);
                        d = zeros(globalLength,1);
                        G = zeros(globalLength,1);
                        delta = zeros(globalLength,1);
                        lambda = zeros(globalLength,1);
                        medianAux = zeros(globalLength,1);
                        P = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);
                        P(:,:,(adapFiltLength(FFIndex,FBIndex) + max(delayVector2))) = eye(adapFiltLength(FFIndex,FBIndex))*1e-6;

                        x = zeros(feedforwardLength(FFIndex),globalLength);
                        yHat = zeros(feedbackLength(FBIndex),1);

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

                        theta = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;

                        channelIndex = 1;
                        blindFlag = 0;

                        for k = (adapFiltLength(FFIndex,FBIndex) + max(delayVector2)):globalLength

                            if k >= changingIteration
                                if feedforwardLength(FFIndex) > 2
                                    delayVector = delayVector2(2);
                                else
                                    delayVector = delayVector2(2) + 2;
                                end
                                channelIndex = 2;
                            else
                                delayVector = delayVector2(1);
                                channelIndex = 1;
                            end

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
                            y = theta(:,k)'*z;

                            if (k > (adapFiltLength(FFIndex,FBIndex) + max(delayVector2)) + windowLength && k < changingIteration) || k > changingIteration + windowLength
                                medianAux(k) = median(abs(delta(k-1 - windowLength:k-1)));
                                if medianAux(k) <= 2*eta(etaIndex) || blindFlag == 1
                                    d(k) = pamHardThreshold(y);

                                    if ~blindFlag && k < changingIteration
                                        blindIt(index,delay,FFIndex,FBIndex,etaIndex) = k;
                                    elseif ~blindFlag && k > changingIteration
                                        blindIt2(index,delay,FFIndex,FBIndex,etaIndex) = k;
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

                            delta(k) = d(k) - y;

                            absoluteValueError = abs(delta(k));

                            if absoluteValueError > barGamma
                                G(k) = z.'*P(:,:,k)*conj(z);
                                lambda(k) = (1/G(k))*((absoluteValueError/barGamma) - 1);
                                lambda2 = 1/lambda(k);

                                P(:,:,k+1) = lambda(k)*(P(:,:,k) - (P(:,:,k)*conj(z)*z.'*P(:,:,k))/(lambda2+G(k)));

                                theta(:,k+1) = theta(:,k) + P(:,:,k+1)*conj(z)*delta(k);

                                count(k,index) = 1;
                            else
                                theta(:,k+1) = theta(:,k);
                                P(:,:,k+1) = P(:,:,k);
                            end

                        end
                        wIndex(:,:,index) = conj(theta(:,1:globalLength));
                        e2(:,index) = abs(delta).^2;
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

save(['.' filesep 'results' filesep 'SB_results04.mat'],'w4','e4','meanCount2','blindIt','blindIt2');

rmpath(['..' filesep '..' filesep 'simParameters' filesep]);

