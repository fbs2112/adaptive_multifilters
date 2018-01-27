%Volterra SM OBE DFE

clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'VLC_param' filesep]);

load paramDFE_FF_VLC.mat;
load VLC_param01.mat;

numberOfSymbols = 2^numberOfBits;

modulationIndexVector = [0.05 0.075 0.1];

e4 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));
w4 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));
meanCount2 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));

for modulationIndex = 1:length(modulationIndexVector)
    maxVoltage = VDC*(1+modulationIndexVector(modulationIndex));

    for FFIndex = 1:length(feedforwardLength)
        FFIndex
        for FBIndex = 1:length(feedbackLength)
             FBIndex

            delayVector = feedforwardLength(FFIndex)+1;

            e3 = cell(length(delayVector),1);
            w3 = cell(length(delayVector),1);
            meanCount = cell(length(delayVector),1);

            for delay = 1:length(delayVector)     

                delay

                globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + delayVector - 1;

                wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
                e2 = zeros(globalLength,maxIt);

                count = zeros(globalLength,maxIt);

                for index = 1:maxIt
                    index

                    d = zeros(globalLength,1);
                    P = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);
                    P(:,:,adapFiltLength(FFIndex,FBIndex) + delayVector) = eye(adapFiltLength(FFIndex,FBIndex))*1e-6;
                    sigma = zeros(globalLength,1);
                    sigma(adapFiltLength(FFIndex,FBIndex) + delayVector(delay)) = 1;
                    delta = zeros(globalLength,1);
                    lambda = zeros(globalLength,1);
                    G = zeros(globalLength,1);

                    x = zeros(feedforwardLength(FFIndex),globalLength);
                    yHat = zeros(feedbackLength(FBIndex),1);

                    input = randi([0,numberOfSymbols-1],globalLength,1);
                    pilot = pammod(input,pamOrder,0,'gray');

                    pilot2 = pilot.*sqrt(signalPower/var(pilot));

                    xAux = VLC_channel(pilot, modulationIndexVector(modulationIndex), maxVoltage, VDC, SNR);

                    xAux = [zeros(feedforwardLength(FFIndex)-1,1);xAux];

                    theta = zeros(adapFiltLength(FFIndex,FBIndex),globalLength);

                    channelIndex = 1;

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

                        d(k) = (pilot(-delayVector(delay) + k + 1));

                        delta(k) = d(k) - theta(:,k).'*z;

                        if abs(delta(k)) > barGamma
                            G(k) = z.'*P(:,:,k)*conj(z);
                            lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);
                            lambda2 = 1/lambda(k);

                            P(:,:,k+1) = lambda(k)*(P(:,:,k) - (P(:,:,k)*conj(z)*z.'*P(:,:,k))/(lambda2+G(k)));

                            theta(:,k+1) = theta(:,k) + P(:,:,k+1)*conj(z)*delta(k);

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

             meanCount2{FFIndex,FBIndex,modulationIndex} = meanCount;
            w4{FFIndex,FBIndex,modulationIndex} = w3;
            e4{FFIndex,FBIndex,modulationIndex} = e3;

        end

    end

end

save(['.' filesep 'results' filesep 'results_BEACON_VLC_02.mat'],'w4','e4','meanCount2');

rmpath(['..' filesep '..' filesep 'VLC_param' filesep]);
