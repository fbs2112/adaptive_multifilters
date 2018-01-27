%Volterra NLMS DFE

clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'VLC_param' filesep]);

load paramDFE_FF_FB_VLC.mat;
load VLC_param01.mat;

numberOfSymbols = 2^numberOfBits;
modulationIndexVector = [0.05 0.075 0.1];

e4 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));
w4 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));
meanCount2 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));
maxIt = 50;

for modulationIndex = 1:length(modulationIndexVector)
    
    maxVoltage = VDC*(1+modulationIndexVector(modulationIndex));

    for FFIndex = 1:length(feedforwardLength)
        FFIndex
        for FBIndex = 1:1%length(feedbackLength)
             FBIndex

            delayVector = feedforwardLength(FFIndex)+1;

            e3 = cell(length(delayVector),1);
            w3 = cell(length(delayVector),1);
            meanCount = cell(length(delayVector),1);

            for delay = 1:length(delayVector)
                globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + delayVector - 1;

                wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
                e2 = zeros(globalLength,maxIt);
                count = zeros(globalLength,maxIt);

                for index = 1:maxIt
                    index

                    mu = zeros(globalLength,1);
                    d = zeros(globalLength,1);
                    e = zeros(globalLength,1);
                    G = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);

                    x = zeros(feedforwardLength(FFIndex),globalLength);
                    yHat = zeros(feedbackLength(FBIndex),1);

                    input = randi([0,numberOfSymbols-1],globalLength,1);

                    pilot = pammod(input,pamOrder,0,'gray');

                    pilot2 = pilot.*sqrt(signalPower/var(pilot));

                    xAux = VLC_channel(pilot, modulationIndexVector(modulationIndex), maxVoltage, VDC, SNR);

                    xAux = [zeros(feedforwardLength(FFIndex)-1,1);xAux];

                    w = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;

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

                        e(k) = d(k) - w(:,k)'*z;

                        absoluteValueError = abs(e(k));

                        if absoluteValueError > barGamma
                            mu(k) = 1 - barGamma/absoluteValueError;
                            G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength(FFIndex,FBIndex)) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                            w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*z*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
                            count(k,index) = 1;
                        else
                            mu(k) = 0;
                            w(:,k+1) = w(:,k);
                            G(:,:,k) = eye(adapFiltLength(FFIndex,FBIndex));
                        end

                    end
                    wIndex(:,:,index) = conj(w(:,1:globalLength));
                    e2(:,index) = abs(e).^2;
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

save(['.' filesep 'results' filesep 'resultsTestDFE_FF_FB.mat'],'w4','e4','meanCount2');

rmpath(['..' filesep '..' filesep 'VLC_param' filesep]);

