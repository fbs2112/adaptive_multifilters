%Volterra NLMS DFE

clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'simParameters' filesep]);

load paramDFE_FF.mat;

numberOfSymbols = 2^numberOfBits;

e4 = cell(length(feedforwardLength),length(feedbackLength));
w4 = cell(length(feedforwardLength),length(feedbackLength));
meanCount2 = cell(length(feedforwardLength),length(feedbackLength));


barGammaLin = 4*sqrt(5*noisePower);

barGammaNonLin = 2*barGammaLin;

for FFIndex = 5:length(feedforwardLength)
    FFIndex
    
    for FBIndex = 5:length(feedbackLength)
         FBIndex
         
        if ~volterraFFFlag
            feedforwardLengthNonLin = 0;
        else
            feedforwardLengthNonLin = adaptfiltFF(FFIndex) - feedforwardLength(FFIndex);
        end
        
        if ~volterraFBFlag
            feedbackLengthNonLin = 0;
        else
            feedbackLengthNonLin = adaptfiltFB(FBIndex) - feedbackLength(FBIndex);
        end
 
        feedforwardLengthLin = feedforwardLength(FFIndex);
        feedbackLengthLin = feedbackLength(FBIndex);
%         delayVector = 1:feedforwardLength+length(h);%adapFiltLength + 10;

        delayVector = feedforwardLength(FFIndex)+1;

        e3 = cell(length(delayVector),1);
        w3 = cell(length(delayVector),1);
        meanCount = cell(length(delayVector),1);
        
        for delay = 1:length(delayVector)
            delayVector2 = [feedforwardLength(FFIndex)+1 feedforwardLength(FFIndex)-2];
            globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + max(delayVector2) - 1;

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

                w = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;

                channelIndex = 1;

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

                        auxFF = zeros((feedforwardLength(FFIndex)^2+feedforwardLength(FFIndex))/2,1);

                        for lIndex = 1:length(l1FF{FFIndex})
                            auxFF(lIndex,1) = x(l1FF{FFIndex}(lIndex),k)*(x(l2FF{FFIndex}(lIndex),k));
                        end
                        xConc = [x(:,k);auxFF];
                    else
                        xConc = x(:,k);
                    end


                    if volterraFBFlag
                        auxFB = zeros((feedbackLength(FBIndex)^2+feedbackLength(FBIndex))/2,1);
                        for lIndex = 1:length(l1FB{FBIndex})
                            auxFB(lIndex,1) = yHat(l1FB{FBIndex}(lIndex),k)*(yHat(l2FB{FBIndex}(lIndex),k));
                        end

                        yHatConc = [yHat(:,k);auxFB];
                    else
                        yHatConc = yHat(:,k);
                    end

                    if ~volterraFFFlag && ~volterraFBFlag 
                        xConc = x(:,k);
                        yHatConc = yHat(:,k);
                    end

                    z = [xConc;yHatConc];
                    
                    zLin = [xConc(1:feedforwardLength(FFIndex));yHatConc(1:feedbackLength(FBIndex))];
                    zNonLin = [xConc(feedforwardLength(FBIndex)+1:end);yHatConc(feedbackLength(FBIndex) + 1:end)];

                    d(k) = (pilot(-delayVector(delay) + k + 1));

                    e(k) = d(k) - w(:,k)'*z;
                    
                    absoluteValueError = abs(e(k));

                    auxLength1 = feedforwardLengthLin + feedforwardLengthNonLin + 1;
                    auxLength2 = auxLength1 + feedbackLengthLin - 1;
                    auxLength3 = feedforwardLengthLin + 1;
                    auxLength4 = auxLength3 + feedforwardLengthNonLin - 1; 
                    
                    if absoluteValueError > barGammaLin
                        muLin(k) = 1 - barGammaLin/absoluteValueError;
                        G(1:feedforwardLengthLin,1:feedforwardLengthLin,k) = diag(((1 - kappa*muLin(k))/feedforwardLengthLin) + (kappa*muLin(k)*abs(w(1:feedforwardLengthLin,k))/norm(w(1:feedforwardLengthLin,k),1)));
                        
                        
                        
                        
                        
                        G(auxLength1:auxLength2,auxLength1:auxLength2,k) = diag(((1 - kappa*muLin(k))/(auxLength2 - auxLength1)) + (kappa*muLin(k)*abs(w(auxLength1:auxLength2,k))/...
                            norm(w(auxLength1:auxLength2,k),1)));
                                               
                        w(1:feedforwardLengthLin,k+1) = w(1:feedforwardLengthLin,k) + muLin(k)* G(1:feedforwardLengthLin,1:feedforwardLengthLin,k)*x(:,k)*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
                        
                        w(auxLength1:auxLength2,k+1) = w(auxLength1:auxLength2,k)...
                            + muLin(k)* G(auxLength1:auxLength2,auxLength1:auxLength2,k)...
                            *yHat(:,k)*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));

                        countLin(k,index) = 1;
                    else
                        muLin(k) = 0;
                        w(1:feedforwardLengthLin,k+1) = w(1:feedforwardLengthLin,k);
                        w(auxLength1:auxLength2,k+1) = w(auxLength1:auxLength2,k);
                        G(:,:,k) = eye(adapFiltLength(FFIndex,FBIndex));
                    end
                    
                    
                    
                    if absoluteValueError > barGammaNonLin
                        muNonLin(k) = 1 - barGammaLin/absoluteValueError;
                        G(auxLength3:auxLength4,auxLength3:auxLength4,k) = diag(((1 - kappa*muNonLin(k))/N(5)) + (kappa*muNonLin(k)*abs(w(auxLength3:auxLength4,k))/norm(w(auxLength3:auxLength4,k),1)));
                        
%                         G(auxLength2+1:end,auxLength2+1:end,k) = ...
%                             diag(((1 - kappa*muNonLin(k))/N(5)) + (kappa*muNonLin(k)*abs(w(auxLength2+1:end,k))/norm(w(auxLength2+1:end,k),1)));
                                               
                        w(auxLength3:auxLength4,k+1) = w(auxLength3:auxLength4,k) + muNonLin(k)* G(auxLength3:auxLength4,auxLength3:auxLength4,k)*auxFF*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
                        
%                         w(auxLength2+1:end,k+1) = w(auxLength2+1:end,k)...
%                             + muNonLin(k)* G(auxLength2+1:end,auxLength2+1:end,k)*auxFB*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));

                        countLin(k,index) = 1;
                    else
                        muLin(k) = 0;
                        w(auxLength3:auxLength4,k+1) = w(auxLength3:auxLength4,k);
%                         w(feedforwardLengthLin + feedforwardLengthNonLin + 1:feedbackLengthLin,k+1) = w(feedforwardLengthLin + feedforwardLengthNonLin + 1:feedbackLengthLin,k);
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
        meanCount2{FFIndex,FBIndex} = meanCount;
        w4{FFIndex,FBIndex} = w3;
        e4{FFIndex,FBIndex} = e3;
    end
end

save(['.' filesep 'results' filesep 'results43.mat'],'w4','e4','meanCount2');

rmpath(['..' filesep '..' filesep 'simParameters' filesep]);

