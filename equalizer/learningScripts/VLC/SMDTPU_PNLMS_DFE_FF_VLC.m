%Volterra NLMS DFE

clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'VLC_param' filesep]); 
addpath(['..' filesep '..' filesep 'VLC_param' filesep 'Utils' filesep]); 
addpath(['..' filesep '..' filesep 'VLC_param' filesep 'LED Parameters' filesep]); 

load paramDFE_FF_VLC.mat;
load VLC_param01.mat;

numberOfSymbols = 2^numberOfBits;
lambdaUp = 0.5;

barGammaLin = 4*sqrt(5*noisePower);
barGammaNonLin = (1:0.5:4)*barGammaLin;
modulationIndexVector = [0.05 0.075 0.1];

e4 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));
w4 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));
meanCountLin2 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));
meanCountNonLin2 = cell(length(feedforwardLength),length(feedbackLength),length(modulationIndexVector));

for modulationIndex = 1:length(modulationIndexVector)

    maxVoltage = VDC*(1+modulationIndexVector(modulationIndex));


    for FFIndex = 1:length(feedforwardLength)
        FFIndex
        for FBIndex = 1:length(feedbackLength)
            FBIndex

            CLin = diag([ones(feedbackLength(FFIndex),1).' zeros(adaptfiltFF(FFIndex) - feedbackLength(FFIndex),1).'...
                ones(feedbackLength(FBIndex),1).' zeros(adaptfiltFB(FBIndex) - feedbackLength(FBIndex),1).'].');
            CNonLin = diag([zeros(feedbackLength(FFIndex),1).' ones(adaptfiltFF(FFIndex) - feedbackLength(FFIndex),1).'...
                zeros(feedbackLength(FBIndex),1).' ones(adaptfiltFB(FBIndex) - feedbackLength(FBIndex),1).'].');  

            delayVector = feedforwardLength(FFIndex)+1;

            e3 = cell(length(delayVector),1);
            w3 = cell(length(delayVector),1);
            meanCountLin = cell(length(delayVector),length(barGammaNonLin));
            meanCountNonLin = cell(length(delayVector),length(barGammaNonLin));

            for delay = 1:length(delayVector)
                globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + delayVector - 1;



                for barGammaNonLinIndex = 1:length(barGammaNonLin)
                    barGammaNonLinIndex
                    wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
                    e2 = zeros(globalLength,maxIt);

                    countLin = zeros(globalLength,maxIt);
                    countNonLin = zeros(globalLength,maxIt);
                    countTotal = zeros(globalLength,maxIt);

                    for index = 1:maxIt
                        index
    
                        mu = zeros(globalLength,1);
                        d = zeros(globalLength,1);
                        e = zeros(globalLength,1);
                        muLin = zeros(globalLength,1);
                        muNonLin = zeros(globalLength,1);
                        G = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);
                        G1 = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);
                        G2 = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);

                        x = zeros(feedforwardLength(FFIndex),globalLength);
                        yHat = zeros(feedbackLength(FBIndex),1);

                        input = randi([0,numberOfSymbols-1],globalLength,1);

                        pilot = pammod(input,pamOrder,0,'gray');

                        pilot2 = pilot.*sqrt(signalPower/var(pilot));

                        xAux = VLC_channel(pilot2, modulationIndexVector(modulationIndex), maxVoltage, VDC, SNR);

                        xAux = [zeros(feedforwardLength(FFIndex)-1,1);xAux];

                        w = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;

                        w1 = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;
                        w2 = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;

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

                            if absoluteValueError > barGammaLin
                                muLin(k) = 1 - barGammaLin/absoluteValueError;
                                G1(:,:,k) = diag(((1 - kappa*muLin(k))/adapFiltLength(FFIndex,FBIndex)) + (kappa*muLin(k)*abs(w1(:,k))/norm(w1(:,k),1)));

                                w1(:,k+1) = w1(:,k) + muLin(k) * CLin*G1(:,:,k)*z*((z'*CLin*G1(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
                                countLin(k,index) = 1;
                            else
                                w1(:,k+1) = w1(:,k);
                            end


                            if absoluteValueError > barGammaNonLin(barGammaNonLinIndex)
                                muNonLin(k) = 1 - barGammaNonLin(barGammaNonLinIndex)/absoluteValueError;

                                G2(:,:,k) = diag(((1 - kappa*muNonLin(k))/adapFiltLength(FFIndex,FBIndex)) + (kappa*muNonLin(k)*abs(w2(:,k))/norm(w2(:,k),1)));

                                w2(:,k+1) = w2(:,k) + muNonLin(k)*CNonLin*G2(:,:,k)*z*((z'*CNonLin*G2(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
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
            meanCountLin2{FFIndex,FBIndex,modulationIndex} = meanCountLin;
            meanCountNonLin2{FFIndex,FBIndex,modulationIndex} = meanCountNonLin;
            w4{FFIndex,FBIndex,modulationIndex} = w3;
            e4{FFIndex,FBIndex,modulationIndex} = e3;
        end
    end
    
end

save(['.' filesep 'results' filesep 'results_SMDTPU_VLC_02.mat'],'w4','e4','meanCountLin2','meanCountNonLin2');

rmpath(['..' filesep '..' filesep 'VLC_param' filesep]); 
rmpath(['..' filesep '..' filesep 'VLC_param' filesep 'Utils' filesep]); 
rmpath(['..' filesep '..' filesep 'VLC_param' filesep 'LED Parameters' filesep]); 
