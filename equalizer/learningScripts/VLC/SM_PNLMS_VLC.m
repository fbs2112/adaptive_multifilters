%Volterra Data Reuse Equalyzer 

clear;
clc;
close all;


addpath(['..' filesep '..' filesep 'VLC_param' filesep]); 
addpath(['..' filesep '..' filesep 'VLC_param' filesep 'Utils' filesep]); 
addpath(['..' filesep '..' filesep 'VLC_param' filesep 'LED Parameters' filesep]); 


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
        NIndex

        delayVector = N(NIndex)+1;%adapFiltLength + 10;
        e3 = cell(length(delayVector),1);
        w3 = cell(length(delayVector),1);
        meanCount = cell(length(delayVector),1);

        for delay = 1:length(delayVector)


            globalLength = maxRuns + adapFiltLength(NIndex) + delayVector - 1;

            wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
            e2 = zeros(globalLength,maxIt);

            count = zeros(globalLength,maxIt);

            for index = 1:maxIt
                index


                mu = zeros(globalLength,1);
                d = zeros(globalLength,1);
                e = zeros(globalLength,1);
                G = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);

                x = zeros(N(NIndex),globalLength);

                input = randi([0,numberOfSymbols-1],globalLength,1);

                pilot = pammod(input,pamOrder,0,'gray');

                pilot2 = pilot.*sqrt(signalPower/var(pilot));


                xAux = VLC_channel(pilot2, modulationIndexVector(modulationIndex), maxVoltage, VDC, SNR);

                xAux = [zeros(N(NIndex)-1,1);xAux];

                w = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;

                channelIndex = 1;

               for k = (adapFiltLength(NIndex) + delayVector):globalLength

                    x(:,k) = xAux(k:-1:k-N(NIndex)+1,channelIndex);

                    xTDLAux = zeros(length(l1{NIndex}),1);

                    for lIndex = 1:length(l1{NIndex})
                        xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex),k)*(x(l2{NIndex}(lIndex),k));
                    end

                    xAP = [x(:,k);xTDLAux];

                    d(k) = (pilot(-delayVector(delay) + k + 1)); 

                    e(k) = d(k) - w(:,k)'*xAP(:,1);

                    absoluteValueError = abs(e(k));

                    if absoluteValueError > barGamma
                        mu(k) = 1 - barGamma/absoluteValueError;
                        G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength(NIndex)) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                        w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                        count(k,index) = 1;
                    else
                        mu(k) = 0;
                        w(:,k+1) = w(:,k);
                        G(:,:,k) = eye(adapFiltLength(NIndex));
                    end

                end
                wIndex(:,:,index) = conj(w(:,1:globalLength));
                e2(:,index) = abs(e).^2;
            end

            meanCount{delay} = mean(count,2);

            w3{delay} = mean(wIndex,3);

            e3{delay} = mean(e2,2);

        end

        meanCount2{NIndex,modulationIndex} = meanCount;
        w4{NIndex,modulationIndex} = w3;
        e4{NIndex,modulationIndex} = e3;

    end
end

save(['.' filesep 'results' filesep 'results_SMPNLMS_VLC_01.mat'],'w4','e4','meanCount2');

rmpath(['..' filesep '..' filesep 'VLC_param' filesep]); 
rmpath(['..' filesep '..' filesep 'VLC_param' filesep 'Utils' filesep]); 
rmpath(['..' filesep '..' filesep 'VLC_param' filesep 'LED Parameters' filesep]); 


