%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers

clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'learningScripts' filesep 'VLC' filesep 'results']);
addpath(['..' filesep 'berParameters']);
addpath(['..' filesep 'VLC_param' filesep]);
addpath(['..' filesep 'VLC_param' filesep 'Utils' filesep]);
addpath(['..' filesep 'VLC_param' filesep 'LED Parameters' filesep]);
addpath(['.' filesep 'Misc']);


filterFile = 'results_SMDTPU_VLC_02';


load VLC_param01.mat;
load paramDFE_FF_VLC.mat;
load param_feedforwardEq.mat;

load([filterFile '.mat']);

modulationIndexVector = [0.05 0.075 0.1];


ber = zeros(length(SNR),size(e4,1),size(e4,2),size(e4,3),size(w4{1,1},2));

for SNRIndex = 1:length(SNR)
    for FFIndex = 1:size(e4,1)
        FFIndex
        for FBIndex = 1:size(e4,2)
            FBIndex
            
            for modulationIndex = 1:length(modulationIndexVector)
                maxVoltage = VDC*(1+modulationIndexVector(modulationIndex));
                
                for barGammaNonLinIndex = 1:size(w4{1,1},2)
                    equalizerFilter = [];
                    berAux = zeros(monteCarloLoops,1);
                    equalizerFilter(:,1) = squeeze(w4{FFIndex,FBIndex,modulationIndex}{1,barGammaNonLinIndex}(:,end));
                    barGammaNonLinIndex
                    for index = 1:monteCarloLoops
%                         index
                        equalizedSignal = zeros(numberOfSymbols,1);
                        
                        binaryInputData = randi([0,1],blockLength + 100,1);
                        binaryInputData = reshape(binaryInputData,[],numberOfBits);
                        deciInputData = bi2de(binaryInputData);
                        pilot = pammod(deciInputData,2^numberOfBits,0,'gray');
                        pilotVar = var(pilot);
                        pilot = pilot.*sqrt(signalPower/var(pilot));
                        
                        xAux = VLC_channel(pilot, modulationIndexVector(modulationIndex), maxVoltage, VDC, db2pow(SNR(SNRIndex)));
                        
                        channelIndex = 1;
                        
                        xAux = [zeros(feedforwardLength(FFIndex)-1,1);xAux];
                        outputFF = zeros(length(pilot),1);
                        outputFB = zeros(length(pilot),1);
                        
                        
                        for k = feedforwardLength(FFIndex) + feedbackLength(FBIndex) + 1:length(pilot)
                            
                            x = xAux(k:-1:k-feedforwardLength(FFIndex)+1,channelIndex);
                            
                            if volterraFFFlag
                                
                                aux = zeros((feedforwardLength(FFIndex)^2+feedforwardLength(FFIndex))/2,1);
                                
                                for lIndex = 1:length(l1FF{FFIndex})
                                    aux(lIndex,1) = x(l1FF{FFIndex}(lIndex),1)*(x(l2FF{FFIndex}(lIndex),1));
                                end
                                xConc = [x(:,1);aux];
                            else
                                xConc = x(:,1);
                            end
                            outputFF(k) =  (equalizerFilter(1:adaptfiltFF(FFIndex),channelIndex))'*xConc;
                            equalizedSignal(k,1) = pamHardThreshold(outputFF(k) + outputFB(k-1));
                            
                            inputFB = equalizedSignal(k:-1:k-feedbackLength(FBIndex) + 1);
                            
                            if volterraFBFlag
                                aux = zeros((feedbackLength(FBIndex)^2+feedbackLength(FBIndex))/2,1);
                                for lIndex = 1:length(l1FB{FBIndex})
                                    aux(lIndex,1) = inputFB(l1FB{FBIndex}(lIndex),1)*(inputFB(l2FB{FBIndex}(lIndex),1));
                                end
                                
                                yHatConc = [inputFB(:,1);aux];
                            else
                                yHatConc = inputFB(:,1);
                            end
                            
                            if ~volterraFFFlag && ~volterraFBFlag
                                xConc = x(:,1);
                                yHatConc = inputFB;
                            end
                            outputFB(k) = (equalizerFilter(adaptfiltFF(FFIndex)+1:end,channelIndex))'*yHatConc;
                            
                            
                        end
                        [corr,lags] = xcorr(equalizedSignal,xAux(feedforwardLength(FFIndex):end,1));
                        [~,idx] = max(abs(corr));
                        delay = abs(lags(idx));
                        
                        decDemodSignal = pamdemod(equalizedSignal,pamOrder,0,'gray');
                        
                        binaryOutputData = de2bi(decDemodSignal,numberOfBits);
                        
                        berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
                    end
                    
                    ber(SNRIndex,FFIndex,FBIndex,modulationIndex,barGammaNonLinIndex) = mean(berAux);
                end
            end
        end
    end
end

% 
% save(['.' filesep 'results' filesep  filterFile '_BER.mat'],'SNR','ber');

rmpath(['..' filesep '.' filesep 'learningScripts' filesep 'VLC' filesep 'results']);
rmpath(['..' filesep 'berParameters']);
rmpath(['..' filesep 'VLC_param' filesep]);
rmpath(['..' filesep 'VLC_param' filesep 'Utils' filesep]);
rmpath(['..' filesep 'VLC_param' filesep 'LED Parameters' filesep]);
rmpath(['.' filesep 'Misc']);


