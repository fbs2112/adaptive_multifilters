%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers

clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'learningScripts' filesep 'VLC' filesep 'semi-blind' filesep 'results']);
addpath(['..' filesep '..' filesep 'berParameters']);
addpath(['..' filesep '..' filesep 'VLC_param' filesep]);
addpath(['..' filesep '..' filesep 'VLC_param' filesep 'Utils' filesep]);
addpath(['..' filesep '..' filesep 'VLC_param' filesep 'LED Parameters' filesep]);



filterFile = 'resultsSB_VLC_02';

load VLC_param01.mat;
load paramEqVLC.mat;
load param_feedforwardEq.mat;

load([filterFile '.mat']);

modulationIndexVector = [0.05 0.075 0.1];
signalPower = 1;
ber = zeros(length(SNR),size(e4,1),size(e4,3),size(e4,2));


for SNRIndex = 1:length(SNR)
    
    for NIndex = 1:size(e4,1)
        for modulationIndex = 1:size(e4,3)
            modulationIndex
            maxVoltage = VDC*(1+modulationIndexVector(modulationIndex));
            
            for etaIndex = 1:size(e4,2)
                equalizerFilter = [];
                berAux = zeros(monteCarloLoops,1);
                equalizerFilter(:,1) = squeeze(w4{NIndex,etaIndex,modulationIndex}{1,1}(:,end));
                
                for index = 1:monteCarloLoops
                    equalizedSignal = zeros(numberOfSymbols,1);
                    
                    binaryInputData = randi([0,1],blockLength + 500,1);
                    binaryInputData = reshape(binaryInputData,[],numberOfBits);
                    deciInputData = bi2de(binaryInputData);
                    pilot = pammod(deciInputData,2^numberOfBits,0,'gray');
                    pilot = pilot.*sqrt(signalPower/var(pilot));
                    
                    xAux = VLC_channel(pilot, modulationIndexVector(modulationIndex), maxVoltage, VDC, db2pow(SNR(SNRIndex)));
                    channelIndex = 1;
                    
                    xAux = [zeros(N(NIndex)-1,1);xAux];
                    
                    for k = N(NIndex):length(pilot)
                        
                        x = xAux(k:-1:k-N(NIndex)+1,channelIndex);
                        
                        xTDLAux = zeros(length(l1{NIndex}),1);
                        
                        for lIndex = 1:length(l1{NIndex})
                            xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex))*(x(l2{NIndex}(lIndex)));
                        end
                        
                        xAP = [x;xTDLAux];
                        
                        equalizedSignal(k,1) =  equalizerFilter(:,channelIndex)'*xAP;
                        
                    end
                    [corr,lags] = xcorr(equalizedSignal,xAux(N(NIndex):end,1));
                    [~,idx] = max(abs(corr));
                    delay = abs(lags(idx));
                    decDemodSignal = pamdemod(equalizedSignal,pamOrder,0,'gray');
                    binaryOutputData = de2bi(decDemodSignal,numberOfBits);
                    
                    berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
                end
                
                ber(SNRIndex,NIndex,modulationIndex,etaIndex) = mean(berAux);
            end
        end
    end
    
end


save(['.' filesep 'results' filesep  filterFile '_BER.mat'],'SNR','ber');

rmpath(['..' filesep '..' filesep '.' filesep 'learningScripts' filesep 'VLC' filesep 'semi-blind' filesep 'results']);
rmpath(['..' filesep '..' filesep 'berParameters']);
rmpath(['..' filesep '..' filesep 'VLC_param' filesep]);
rmpath(['..' filesep '..' filesep 'VLC_param' filesep 'Utils' filesep]);
rmpath(['..' filesep '..' filesep 'VLC_param' filesep 'LED Parameters' filesep]);
