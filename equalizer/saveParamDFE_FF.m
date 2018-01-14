clear;
clc;
close all;


maxRuns = 10000; % max runs in a single independent trial
maxIt = 1000;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

h(:,1) = [0.5 3 5 0 0.3 0 0 1.2 0].';
h(:,2) = [2 0 0 0.2 0.3 -0.7 0 0 0].'; 
%-------------------------------------------------------------------------%
%DFE set-membership

kappa = 0.5;
gamma = 1e-8;
lambdaUp = 0.5;
alpha = 0.6;
beta = 70;

memoryChannelLength = 3;

volterraFFFlag = 1;
volterraFBFlag = 0;


feedforwardLength = 1:5;
feedbackLength = 1:5;


% feedforwardLength = 6;
% feedbackLength = 6;

adaptfiltFF = zeros(length(feedforwardLength),1);
l1FF = cell(length(feedforwardLength),1);
l2FF = cell(length(feedforwardLength),1);

l1FB = cell(length(feedbackLength),1);
l2FB = cell(length(feedbackLength),1);

adaptfiltFB = zeros(length(feedbackLength),1);
adapFiltLength = zeros(length(feedforwardLength),length(feedbackLength));

for i = 1:length(feedforwardLength)
    adaptfiltFF(i) = (feedforwardLength(i)^2+feedforwardLength(i))/2 + feedforwardLength(i);
    auxMatrix = triu(ones(feedforwardLength(i)));
    [l1FF{i},l2FF{i}] = find(auxMatrix);
    for j = 1:length(feedbackLength)

        
        adaptfiltFB(j) = (feedbackLength(j)^2+feedbackLength(j))/2 + feedbackLength(j);
        

        auxMatrix = triu(ones(feedbackLength(j)));
        [l1FB{j},l2FB{j}] = find(auxMatrix);

        

        if ~volterraFFFlag
            adaptfiltFF(i) = feedforwardLength(i);
        end

        if ~volterraFBFlag
            adaptfiltFB(j) = feedbackLength(j);
        end


        adapFiltLength(i,j) = adaptfiltFF(i) + adaptfiltFB(j);
    end
end

auxMatrix = triu(ones(memoryChannelLength));
[l1Pilot,l2Pilot] = find(auxMatrix);


barGamma = 4*sqrt(5*noisePower); %threshold for set-membership purposes


numberOfBits = 2;
changingIteration = 5000;

pamOrder = 2^numberOfBits;

SNR = db2pow(30);

save(['.' filesep 'simParameters' filesep 'paramDFE_FF.mat']);