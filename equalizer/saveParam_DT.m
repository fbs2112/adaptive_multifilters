clear;
clc;
close all;


maxRuns = 7000; % max runs in a single independent trial
maxIt = 1000;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

h(:,1) = [0.5 3 5 0 0.3 0 0 1.2 0].';
h(:,2) = [2 0 0 0.2 0.3 -0.7 0 0 0].'; 

%-------------------------------------------------------------------------%
%Volterra set-membership

N = 7:9;

kappa = 0.5;
gamma = 1e-8;
lambdaUp = 0.5;
alpha = 0.6;
beta = 70;

% N = 6;

l1 = cell(length(N),1);
l2 = cell(length(N),1);
adapFiltLength = zeros(length(N),1);

for i = 1:length(N)
    auxMatrix = triu(ones(N(i)));
    [l1{i},l2{i}] = find(auxMatrix);
    adapFiltLength(i) = (N(i)^2+N(i))/2 + N(i);
end

barGamma = 4*sqrt(5*noisePower); %threshold for set-membership purposes

numberOfBits = 2;

pamOrder = 2^numberOfBits;

memoryChannelLength = 3;

auxMatrix = triu(ones(memoryChannelLength));
[l1Pilot,l2Pilot] = find(auxMatrix);


changingIteration = 10000;


SNR = db2pow(30);

save(['.' filesep 'simParameters' filesep 'paramEq_DT.mat']);
