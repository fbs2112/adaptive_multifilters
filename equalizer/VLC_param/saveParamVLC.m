clear;
clc;
close all;


maxRuns = 5000; % max runs in a single independent trial
maxIt = 1000;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

%-------------------------------------------------------------------------%
%Volterra set-membership

N = 7:9;

kappa = 0.5;
gamma = 1e-8;
lambdaUp = 0.5;
alpha = 0.6;
beta = 70;

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

SNR = db2pow(30);

save('paramEqVLC.mat');