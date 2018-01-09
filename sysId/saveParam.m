clear;
clc;
close all;


maxRuns = 3000; % max runs in a single independent trial
maxIt = 1000;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

alpha = 0.01;      %forgetting factor of the correlation matrix in SML case

K = 2;             %number of products in the SML case
M = 10;            %length of the adaptiv filter in SML case
mu = 0.8;         %step size

h1 = [0.544 -0.252 0.593 0.236 -0.077 0.156 -0.5 0.025 -0.023 0.099].';
h2 = [-0.204 0.274 0.023 0.024 0.022 -0.274 -0.321 -0.070 0.712 0.433].';

ho = kron(h1,h2); %unknown system in SML case


%-------------------------------------------------------------------------%
%Volterra set-membership

N = 3;

auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);
adapFiltLength = (N^2+N)/2 + N;
kappa = 0.5;
gamma = 1e-12;

wo(:,1) = [0.5 3 5 0 0.3 0 0 1.2 0].';
wo(:,2) = [0.5 3 0 0.5 0.001 0.3 0 0 0].';
wo(:,2) = [2 0 0 0.2 0.3 -0.7 0 0 0].'; 


% wo(:,2) = [1 -2.5 0 0.01 0.007 0.2 0 0 0].';


pamOrder = 4;

changingIteration = 1000;

barGamma = sqrt(5*noisePower); %threshold for set-membership purposes


save(['.' filesep 'simParameters' filesep 'param01.mat']);