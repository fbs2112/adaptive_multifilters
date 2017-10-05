clear;
clc;
close all;


maxRuns = 8000; % max runs in a single independent trial
maxIt = 1000;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

alpha = 0.01;      %forgetting factor of the correlation matrix in SML case

K = 2;             %number of products in the SML case
M = 10;            %length of the adaptiv filter in SML case
mu = 0.1;         %step size

% h = [1 0.2 -0.3].';


h(:,1) = [0.5 3 5 0 0.3 0 0 1.2 0].';
h(:,2) = [0.5 3 0 0.5 0.001 0.3 0 0 0].';


h1 = [0.544 -0.252 0.593 0.236 -0.077 0.156 -0.5 0.025 -0.023 0.099].';
h2 = [-0.204 0.274 0.023 0.024 0.022 -0.274 -0.321 -0.070 0.712 0.433].';

ho = kron(h1,h2); %unknown system in SML case


%-------------------------------------------------------------------------%
%DFE set-membership

kappa = 0.5;
gamma = 1e-12;

feedforwardLength = 4;
feedbackLength = 4;
   
adapFiltLength = feedforwardLength + feedbackLength;


barGamma = 4*sqrt(5*noisePower); %threshold for set-membership purposes


numberOfBits = 2;
pamOrder = 2^numberOfBits;
changingIteration = 4000;


SNR = db2pow(30);

save(['.' filesep 'simParameters' filesep 'paramDFE.mat']);