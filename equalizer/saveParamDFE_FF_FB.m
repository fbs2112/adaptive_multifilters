clear;
clc;
close all;


maxRuns = 10000; % max runs in a single independent trial
maxIt = 1000;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

alpha = 0.01;      %forgetting factor of the correlation matrix in SML case

K = 2;             %number of products in the SML case
M = 10;            %length of the adaptiv filter in SML case
mu = 0.5;         %step size

h(:,1) = [0.5 3 5 0 0.3 0 0 1.2 0].';
h(:,2) = [0.5 3 0 0.5 0.001 0.3 0 0 0].';

h(:,2) = [2 0 0 0.2 0.3 -0.7 0 0 0].'; 

h1 = [0.544 -0.252 0.593 0.236 -0.077 0.156 -0.5 0.025 -0.023 0.099].';
h2 = [-0.204 0.274 0.023 0.024 0.022 -0.274 -0.321 -0.070 0.712 0.433].';


ho = kron(h1,h2); %unknown system in SML case


%-------------------------------------------------------------------------%
%DFE set-membership

kappa = 0.5;
gamma = 1e-12;


memoryChannelLength = 3;

volterraFFFlag = 1;
volterraFBFlag = 1;

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

save(['.' filesep 'simParameters' filesep 'paramDFE_FF_FB.mat']);