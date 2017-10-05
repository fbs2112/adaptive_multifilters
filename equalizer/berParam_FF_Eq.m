clear;
clc;
close all;



numberOfSymbols = 10000;
numberOfBits = 2;
blockLength = numberOfSymbols*numberOfBits;
monteCarloLoops = 1000;
SNR = 0:5:30;


save(['.' filesep 'berParameters' filesep 'param_feedforwardEq.mat']);

