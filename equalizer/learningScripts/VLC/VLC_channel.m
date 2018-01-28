function [receivedVoltageSignal] = VLC_channel(signal, modulationIndex, maxVoltage, VDC, SNR)


load whiteLED_334-15_Param.mat;
load whiteLED_334-15.mat;
load VLC_param01.mat;

convLength = length(signal) + LEDfreqRespPoints -1;
NFFT = 2^nextpow2(convLength);

signalFreq = fft(signal,NFFT);

f = fs/2*linspace(0,1,NFFT/2 + 1)*2*pi;

w = [-fliplr(f(2:end-1)) f];

LEDResp = freqRespLED(w);

filteredVinAux = real(ifft(signalFreq.*fftshift(LEDResp)));

filteredVin = filteredVinAux(1:length(signal));

VoltageConstant = modulationIndex*maxVoltage/((1+modulationIndex)*max(filteredVin));

filteredVin = filteredVin*VoltageConstant + VDC;

iLEDOutput = I_V_Fun(filteredVin,VT,nLED,ISat);

eletricalPowerOutput = filteredVin.*iLEDOutput;

opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,kNonLinearity,maxLuminousIntensityLED);

opticalPowerOutputConvolved = opticalPowerOutput*H_0;

n = randn(length(opticalPowerOutputConvolved),1); %noise signal

receivedCurrentSignal = opticalPowerOutputConvolved*R*A;
receivedCurrentSignalAC = receivedCurrentSignal - mean(receivedCurrentSignal);
receivedCurrentSignalPower = receivedCurrentSignalAC'*receivedCurrentSignalAC/length(receivedCurrentSignal);

powerNoiseAux = n'*n/(length(n));
powerNoise = receivedCurrentSignalPower/SNR;
n = n.*sqrt(powerNoise/powerNoiseAux);

receivedVoltageSignalAux = (receivedCurrentSignal + n);
receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(signal)/var(receivedVoltageSignalAux));




