function [decidedSymbol] = signal2Symb(signal,signalPower)

qamSymbols = [1+1i 1-1i -1+1i -1-1i].';


powerqamSymbols = qamSymbols'*qamSymbols/(length(qamSymbols));


qamSymbols = qamSymbols.*sqrt(signalPower/powerqamSymbols);

qamSymbols2 = [real(qamSymbols) imag(qamSymbols)].';

signalAux = [real(signal);imag(signal)];

aux = repmat(signalAux,1,4);

distance = sqrt(sum((qamSymbols2 - aux).^2));


[~,minDistanceIndex] = min(distance);

decidedSymbol = qamSymbols(minDistanceIndex);




% realSymbol = (real(signal) - real(qamSymbols)).^2;
% 
% imagSymbol = (imag(signal) - imag(qamSymbols)).^2;
% 
% 
% [~,indexReal] = sort(realSymbol);
% 
% [~,indexImag] = sort(imagSymbol);


% [~,decidedSymbolIndex] = min([realSymbol imagSymbol]);

