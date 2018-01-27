function IF = I_V_Fun(VF,VT,n,ISat)
%This function returns the LED current as a function of the voltage
%according to the Shockley Equation


aux = VF > 0;
IF = ISat.*(exp(VF./(n*VT)) - 1) .* aux;
