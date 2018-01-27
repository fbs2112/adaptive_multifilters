function aux = pamHardThreshold(x)
%This function works as the decicion device in DFE for the case of a 4-PAM
%modulation using gray code

aux = x;

aux(aux >= 0 & aux<2) = 1;

aux(aux >=2) = 3;

aux(aux<0 & aux>= -2) = -1;

aux(aux < -2) = -3;


