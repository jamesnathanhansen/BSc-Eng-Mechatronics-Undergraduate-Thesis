function equ_H = predict_hj(m,g,l,dt,th,dth,ddth)
%PREDICT_HJ
%    EQU_H = PREDICT_HJ(M,G,L,DT,TH,DTH,DDTH)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    09-Nov-2020 00:15:33

equ_H = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,3]);
