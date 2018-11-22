function [y]=sfun(x,T)
phi=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
y=phi*x;