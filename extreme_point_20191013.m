%*********************************************************
% This program is to solve the zero point of a function .
% Means the Extreme point of HHmodel
%*********************************************************

clear;
clc;

tm = 0.002;
th = 0.05;
m0 = 0;
h0 = 1;
Imax = 200;
k = 2;
n = 3;
mgig = 0.8;
hgig = 0.5;

%----------------------------------------------------
% Fail to find the zero point of derivative
%syms t;
%eqn = (1-hgig)*(-1/th)*exp(-t/th)+(hgig/tm)*exp(-t/tm)-(1/tm+1/th)*(hgig-1)*exp(-(1/tm+1/th)*t) == 0;
%t0 = solve(eqn,t);

figure;hold on;
t = 0:0.01:1;
for i=1:100
    tm = i/1000;
    I = Isimulate(mgig,hgig,m0,h0,k,n,tm,th,Imax,t);
    plot(t,I);
end




