%##############################################################
% Read the Imax data caused by force with different V &D
% Saved as a matrix file
% 
%##############################################################
clear;
clc;

[d,scale,hi] = abfload('./data/Imax/2019_05_14_0108.abf');
data(1).currentData = d;
data(1).si = scale;
data(1).h = hi;
[d,scale,hi] = abfload('./data/Imax/2019_05_14_0109.abf');
data(2).currentData = d;
data(2).si = scale;
data(2).h = hi;

for i = 10:16
    a = num2str(i);
    [d,scale,hi] = abfload(strcat('./data/Imax/2019_05_14_01',a,'.abf'));
    data(i-7).currentData = d;
    data(i-7).si = scale;
    data(i-7).h = hi;
end
current_withF_VandD = data;
save current_withF_VandD

