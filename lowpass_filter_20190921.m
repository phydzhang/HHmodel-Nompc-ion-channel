clear
clc

data1 = abfload('2019_07_20_0202.abf');
data1_lowpass = lowpass(data1(:,1),10,500);

windowSize = 15; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

data1_lowpass_mean = filter(b,a,data1_lowpass);

figure;
hold on;
%plot(data1(:,1));
plot(data1_lowpass-10);
plot(data1_lowpass_mean);


