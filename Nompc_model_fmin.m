%******************************************************************************
%This program is to use HH model to simulate the Mechanosensory nompc channel.
%From the basic of an ion in HH model like g=gmax*m^k*h^n£»
%We propose the form of Mechanical force induced current is like the
%conductance.
% ! Optimizing method is fmincon !
%******************************************************************************

clear
clc

%-------------------------------------------------------------------
% Loading the experimental data.
distance = zeros(44);
current_Data = zeros(4982,44);
[distance,current_Data] = load_currentData();

%-------------------------------------------------------------------
% Moving the current_Data to the same starting line
for i=1:44
    move_d = current_Data(1,i)-current_Data(1,1);
    current_Data(:,i) = current_Data(:,i)-move_d;
end

%--------------------------------------------------------------------
% Setting the constant variables
k    = 2;
n    = 3;
m0   = 0;
h0   = 1;
Imax = 171.5088;

%----------------------------------------------------------------------
%Definiting the timeline t & force change over time
t = 0:2e-4:0.4582;
%force = 0.*(t<t(200))+2*t.*(t>t(200)&t<t(2000))+0.*(t>t(2000));
%force = 0.*(t<t(200))+10.*(t>t(200));

% Variables Changing With Force
%tm   = x(3)  0.12;
%th   = x(4)  0.5;
%Ib   = x(5)  -1032.7;
%Ib   = -1032.7;

%-----------------------------------------------------------
% Parameter Optimizing

A=[];
b=[];
Aeq=[];
beq=[];
%lb=[0,0,0,0,-1500];
%ub=[1,1,1,1,0];
lb=[0,0,0,0];
ub=[1,1,1,1];

result = zeros(44,5);

for i=1:44
        
    x0 = [0.9,0.5,0.001,0.5];
    [x,fval] = fmincon(@(x)Rmse(x,m0,h0,k,n,Imax,t,current_Data,i),x0,A,b,Aeq,beq,lb,ub);
    result(i,1)=x(1);
    result(i,2)=x(2);
    result(i,3)=x(3);
    result(i,4)=x(4);
%    result(i,5)=x(5);
    result(i,5)=fval;
    
end

%--------------------------------------------------------------------------
% Finding the mean value of optimized parameter from different cell datas

mean_result = zeros(10,5);
num = zeros(10);

for i=1:44
    mean_result(distance(i),:) = mean_result(distance(i),:)+result(i,:);
    num(distance(i)) = num(distance(i))+1;
end

for i=1:10
    mean_result(i,:)=mean_result(i,:)/num(i);
end

figure;
hold on;

for i=1:10
    tm = mean_result(i,3);
    th = mean_result(i,4);
    mgig = mean_result(i,1);
    hgig = mean_result(i,2);
    
    plot(t,Isimulate(mgig,hgig,m0,h0,k,n,tm,th,Imax,t));
end
%plot(t,current_Data(200:2491,4:13))

plot(reshape(t,2292,1),current_Data(200:2491,4:13));




