%%******************************************************************************
%This program is to use HH model to simulate the Mechanosensory nompc channel.
%From the basic of an ion in HH model like g=gmax*m^k*h^n£»
%We propose the form of Mechanical force induced current is like the
%conductance.
% ! Optimizing method is fit !
% ! Different cell with different Imax !
%******************************************************************************

clear
clc

%% -------------------------------------------------------------------
% Loading the experimental data.
distance = zeros(44);
current_Data = zeros(4982,44);
[distance,current_Data] = load_currentData();

%% -------------------------------------------------------------------
% Moving the current_Data to the same starting line
for i=1:44
    move_d = current_Data(1,i)-current_Data(1,1);
    current_Data(:,i) = current_Data(:,i)-move_d;
end

%% --------------------------------------------------------------------
% Setting the constant variables
k    = 2;
n    = 3;
m0   = 0;
h0   = 1;
Imax = 171.5088;

%% ----------------------------------------------------------------------
%Definiting the timeline t & force change over time
%t = 0:2e-4:0.9962;
t = 0:2e-4:0.4582;
%force = 0.*(t<t(200))+2*t.*(t>t(200)&t<t(2000))+0.*(t>t(2000));
%force = 0.*(t<t(200))+10.*(t>t(200));

%% -------------------------------------------------------------
% Variables Changing With Force
tm   = 0.002;
th   = 0.05;
%Ib   = x(5)  -1032.7;
%Ib   = -1032.7;

%% -----------------------------------------------------------
% Parameter Optimizing
result = zeros(44,2);
error  = zeros(44,3);

for i=1:44

    Imax = choose_Imax(i);
    f = fittype(@(mgig,hgig,t)-Imax.*(mgig-(mgig-m0).*exp(-t/tm)).^k.*(hgig-(hgig-h0).*exp(-t/th)).^n-1032,'independent','t',...
        'coefficients',{'mgig','hgig'});
    
    
    t_reshape = reshape(t,2292,1);
    [cfun,gof] = fit(t_reshape,current_Data(200:2491,i),f,'Lower',[0,0],'Upper',[1,1],'StartPoint',[0.9,0.1]);
    
    result(i,1)=cfun.mgig;
    result(i,2)=cfun.hgig;
    error(i,1) =gof.rsquare;
    error(i,2) =gof.dfe;
    error(i,3) =gof.rmse;
end


%% --------------------------------------------------------------------------
% Finding the mean value of optimized parameter from different cell datas

mean_result = zeros(10,2);
num = zeros(10);

for i=1:44
    mean_result(distance(i),:) = mean_result(distance(i),:)+result(i,:);
    num(distance(i)) = num(distance(i))+1;
end

for i=1:10
    mean_result(i,:)=mean_result(i,:)/num(i);
end

%% -----------------------------------------------------
% Calculate the rmse with dirrerent parameter(mean or individual)
para_rmse       = zeros(44,1);
mean_para_rmse  = zeros(44,1);
for i=1:44
    Imax = choose_Imax(i);
    mgig=result(i,1);
    hgig=result(i,2);
    para_rmse(i,1)=calRmse(mgig,hgig,m0,h0,tm,th,k,n,Imax,t,current_Data,i); 
end    

for i=1:44
    Imax = choose_Imax(i);
    mgig=mean_result(distance(i),1);
    hgig=result(distance(i),2);
    mean_para_rmse(i,1)=calRmse(mgig,hgig,m0,h0,tm,th,k,n,Imax,t,current_Data,i); 
end    

%figure;hold on;
%plot(para_rmse(:,1));
%plot(mean_para_rmse(:,1));



%% ------------------------------------------------
% Plot the simulated curve

%figure;hold on;
for i=4:13
    Imax = choose_Imax(i);
    mgig = result(i,1);
    hgig = result(i,2);
    
    %plot(t,Isimulate(mgig,hgig,m0,h0,k,n,tm,th,Imax,t));
end
%plot(t,current_Data(200:2491,4:13))


%plot(t_reshape,current_Data(200:2491,4:13));


%% ------------------------------------------------------------------------
% Notice that there are same flutuations in each current curve !

figure;hold on;
%t_new = 0:2e-4:0.08;
%plot(reshape(t_new,401,1),current_Data(200:600,4:13));
%xlabel('cell4');

t_new = 0:2e-4:0.05;
for i=1:44
    if(distance(i)==8)
        plot(reshape(t_new,251,1),current_Data(200:450,i));
    end
end

%xlabel('Force=8um');



