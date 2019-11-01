%% ******************************************************************************
%This program is to use HH model to simulate the Mechanosensory nompc channel.
%From the basic of an ion in HH model like g=gmax*m^k*h^n£ª
%We propose the form of Mechanical force induced current is like the conductance.
% ! Optimizing !
% ! To simulate the triangle force and response £°
%force_start=[623 623 623 623 623 623 623  623  623]
%force_end  =[630 640 650 674 721 873 1115 3090 5625]
%******************************************************************************
clear
clc

%% Load the current Data caused by force with different V&D
load('./data/current_withF_VandD.mat');
%  Moving the current Data to the same resting value
sameRestingCurrent = zeros(7500,3,9);  % 2 means raw data and smoothing data, 9 means 9 different force type
for i=1:9
    mean1 = mean(current_withF_VandD(i).currentData(500:800,1,1));
    mean2 = mean(current_withF_VandD(i).currentData(500:800,2,1));
    sameRestingCurrent(:,1,i) = current_withF_VandD(i).currentData(:,1,1)-mean1;
    sameRestingCurrent(:,2,i) = current_withF_VandD(i).currentData(:,2,1)-mean2;
    sameRestingCurrent(:,3,i) = current_withF_VandD(i).currentData(:,3,1);
end
% Intercept 501 to 7500 data of the sameRestingCurrent
InterceptData = zeros(7000,3,9);
for i=1:9
    InterceptData(:,:,i)=sameRestingCurrent(501:7500,:,i);
end

%% --------------------------------------------------------------------
% Setting the constant variables
k    = 4;
n    = 4;
m0   = 0;
h0   = 1;
tm   = 0.002;
th   = 0.05;
Imax = 200;

%% ----------------------------------------------------------------------
%Definiting the timeline t & other variations.
t = 0.0002:2e-4:1.4;
force_start=[623 623 623 623 623 623 623  623  623];
force_end  =[630 640 650 674 721 873 1115 3090 5625];
result_km = zeros(9,100);
result_kh = zeros(9,100);
result_betam = zeros(9,100);
result_betah = zeros(9,100);
result_gof = zeros(9,5,100);

%% ----------------------------------------------------------------------
% ! Defining the model !
% ±È¿˙Imax¥”5µΩ500
%for num_Imax=1:100

for num_Imax=1:1
    Imax = 200;
    for j=1:9
        start_force = force_start(j);
        Nt_force = force_end(j)+100;
        K_force  = 10/((Nt_force-start_force)*2e-4);
        force = 0.*(t<t(start_force))+K_force*(t-0.0002*start_force).*(t>=t(start_force)&t<=t(Nt_force))+0.*(t>t(Nt_force));
        t_former = zeros(1,start_force-1);
        t_mid = 0.0002:2e-4:0.0002*(Nt_force-start_force+1);
        t_latter = zeros(1,7000-Nt_force);
        t_correct = [t_former t_mid t_latter];
        
        f = fittype(@(km,betam,kh,betah,force,t)-Imax.*((exp(km*force.^2)./(betam+exp(km*force.^2)))...
            -((exp(km*force.^2)./(betam+exp(km*force.^2)))-m0).*exp(-t/tm)).^k.*((1-exp(kh*force.^2)./...
            (betah+exp(kh*force.^2)))-((1-exp(kh*force.^2)./(betah+exp(kh*force.^2)))-h0).*exp(-t/th)).^n,...
            'problem','force','independent','t','coefficients',{'km','betam','kh','betah'});
        
        force_reshape = reshape(force,7000,1);
        t_reshape = reshape(t_correct,7000,1);
        force_problem = force_reshape(start_force:Nt_force);
        t_active = t_reshape(start_force:Nt_force);
        [cfun,gof] = fit(t_active,InterceptData(start_force:Nt_force,2,j),f,'Lower',[0,0,0,0],'Upper',...
            [1,100,1,100],'StartPoint',[0.5,5,0.02,10],'problem',force_problem);
        
        result_km(j,num_Imax) = cfun.km;
        result_kh(j,num_Imax) = cfun.kh;
        result_betam(j,num_Imax) = cfun.betam;
        result_betah(j,num_Imax) = cfun.betah;
        result_gof(j,1,num_Imax) = gof.sse;
        result_gof(j,2,num_Imax) = gof.rsquare;
        result_gof(j,3,num_Imax) = gof.dfe;
        result_gof(j,4,num_Imax) = gof.adjrsquare;
        result_gof(j,5,num_Imax) = gof.rmse;
        
        mgig = (exp(cfun.km*force.^2)./(cfun.betam+exp(cfun.km*force.^2)));
        hgig = (1-exp(cfun.kh*force.^2)./(cfun.betah+exp(cfun.kh*force.^2)));
        m = mgig-(mgig-m0).*exp(-t_correct/tm);
        h = hgig-(hgig-h0).*exp(-t_correct/th);
        I = -Imax*(m.^k).*(h.^n);
        
%        plot(t,I,'LineWidth',3);
%        plot(t,force+5,'LineWidth',2);
%        plot(t,InterceptData(:,2,j),'LineWidth',1.5);
%        xlabel('Time(S)');
%        ylabel('Current(pA)');
%        box on;
%            force = 0:0.001:10;
%            mgig = (exp(cfun.km*force.^2)./(cfun.betam+exp(cfun.km*force.^2)));
%            hgig = (1-exp(cfun.kh*force.^2)./(cfun.betah+exp(cfun.kh*force.^2)));
%            plot(force,mgig);
%            plot(force,hgig);
    end
end
%% ------------------------------------------------
% Plot the simulated curve
km = mean(result_km(5:9));
kh = mean(result_kh(5:9));
betam = mean(result_betam(5:9));
betah = mean(result_betah(5:9));
for nu=1:9
    start_force = force_start(nu);
    Nt_force = force_end(nu)+100;
    K_force  = 10/((Nt_force-start_force)*2e-4);
    force = 0.*(t<t(start_force))+K_force*(t-0.0002*start_force).*(t>t(start_force)&t<t(Nt_force))+0.*(t>t(Nt_force));

    mgig = (exp(km.*force.^2)./(betam+exp(km.*force.^2)));
    hgig = (1-exp(kh.*force.^2)./(betah+exp(kh.*force.^2)));
    m = mgig-(mgig-m0).*exp(-t/tm);
    h = hgig-(hgig-h0).*exp(-t/th);
    I = -Imax*(m.^k).*(h.^n);

%    figure;hold on;
%    plot(t,I);
%    plot(t,force);
%    plot(t,InterceptData(:,2,nu));
end

figure;
plot(result_km(2:9,1),'linewidth',3);box on;xlabel('Different Velocity');ylabel('km');
figure;
plot(result_kh(2:9,1),'linewidth',3);box on;xlabel('Different Velocity');ylabel('kh');
figure;
plot(result_betam(2:9,1),'linewidth',3);box on;xlabel('Different Velocity');ylabel('betam');
figure;
plot(result_betah(2:9,1),'linewidth',3);box on;xlabel('Different Velocity');ylabel('betah');













