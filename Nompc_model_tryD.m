%% ******************************************************************************
%This program is to use HH model to simulate the Mechanosensory nompc channel.
%From the basic of an ion in HH model like g=gmax*m^k*h^n；
%We propose the form of Mechanical force induced current is like the conductance.
% ! Optimizing !
%******************************************************************************
clear
clc

%% -------------------------------------------------------------------
% Loading the experimental data.
distance = zeros(44);
current_Data = zeros(4982,44);
[distance,current_Data] = load_currentData();

% Moving the current_Data to the same starting line
for i=1:44
    move_d = mean(current_Data(1:100,i));
    current_Data(:,i) = current_Data(:,i)-move_d;
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
t = 0:2e-4:0.4582;
result_km = zeros(44,1);
result_kh = zeros(44,1);
result_betam = zeros(44,1);
result_betah = zeros(44,1);

%% ----------------------------------------------------------------------
% ! Defining the model !
for j=1:44
    force = 0.*(t<t(200))+distance(j).*(t>=t(200));
    t_former = zeros(1,199);
    t_mid = 0.0002:2e-4:0.0002*2093;
    t_correct = [t_former t_mid];
    
    f = fittype(@(km,betam,kh,betah,force,t)-Imax.*((exp(km*force.^2)./(betam+exp(km*force.^2)))...
        -((exp(km*force.^2)./(betam+exp(km*force.^2)))-m0).*exp(-t/tm)).^k.*((1-exp(kh*force.^2)./...
        (betah+exp(kh*force.^2)))-((1-exp(kh*force.^2)./(betah+exp(kh*force.^2)))-h0).*exp(-t/th)).^n,...
        'problem','force','independent','t','coefficients',{'km','betam','kh','betah'});

%    f = fittype(@(km,betam,hg,force,t)-Imax.*((exp(km*force.^2)./(betam+exp(km*force.^2)))...
%        -((exp(km*force.^2)./(betam+exp(km*force.^2)))-m0).*exp(-t/tm)).^k.*(hg-(hg-h0).*exp(-t/th)).^n,...
%        'problem','force','independent','t','coefficients',{'km','betam','hg'});
    
    force_reshape = reshape(force,length(force),1);
    t_reshape = reshape(t_correct,length(t_correct),1);
    [cfun,gof] = fit(t_reshape(1:2292,1),current_Data(1:2292,j),f,'Lower',[0,0,0,0],'Upper',...
        [1,100,1,100],'StartPoint',[0.5,5,0.5,10],'problem',force_reshape);
    
    result_km(j) = cfun.km;
    result_kh(j) = cfun.kh;
    result_betam(j) = cfun.betam;
    result_betah(j) = cfun.betah;
    
    mgig = (exp(cfun.km*force.^2)./(cfun.betam+exp(cfun.km*force.^2)));
    hgig = (1-exp(cfun.kh*force.^2)./(cfun.betah+exp(cfun.kh*force.^2)));

    m = mgig-(mgig-m0).*exp(-t_correct/tm);
    h = hgig-(hgig-h0).*exp(-t_correct/th);
    I = -Imax*(m.^k).*(h.^n);
%    figure;hold on;
%    plot(t,I);
%    plot(t,force);
%    plot(t,current_Data(1:2292,j));
%    xlabel('Time(S)');
%    ylabel('Current(pA)');
%    box on;
end
%plot(flipud(result_betah(38:44,1)),'Linewidth',3);
%xlabel('Distance(μm)');
%ylabel('βh');


%% ------------------------------------------------
% Plot the simulated curve
km = mean(result_km);
kh = mean(result_kh);
betam = mean(result_betam);
betah = mean(result_betah);

%force = 0.*(t<t(200))+20*t.*(t>=t(200)&t<t(1000))+(3.6+40*(t-t(1000))).*(t>=t(1000));
force = 0.*(t<t(200))+5.*(t>=t(200)&t<t(1000))+10.*(t>=t(1000));
t_former = zeros(1,199);
t_mid = 0.0002:2e-4:0.0002*801;
t_latter = 0.0002:2e-4:0.0002*1292;
t_correct = [t_former t_mid t_latter];
mgig = (exp(cfun.km*force.^2)./(cfun.betam+exp(cfun.km*force.^2)));
hgig = (1-exp(cfun.kh*force.^2)./(cfun.betah+exp(cfun.kh*force.^2)));
m = mgig-(mgig-m0).*exp(-t_correct/tm);
h = hgig-(hgig-h0).*exp(-t_correct/th);
I = -Imax*(m.^k).*(h.^n);
%figure;hold on;
%plot(t,I);
%plot(t,force);
%plot(t,m);


%----------------------------------------------------------------------
% 画出m无穷和h无穷关于force的关系
%force = 0:0.001:10;
%km = 0.3;
%betam = 80;
%kh = 0.08;
%betah = 50;
%mgig = (exp(km*force.^2)./(betam+exp(km*force.^2)));
%hgig = (1-exp(kh*force.^2)./(betah+exp(kh*force.^2)));
%figure;hold on;
%plot(force,mgig);
%plot(force,hgig,'--');
%xlabel('D(μm)');
%ylabel('Open Probability');
%box on;
