%******************************************************************************
%This program is to use HH model to simulate the Mechanosensory nompc channel.
%From the basic of an ion in HH model like g=gmax*m^k*h^n£»
%We propose the form of Mechanical force induced current is like the
%conductance.
% ! Optimizing method is fit !
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


%-------------------------------------------------------
% Ploting the Imax of different cell datas

Im = min(current_Data-current_Data(1,1));
cell(:,1) = fliplr(Im(4:13));
cell(:,2) = fliplr([Im(14:15) Im(17:24)]);
cell(:,3) = fliplr(Im(25:34));
cell(:,4) = fliplr(Im(35:44));
x = reshape([1:10],10,1);

result = zeros(4,2);

figure;
hold on;  
for i=1:4
    
%    f = fittype('a*x+b','independent','x','coefficients',{'a','b'});
    cfun = fit(x,cell(:,i),'poly1');
    result(i,1)=cfun.p1;
    result(i,2)=cfun.p2;
   
end

plot(x,result(1,1)*x+result(1,2),'r');
plot(x,cell(:,1),'r*');
plot(x,result(2,1)*x+result(2,2),'g');
plot(x,cell(:,2),'g*');
plot(x,result(3,1)*x+result(3,2),'b');
plot(x,cell(:,3),'b*');
plot(x,result(4,1)*x+result(4,2),'k');
plot(x,cell(:,4),'k*');

plot(1:13,Im(1:13),'ro');plot(14:24,Im(14:24),'go');plot(25:34,Im(25:34),'bo');plot(35:44,Im(35:44),'ko');
xlabel('cell and force');
ylabel('Imax');


%----------------------------------------------------
% Calculate the mean Imax of 4 different cells.
% Then fit the mean_cell_Imax.

mean_cell = zeros(10,1);
for i=1:10
    mean_cell(i,1) = mean(cell(i,:));
end

[cfun_mean, gof] = fit(x,mean_cell(1:10,1),'poly1');

figure;
hold on;
plot(x,cfun_mean.p1*x+cfun_mean.p2);
plot(x,mean_cell(:,1),'o');
xlabel('Force/um');
ylabel('mean_Imax');



