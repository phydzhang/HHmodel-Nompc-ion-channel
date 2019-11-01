clear
clc

distance = zeros(44);
current_Data = zeros(4982,44);
[distance,current_Data] = load_currentData();

%data1 = abfload('2019_07_20_0201.abf');
%L = length(data1);
%fft_data1_1 = fft(data1(:,1));
%P2 = abs(fft_data1_1/L);
%P1 = P2(1:L/2+1);
%P1(2:end-1) = 2*P1(2:end-1);
%f = (0:L/2)*2/3;

figure;hold on;
for i=35:44
    
    data1 = current_Data(200:451,i);
    L = length(data1);
    fft_data1_1 = fft(data1(:,1));
    P2 = abs(fft_data1_1/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (0:L/2)*2/3;
    %plot(f(3:50),P1(3:50));

end
%xlabel('cell3_f(Hz)');
%ylabel('|p1(f)|');


%------------------------------------------------------------------
% See Imax Data(µÝ½øµÄforce´Ì¼¤)

Imax_data = abfload('./data/Imax/2019_05_14_0115.abf');

plot(Imax_data(2000:3800,1,1));
plot(Imax_data(2000:3800,2,1));
plot(Imax_data(2000:3800,2,2));


Imax_data2 = abfload('./data/Imax/2019_05_14_0116.abf');
%plot(Imax_data2(3000:6500,1,1));
%plot(Imax_data2(3000:6500,2,1));
%plot(Imax_data2(3000:6500,2,2));
%plot(Imax_data2(3000:6500,2,3));
%plot(Imax_data2(3000:6500,2,4));



