function [A,B]=load_currentData()

load('./data/currentData_step_T0_200.mat');
A = [10,10,10,10,9:-1:1,10,9,9,8:-1:1,10:-1:2,2,10:-1:1];
B = zeros(length(currentData(1).current),length(currentData));

for i=1:44
    B(:,i) = currentData(i).current;
end


end