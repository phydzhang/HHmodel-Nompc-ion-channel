function output=choose_Imax(i)

Imax_cell = [153.8086 171.5088 122.6807 71.4111];
if (i >= 1) && (i <= 13)
    output = Imax_cell(1);
elseif (i >= 14) && (i <= 24)
    output = Imax_cell(2);
elseif (i >= 25) && (i <= 34)
    output = Imax_cell(3);
else
    output = Imax_cell(4);
end

end