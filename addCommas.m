function  [str] = addCommas(num)
str=num2str(fix(num));
for i=length(str)-2:-3:2
    str(i+1:end+1) = str (i:end);
    str(i) = ',';
end
