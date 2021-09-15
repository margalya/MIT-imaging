function str = vec2str( vec)
d = diff(vec);
if min(d) == max(d)
    str = [ num2str(min(vec)) ' : ' num2str(d(1)) ' : ' num2str(max(vec))];
else
    str = num2str(vec);
end
end