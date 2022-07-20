x = [0:1:100]
c = 1
y = c*x^2;
for n = 1:1:365
    if 0.5 > n > 0.7
        c = 3;
    else 
        c = 0.2;
    end
end
    