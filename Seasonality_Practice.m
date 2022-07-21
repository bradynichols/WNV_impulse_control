x = [0:1:365]
c = 1
y = c*x;
for n = 1:1:365
    if (0 < n < 50)
        c = 3
        plot(x,y)
        hold on;
    elseif (51 < n < 200)
        c = 0.2
        plot(x,y)
        hold on;
    else
        c = 10;
        plot(x,y)
        hold on;
    end
end
    