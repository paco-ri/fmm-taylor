cols = 1:4;
plot(intquadtest(1,cols).*intquadtest(2,cols).*intquadtest(3,cols), ...
    intquadtest(4,cols), 'o-')
hold on
for i = 1:0
    cols = cols + 3;
    plot(intquadtest(1,cols).*intquadtest(2,cols).*intquadtest(3,cols), ...
        intquadtest(4,cols), 'o-')
end