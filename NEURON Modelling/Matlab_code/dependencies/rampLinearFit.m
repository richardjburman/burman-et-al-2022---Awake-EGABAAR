function r = rampLinearFit(vec1V,vec1I)


%Calculating the reversal by fitting a linear polynomial

p1 = polyfit(vec1V,vec1I,1);

r = roots(p1);
