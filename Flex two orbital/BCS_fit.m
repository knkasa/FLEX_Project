function f = BCS_fit(x, xdata)

 f = x(1).*tanh( x(2).*sqrt(1-xdata./x(3) ) )  ;


%f = 0.5*x(1).*(3*xdata.^2-1)+0.125*x(2).*(35*xdata.^4-30*xdata.^2+3)+ ...
 %   0.0625*x(3).*(231*xdata.^6-315*xdata.^4+105*xdata.^2-5);