d = linspace(1e-2,10,100);

[ymin,ymax] = ybounds(d);

figure(1)
loglog(d,-ymin,'r-',d,ymax,'b-','LineWidth',2)
xlabel('d')
ylabel('-y_{min}, y_{max}')

% Compare to approximate values.
c = sqrt(2/pi)*d;
yminapprox = lambertw(4./c.^2)/4;
% First few terms of a small-c expansion.
R = [1 -1/4 3/32 -5/64 231/2048];
cc = [c.^-2;c.^0;c.^2;c.^4;c.^6];
ymaxapprox = R*cc;
hold on
loglog(d,yminapprox,'r--',d,ymaxapprox,'b--','LineWidth',2)
hold off
legend('-y_{min}','y_{max}')

figure(2)
errymin = (yminapprox + ymin)./-ymin;
errymax = (ymaxapprox - ymax)./ymax;
loglog(d,errymin,'r-',d,errymax,'b-','LineWidth',2)
xlabel('d')
ylabel('error y_{min}, y_{max}')
legend('error y_{min}','error y_{max}','Location','NorthWest')
