d = logspace(-4,1,100);

figure(1)
A = contour2area(d);
loglog(d,A,'LineWidth',2)
xlabel('d')
ylabel('A')

% Compare to approximate values.
c = sqrt(2/pi)*d;
Aapprox = 2/3*sqrt(2/3*pi)*c.^-3;
hold on
loglog(d,Aapprox,'--','LineWidth',2)
hold off

figure(2)
% Compare relative error of the approximation.
errA = abs(A - Aapprox)./Aapprox;
loglog(d,errA,'LineWidth',2)
xlabel('d')
ylabel('relative error A')

figure(3)
semilogx(d,A./Aapprox,'.-','LineWidth',1)
xlabel('d')
ylabel('A / Aapprox')
