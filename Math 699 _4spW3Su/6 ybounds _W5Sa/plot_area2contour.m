A = logspace(-2,6,30);

loglog(A,area2contour(A),'LineWidth',2)
xlabel('A')
ylabel('d')

% Compare to approximate values.
d = sqrt(pi/2)*(2/3*sqrt(2/3*pi)./A).^(1/3);
hold on
loglog(A,d,'--','LineWidth',2)
hold off
