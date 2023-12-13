% Compare numerical to approximate contours.

A = linspace(10,1010,10);
N = 200;

for i = 1:length(A)
  d = area2contour(A(i));
  [x,y,xa,ya] = contourxy(d,N);
  plot(x,y,'k-','LineWidth',2)
  hold on
  %axis equal

  % Approximate contour.
  plot(xa,ya,'r--','LineWidth',2)
end

hold off
xlabel('x')
ylabel('y')
