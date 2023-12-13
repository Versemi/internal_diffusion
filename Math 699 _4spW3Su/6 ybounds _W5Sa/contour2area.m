function A = coutour2area(d)
%CONTOUR2AREA   Enclosed area of countour for a source with drift
%   A = CONTOUR2AREA(D) converts a contour value D to the enclosed area A
%   for the function
%
%     exp(Y) K0(R) == D
%
%   where K0 is the 0th modified Bessel function of the second kind, and
%   R = sqrt(X^2 + Y^2).
%
%   See also YBOUNDS, AREA2CONTOUR, CONTOURXY.

if nargin < 1, d = 1; end

A = zeros(size(d));

for i = 1:length(d)
  % Find bounds.
  [ymin,ymax] = ybounds(d(i));
  rmin = -ymin; rmax = ymax;

  if d(i) >= .048
    y0 = @(r) log(d(i)./besselk(0,r));
    dydr0 = @(r) besselk(1,r)./besselk(0,r);
    y1 = y0;
    dydr1 = dydr0;
  else
    y0 = @(r) log(d(i)./besselk(0,r));
    dydr0 = @(r) besselk(1,r)./besselk(0,r);
    % For large r, use a large-r expansion of K0 to avoid underflow.
    y1 = @(r) r + log(sqrt(2/pi)*d(i).*sqrt(r)) - ...
         log(1 - 1/8./r + 9/128./r.^2 - 75/1024./r.^3);
    dydr1 = @(r) (525 - 210*r + 240*r.^2 - 768*r.^3 - 2048*r.^4) ./ ...
            (150*r - 144*r.^2 + 256*r.^3 - 2048*r.^4);
  end

  yminerr = abs(y0(rmin) - ymin)/-ymin;
  if yminerr > 1e-8
    warning('yminerr = %.4e for d = %g',yminerr,d(i))
  end

  ymaxerr = abs(y1(rmax) - ymax)/ymax;
  if ymaxerr > 1e-8
    warning('ymaxerr = %.4e for d = %g',ymaxerr,d(i))
  end

  rswitch = min(600,rmax);

  % Break up integral into [rmin,rswitch] and [rswitch,rmax].
  x0 = @(r) sqrt(r.^2 - y0(r).^2);
  dAdr0 = @(r) x0(r).*dydr0(r);
  A0 = 2*integral(dAdr0,rmin,rswitch);

  if rswitch < rmax
    x1 = @(r) sqrt(r.^2 - y1(r).^2);
    dAdr1 = @(r) x1(r).*dydr1(r);
    A1 = 2*integral(dAdr1,rswitch,rmax);
  else
    A1 = 0;
  end

  A(i) = A0 + A1;
end
