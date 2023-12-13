function [ymin,ymax] = ybounds(d)
%YBOUNDS   Min and max Y values of a countour for a source with drift.
%   [YMIN,YMAX] = YBOUNDS(D) returns the minimum and maximum Y values for
%   the contour
%
%     exp(Y) K0(R) == D
%
%   where K0 is the 0th modified Bessel function of the second kind, and
%   R = sqrt(X^2 + Y^2).
%
%   See also CONTOUR2AREA, AREA2CONTOUR, CONTOURXY.

if nargin < 1, d = 1; end

if ~isscalar(d)
  % Lazy vectorize if d is a vector.
  ymin = zeros(size(d)); ymax = zeros(size(d));
  for i = 1:length(d)
    [ymin(i),ymax(i)] = ybounds(d(i));
  end
  return
end

% Use approximate values as initial guesses.
if d < 5
  % Relate to the earlier constant c, assuming a = 1:
  c = sqrt(2/pi)*d;
  yminguess = lambertw(4/c^2)/4;
  if d < 2
    % When d is small enough, use as good an approximation as possible.
    % First few terms of a small-c expansion.
    % This is enough to get machine precision relative error when
    % d < .048, as required below.
    R = [1 -1/4 3/32 -5/64 231/2048];
    cc = c.^(2*(1:length(R))-4);
    ymaxguess = cc*R';
  else
    ymaxguess = 1/c^2;
  end
else
  % For large d the contours are nearly circular, as for a point source.
  % The contours satisfy log r == -d.
  ymaxguess = exp(-d);
  yminguess = ymaxguess;
end

% This fails if d is too large (> 5).
f1 = @(r) log(d./besselk(0,r)) + r;
ymin = fzero(f1, yminguess);

if d >= .048
  % This fails if d is too small (< .048) or too large (> 5).
  f2 = @(r) log(d./besselk(0,r)) - r;
  ymax = fzero(f2, ymaxguess);
else
  % Use approximate value if d is too small.
  % The relative accuracy is machine precision.
  ymax = ymaxguess;
end

% ymin should be negative.
ymin = -ymin;

end
