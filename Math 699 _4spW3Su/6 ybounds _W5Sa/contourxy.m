function [x,y,xa,ya] = contourxy(d,N)
%CONTOURXY   X and Y coordinates of a contour for source with drift.
%   [X,Y] = CONTOURXY(D) returns two vectors of coordinates X and Y
%   giving contours for the function
%
%     exp(Y) K0(R) == D
%
%   where K0 is the 0th modified Bessel function of the second kind, and
%   R = sqrt(X^2 + Y^2).
%
%   [X,Y,XA,YA] = CONTOURXY(D) also returns an approximation to the
%   contour XA,YA.
%
%   CONTOURXY(D,N) uses 2N-1 points to discretize the contour (default N=100).
%
%   See also YBOUNDS, AREA2CONTOUR, CONTOUR2AREA.

if nargin < 1, d = 1; end
if nargin < 2, N = 100; end

% Find bounds.
[ymin,ymax] = ybounds(d);
rmin = -ymin; rmax = ymax;
r = linspace(rmin,rmax,N);

% Contour.
y = log(d./besselk(0,r));
x = real(sqrt(r.^2 - y.^2)); % first and last points can have small imag part

% Reflect contour about y axis.
i = [N-1:-1:1];
x = [x -x(i)];
y = [y y(i)];

if nargout > 2
  % Compute approximation to contour.
  c = sqrt(2/pi)*d;
  ya = linspace(0,1/c^2,N);
  xa = sqrt(-ya.*log(c^2*ya));
  xa(1) = 0;
  xa = [xa -xa(i)];
  ya = [ya ya(i)];
end
