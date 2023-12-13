function d = area2contour(A)
%AREA2CONTOUR   Contour value from enclosed area for source with drift.
%   D = AREA2CONTOUR(A) converts enclosed area A to contour value D
%   for the function
%
%     exp(Y) K0(R) == D
%
%   where K0 is the 0th modified Bessel function of the second kind, and
%   R = sqrt(X^2 + Y^2).
%
%   See also YBOUNDS, CONTOUR2AREA, CONTOURXY.

if nargin < 1, A = 1; end

if ~isscalar(A)
  % Lazy vectorize if A is a vector.
  d = zeros(size(A));
  for i = 1:length(d)
    d(i) = area2contour(A(i));
  end
  return
end

cguess = (2/3*sqrt(2/3*pi)/A).^(1/3);
dguess = sqrt(pi/2)*cguess;

f = @(d) contour2area(d) - A;
d = fzero(f,dguess);
