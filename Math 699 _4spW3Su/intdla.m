function [varargout] = intdla(Npart)
%INTDLA   Simulate internal diffusion-limited aggregation.
%   INTDLA(N) simulates an internal DLA process with N particles.
%   (If unspecified N defaults to 10000.)
%
%   G = INTDLA(N) also returns the final occupation grid.

if nargin < 1
  Npart = 1000;                % number of particles
end

drift = [0 0];
% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(sqrt(Npart));  % size of grid quadrant
grid0 = Ngrid+1;               % m center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

rng('default')

for i = 1:Npart
  X = [0 0];
  while 1
    % Note: these can move diagonally.
    X = X + (randi(3,1,2)-2) + drift;
    if ~grid(X(1)+grid0,X(2)+grid0)
      grid(X(1)+grid0,X(2)+grid0) = 1;
      break
    end
  end
  if ~mod(i,100)
    fprintf('particle %g\n',i)
    plotgrid(x,grid)
    drawnow
  end
end

plotgrid(x,grid)

if nargout > 0
  varargout{1} = grid;
end

%==========================================================================
function plotgrid(x,grid)

imagesc(x,x,grid.'), axis square, axis xy
xlabel('x'), ylabel('y')
