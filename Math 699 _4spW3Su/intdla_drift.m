function [varargout] = intdla_drift(Npart)
%INTDLA   Simulate internal diffusion-limited aggregation.
%   INTDLA_DRIFT(N) simulates an internal DLA process with N particles with
%   drift [0 1]
%   (If unspecified N defaults to 10000.)
%
%   G = INTDLA_DRIFT(N) also returns the final occupation grid.

if nargin < 1
  Npart = 1000; % number of particles 
end
  
biasedwalk = false;

if biasedwalk
  % Instead of making a "drift" jump at each timestep, bias the walk.
  drift = [0 0]; %#ok<UNRCH>
  w1 = .63;  % probability of moving north (to match "jump" case)
  w = [w1 ones(1,3)*(1-w1)/3];

  if any(w) < 0 || abs(sum(w)-1) > 10*eps
    error('Bad weights vector [%s ].',sprintf(' %g',w));
  end

  % Same formula as for "jump" case, but with Npart/2.
  maxy = ceil(1.28164*(Npart/2)^.66 + 4.96289); miny = -5;
  maxx = 3*ceil(Npart/maxy); minx = -maxx;
else
  % Move with equal probabilities, but take an extra "drift" jump.
  drift = [0 1];
  w = ones(1,4)/4;

  % Grid boundaries with drift.
  % Formulas are empirical.
  maxy = ceil(4*sqrt(Npart)); miny = -maxy;
  % Deduce maxx by area preservation.
  maxx = maxy; minx = -maxx;
end

x = minx:maxx; y = miny:maxy;
grid = zeros(length(x),length(y));
grid0 = [find(x == 0) find(y == 0)];

rng('default')

v_dir = [0 1; 1 0; 0 -1; -1 0;];

if length(w) ~= size(v_dir,1)
  error('Weights vector has length %d rather than %d.',length(w),size(v_dir,1))
end

for i = 1:Npart
  X = [0 0];
  while 1

    % Select direction with weights vector w.
    d = randsample(size(v_dir,1),1,true,w);

    X = X + v_dir(d,:) + drift;

    if ~grid(X(1)+grid0(1),X(2)+grid0(2))
      grid(X(1)+grid0(1),X(2)+grid0(2)) = 1;
      break
    end
  end
  if ~mod(i,100)
    fprintf('particle %g\n',i)
    plotgrid(x,y,grid)
    drawnow
  end
end

% % plot the contour 2/11 4spW4M
% a = .39;
% if biasedwalk
%   A = Npart;
% else
%   A = 2*Npart;
% end
% c = (A*sqrt(a)*3/2*sqrt(3/2/pi))^(-1/3);
% y_max = 1/c^2;
% 
% n = 500; % how many points? 
% y2 = linspace(0,y_max,n);
% 
% x2 = [];
% for i=1:n
%     x2_each = sqrt(-y2(i)/a*log(c^2*y2(i)));
%     x2 = [x2, x2_each];
% end
% 
% 
% hold on
% plot(x2,y2,'c--','LineWidth',3) % y2 is multipled by 1.2
% plot(-x2,y2,'c--','LineWidth',3) 
% hold off


sz = size(grid); size_y = sz(2);
szx = sz(1);
%grid

for j = 1:(size_y-1)
    for i = 1:szx
        if grid(i,j)==1 & grid(:,j+1)==0
            maxy_rn = j;
        end
    end
end

% maxy_rn
% diff = size_y - maxy_rn

if nargout > 0
  g.grid = grid;
  g.Nbugs = Npart;
  varargout{1} = g;
end

%==========================================================================
function plotgrid(x,y,grid)

imagesc(x,y,grid.'), axis square, axis xy
xlabel('x'), ylabel('y')
