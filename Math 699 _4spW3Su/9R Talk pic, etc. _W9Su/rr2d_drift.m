function [varargout] = rr2d_drift(Nbugs,gridstruct0)
%RR2D   Rotor-router walk in two dimensions.
%   GS = RR2D_DRIFT(NBUGS) performs a rotor-router walk for NBUGS bugs and returns
%   a struct GS with members GS.GRID and GS.NBUGS.  GS.GRID is a square
%   matrix giving the number of times each site was visited, with all bugs
%   starting at the central site.
%
%   GS = RR2D_DRIFT(NBUGS,GS) continues the walk for NBUGS more bugs.
%
%   [GS,BUGPATH] = RR2D(1,...) plots the path BUGPATH of the single bug and
%   optionally returns BUGPATH as rows of grid coordinates.
%
%   References:
%
%   M. Kleber, "Goldbug Variations," The Mathematical Intelligencer,
%   Vol. 27, 55-63 (2005).
%
%   See also RR1D.

%   RR2D_DRIFT(NBUGS, BDRY) returns a rotor-router walk for NBUGS bugs with
%   drift [0 1] and with a boundary plotted. 

if nargin < 1, Nbugs = 1000; end

addpath rr2d

% Size of grid quadrant for a given number of bugs.
gridquadsize = @(Nbugs) ceil(4*sqrt(Nbugs));

if nargin < 2
  % No initial grid was specified, so make an empty one.
  Nbugs0 = 0;
  maxy = gridquadsize(Nbugs); miny = -maxy; 
  maxx = maxy; minx = -maxx;
  x = minx:maxx; y = miny:maxy;
  grid = zeros(length(x), length(y));
  grid0 = [find(x==0) find(y==0)];
else
  error('Restarting nonsquare grid not implemented yet.')

  % An initial grid structure was specified.  Continue adding bugs to it.
  Nbugs0 = gridstruct0.Nbugs;
  grid0 = gridstruct0.grid;
  %grid0 = [find(x==0) find(y==0)];
  %if length(size(grid0)) > 2 || size(grid0,1) ~= size(grid0,2)
    %error('Starting grid must be a square matrix.')
  %end
  Ngrid0 = (size(grid0(1),1) - 1)/2;
  if rem(Ngrid0,1)
    error('Starting grid must have odd number of elements.')
  end
  
  maxy = gridquadsize(Nbugs0+Nbugs); miny = -maxy; 
  maxx = maxy; minx = -maxx;
  x = minx:maxx; y = miny:maxy;
  grid = zeros(length(x), length(y));
  % Embed the starting grid in a larger grid, keeping it centered.
  ii = (1:(2*Ngrid0+1)) + (maxy-Ngrid0); 
  grid(ii,ii) = grid0;
end

if Nbugs == 1
  usemex = false;
  bugpath = [];
else
  usemex = false; % don't use mex file, ever
end

if usemex
  % Use mex file rr2d_mex.
  grid = rr2d_mex(Nbugs,grid);
else
  X0 = grid0;                % center of grid

  for i = 1:Nbugs
    if ~mod(i,1), fprintf('bug %8d\n',i); end

    % Start bug at origin.
    X = [X0 X0];
    while true
      if Nbugs == 1, bugpath = [bugpath ; X]; end
      if ~grid(X(1),X(2))
        % Grid cell is empty.
        % Occupy grid by setting grid to 1.
        % (Another interesting case: set to randi(4); still circular.)
        grid(X(1),X(2)) = 1;
        break;
      else
        % Grid cell has an arrow.
        % Rotate arrow ccw then move in that direction.
        grid(X(1),X(2)) = grid(X(1),X(2)) + 1;
        % For internal DLA, just pick direc = randi(4).
        direc = mod(grid(X(1),X(2))-1,4) + 1;
        switch direc
          case 1 % E
            X(1) = X(1) + 1;
          case 2 % N
            X(2) = X(2) + 1;
          case 3 % W
            X(1) = X(1) - 1;
          case 4  % S
            X(2) = X(2) - 1;
        end
        X(2) = X(2) + 1;
      end
    end
  end
end

gridstruct.grid = grid;
gridstruct.Nbugs = Nbugs0+Nbugs;
gridstruct.x = x;
gridstruct.y = y;

if Nbugs == 1
  % Plot the path of the single bug that was added.
  rr2d_plot(grid,bugpath)
  
  if nargout > 1
    varargout{2} = bugpath;
  end
else
  rr2d_plot(gridstruct)
end

% % plot the contour 2/18 4spW4M
% alist = [.38 .39 .4];
% Alist = 2*Nbugs;
% for i = 1:length(alist)
%   a = alist(i);
%   for j = 1:length(Alist)
%     A = Alist(j);
%     c = (A*sqrt(a)*3/2*sqrt(3/2/pi))^(-1/3);
%     y_max = 1/c^2;
% 
%     n = 500; % how many points? 
%     y2 = linspace(0,y_max,n);
% 
%     x2 = [];
%     for i=1:n
%       x2_each = sqrt(-y2(i)/a*log(c^2*y2(i)));
%       x2 = [x2, x2_each];
%     end
% 
%     hold on
%     plot(x2,y2 + 10,'c--','LineWidth',3) 
%     plot(-x2,y2 + 10,'c--','LineWidth',3) 
%   end
% end
% hold off

if nargout > 0
  varargout{1} = gridstruct;
end
