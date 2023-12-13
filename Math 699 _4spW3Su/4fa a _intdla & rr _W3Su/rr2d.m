function [varargout] = rr2d(Nbugs,gridstruct0)
%RR2D   Rotor-router walk in two dimensions.
%   GS = RR2D(NBUGS) performs a rotor-router walk for NBUGS bugs and returns
%   a struct GS with members GS.GRID and GS.NBUGS.  GS.GRID is a square
%   matrix giving the number of times each site was visited, with all bugs
%   starting at the central site.
%
%   GS = RR2D(NBUGS,GS) continues the walk for NBUGS more bugs.
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

if nargin < 1, Nbugs = 1000; end

% Size of grid quadrant for a given number of bugs.
gridquadsize = @(Nbugs) ceil(.6*sqrt(Nbugs));

if nargin < 2
  % No initial grid was specified, so make an empty one.
  Nbugs0 = 0;
  Ngrid = gridquadsize(Nbugs);
  grid = zeros(2*Ngrid+1);
else
  % An initial grid structure was specified.  Continue adding bugs to it.
  Nbugs0 = gridstruct0.Nbugs;
  grid0 = gridstruct0.grid;
  if length(size(grid0)) > 2 || size(grid0,1) ~= size(grid0,2)
    error('Starting grid must be a square matrix.')
  end
  Ngrid0 = (size(grid0,1) - 1)/2;
  if rem(Ngrid0,1)
    error('Starting grid must have odd number of elements.')
  end
  Ngrid = gridquadsize(Nbugs0+Nbugs);
  % Embed the starting grid in a larger grid, keeping it centered.
  grid = zeros(2*Ngrid+1);
  ii = (1:(2*Ngrid0+1)) + (Ngrid-Ngrid0);
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
  X0 = Ngrid+1;                % center of grid

  for i = 1:Nbugs
    if ~mod(i,1000), fprintf('bug %8d\n',i); end

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
      end
    end
  end
end

gridstruct.grid = grid;
gridstruct.Nbugs = Nbugs0+Nbugs;

if Nbugs == 1
  % Plot the path of the single bug that was added.
  rr2d_plot(grid,bugpath)
  if nargout > 1
    varargout{2} = bugpath;
  end
else
  rr2d_plot(gridstruct)
end

if nargout > 0
  varargout{1} = gridstruct;
end
