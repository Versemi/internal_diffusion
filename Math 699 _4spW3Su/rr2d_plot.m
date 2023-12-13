function rr2d_plot(gridstruct,bugpath)
%RR2D_PLOT   Plot RR2D grid.
%
%   See also RR2D, RR2D_ARROWS.

if isstruct(gridstruct)
  grid = gridstruct.grid;
else
  grid = gridstruct;
end
Ngrid = (size(grid,1)-1)/2;
x = -Ngrid:Ngrid;

grid = rr2d_arrows(grid);

imagesc(x,x,grid.'), axis square, axis xy
xlabel('x'), ylabel('y')
axis off

%map = [
%    0.0 0.0 0.0  % black (empty site)
%    1.0 0.0 0.0  % red
%    1.0 1.0 0.0  % yellow
%    0.0 1.0 1.0  % cyan
%    0.0 0.0 1.0  % blue
%      ];

% To see the same patterns as in Kleber's fig. 5, need to choose colors
% carefully.  The colors for 1/3 (E/W) and 2/4 (N/S) must be close to
% each other, but different enough.
% map = [
%     0.0 0.0 0.0
%     0.1 0.1 0.1 % E c1
%     0.7 0.7 0.7 % N c2
%     0.4 0.4 0.4 % W color close to c1
%     1.0 1.0 1.0 % S color close to c2
%       ];
  
map = [
0.0 0.0 0.0
0.4 0.4 0.4 
0.4 0.4 0.4 
0.4 0.4 0.4 
0.4 0.4 0.4 
  ];

colmap = map;
%colmap = copper;
colormap(colmap)

if nargin > 1
  % Plot the path of a single bug.
  hold on
  %plot(bugpath(:,1)-Ngrid,bugpath(:,2)-Ngrid,'y.','MarkerSize',10)
  plot(bugpath(:,1)-Ngrid,bugpath(:,2)-Ngrid,'y-','LineWidth',1)
  hold off
end
