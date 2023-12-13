function grid = rr2d_arrows(gridstruct)
%RR2D_ARROWS   Convert RR2D grid to arrow directions.
%   G = RR2D_ARROWS(GS) converts the struct GS containing the number of
%   times each site was visited to arrows.  The arrows correspond to
%
%     0 = site never visited
%     1 = arrow E
%     2 = arrow N
%     3 = arrow W
%     4 = arrow S
%
%   See also RR2D, RR2D_PLOT.

if isstruct(gridstruct)
  grid = gridstruct.grid;
else
  grid = gridstruct;
end

% Map 0 to 0, and occupied sites to directions 1,2,3,4.
grid(grid == 0) = NaN;   % Hide the 0 sites from the mod operation.
grid = mod(grid-1,4)+1;
grid(isnan(grid)) = 0;   % Restore 0 sites.
