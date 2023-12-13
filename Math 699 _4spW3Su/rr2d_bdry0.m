function [varargout] = rr2d_bdry(Nbugs,gridstruct0) %flag
%RR2D   Rotor-router walk in two dimensions.
%   GS = rr2d_bdry(Nbugs) gives 2 plots of the distance between center of
%   each occupied grid and the origin X0 versus Nbugs/Npart, and the standard
%   deviation of that versus Nbugs/Npart.
%
%   GS = rr2d_bdry(Nbugs,GS) continues the walk for more Nbugs/Npart and
%   returns 2 plots.

if nargin < 1, Nbugs = 1000; end

% Size of grid quadrant for a given number of bugs.
gridquadsize = @(Nbugs) ceil(.6*sqrt(Nbugs)); %?

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


all_dist=[]; avg_dist=[]; all_avg=[]; all_std=[]; % array to store the average distance and std

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
    
    gridstruct.grid = grid; 
    gridstruct.Nbugs = Nbugs0+Nbugs;
    
% select out the boundary 
        gridB = gridstruct.grid;

        for k = 2:(2*Ngrid)
            for l = 2:(2*Ngrid)            
                if ((gridstruct.grid(k, l)~=0) && (gridstruct.grid(k, l+1)~=0) && (gridstruct.grid (k, l-1)~=0) && (gridstruct.grid (k+1, l)~=0) && (gridstruct.grid (k-1, l)~=0))
                    gridB(k, l) = 0; 
                end 
            end 
        end

        for kk = 2:(2*Ngrid)
            for ll = 2:(2*Ngrid)            
                if gridB(kk,ll)~=0
                    gridB(kk, ll) = 1; 
                end 
            end 
        end

        
        % calculate the distance between center of each grid and X0        
        for m = 1:(2*Ngrid+1)
            for n = 1:(2*Ngrid+1)
                if (gridB(m,n) == 1)
                    dist = sqrt((m - X0).^2 + (n - X0).^2); 
                    all_dist = [all_dist, dist];               
                end
            end 
        end
        
        avg_dist = mean(all_dist);
        std_dist = std(all_dist);

        all_avg = [all_avg avg_dist];
        all_std = [all_std std_dist];  
  end
end

sz_all_avg = size(all_avg);

disp("average distant of roter-router = "+all_avg(sz_all_avg(2)))
disp("std of roter-router = "+all_std(sz_all_avg(2)))

all_std_dist = all_std./avg_dist;

%disp("all_std = "+all_std)


%%%%%%
% IDLA boundary

Npart1 = Nbugs; 

drift1 = [0 0];
% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid1 = ceil(1.2*sqrt(Npart1));  % size of grid quadrant

grid01 = Ngrid1+1;                % center of grid
grid1 = zeros(2*Ngrid1+1);
x1 = -Ngrid1:Ngrid1;

rng('default')

v_dir = [1 0; 0 1; -1 0; 0 -1];
n_dir = 4;


all_avg1 = [];
all_std1 = [];

for i = 1:Npart1

  all_dist1 = [];
  
  X1 = [0 0];
  while 1
    
    d = randi(n_dir);
    
    X1(1) = X1(1) + v_dir(d, 1) + drift1(1); 
    X1(2) = X1(2) + v_dir(d, 2) + drift1(2);
    
    % record the number of steps each particles takes before it reaches an empty grid, and store in a vector.  
    
    if X1(2)+grid01 > 2*Ngrid1+1
        xPos = 2*Ngrid1+1;
    else
        xPos = X1(2)+grid01; 
    end
     
    if X1(1)+grid01 > 2*Ngrid1+1
        yPos = 2*Ngrid1+1;
    else
        yPos = X1(1)+grid01;
    end
    
    if ~grid1(yPos,xPos)
      grid1(yPos,xPos) = 1;
      break
    end
    
  end  

  
  % algr 2
    gridB1 = grid1;
 
    for k = 2:(2*Ngrid1)
        for l = 2:(2*Ngrid1)
            if grid1(k, l) == 1 && grid1(k, l+1) == 1 && grid1 (k, l-1) == 1 && grid1 (k+1, l) == 1 && grid1 (k-1, l) == 1
                gridB1(k, l) = 0; 
            end 
        end 
    end
   
    % calculate the distances in this step: r = ((x^2)+y^2)^(1/2)
    for m = 1:(2*Ngrid1+1)
        for n = 1:(2*Ngrid1+1)
            if (gridB1(m,n) == 1)
                dist1 = sqrt((m - grid01).^2 + (n - grid01).^2); 

                all_dist1 = [all_dist1, dist1];
                
            end
        end 
    end 
    
    avg_dist1 = mean(all_dist1);
    std1 = std(all_dist1);
 
    all_avg1 = [all_avg1, avg_dist1];
    all_std1 = [all_std1, std1];
    

end 

%sz_all_avg1 = size(all_avg1)
disp("average distant of IDLA = "+all_avg1(Npart1))
disp("std of IDLA = "+all_std1(Npart1))


%%%%%%
%plot 
f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = ''; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';

subplot(1,2,1,'Parent',p) 
title('Average Distance versus Npart & Nbugs') 
xlabel('Npart & Nbugs') % x-axis label
ylabel('D') % y-axis label

hold on
loglog(all_avg1)
loglog(all_avg)
legend({'IDLA','Roter-Router'},'Location','southwest')
hold off


subplot(1,2,2,'Parent',p)
title('Standard Deviation versus Npart & Nbugs')
xlabel('Npart & Nbugs') % x-axis label
ylabel('Standard Deviation') % y-axis label

hold on
loglog(all_std1)
loglog(all_std)
loglog(all_std_dist)
legend({'IDLA','Roter-Router','Roter-Router Scaled by D'},'Location','southwest')
hold off

if nargout > 0
  varargout{1} = gridstruct;
end
