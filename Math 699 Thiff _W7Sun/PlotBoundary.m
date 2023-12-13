% a combination of all three algorithms 

% _sppW18T
% _sppW19M

function [varargout] = PlotBoundary(Npart) 

notFinish = 1;  

while (notFinish == 1) 
    
    x = input('What do we do? \n 0. IDLA with 4 directions; \n 1. Build entire boundary; \n 2. Build boundary incrementally; \n 3. Build boundary by discrete laplace; \n 4. Compare the boundaries obtained from different algorithms; \n 5  Plot radius versus particle number; \n 6. Do nothing. \n \n');
    
    if (x==0)
        Npart = input('\n How many particles in total? \n');
        IDLA_fourDt(Npart)
    end 
    
    if (x==1)
        Npart = input('\n How many particles in total? \n');
        intdlaBuildEntBdry_suW17T(Npart)
    end 

    if (x==2) 
        Npart = input('\n How many particles in total? \n');
        intdlaBuildBdryIncrem_suW17T(Npart)
    end 
    
    if (x==3)
        Npart = input('\n How many particles in total? \n');
        intdlaDiscreteLaplace_suW17W(Npart)
    end
    
    if (x==4)
        Npart = input('\n How many particles in total? \n');
        compareBdry(Npart) 
    end 
    
    if (x==5)
        Npart = input('\n How many particles in total? \n');
        PlotRadiusVersusNumber(Npart)
    end
        
    if (x==6)
        notFinish = 0;
    end 
end 


% algr 1: Build Entire Boundary 
function intdlaBuildEntBdry_suW17T(n_particle)

Npart = n_particle; 

drift = [0 0];
% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

rng('default')


v_dir = [1 0; 0 1; -1 0; 0 -1];
n_dir = 4;

for i = 1:Npart
  X = [0 0];
  while 1
    
    d = randi(n_dir);
    
    X(1) = X(1) + v_dir(d, 1) + drift(1); 
    X(2) = X(2) + v_dir(d, 2) + drift(2);
    
    % record the number of steps each particles takes before it reaches an
    % empty spot, and store in a vector.  
    
    
    if X(2)+grid0 > 2*Ngrid+1
        xPos = 2*Ngrid+1;
    else
        xPos = X(2)+grid0; 
    end
     
    if X(1)+grid0 > 2*Ngrid+1
        yPos = 2*Ngrid+1;
    else
        yPos = X(1)+grid0;
    end
    
    if ~grid(yPos,xPos)
      grid(yPos,xPos) = 1;
      break
    end
    
  end
 
  
  % find the first boundary point
  
  y_loc = grid0;
  x_loc = 0; 
  
  isEnd = 0;
  
  for kk = 0:(grid0-1)
      x_loc2 = grid0 + kk; % start from the center of the circle and examine the grids outwards
      
      if (grid(y_loc, x_loc2) == 0 && isEnd == 0) % find the first grid which is not occupied
          x_loc = x_loc2 - 1;
          isEnd = 1;
      end 
  end
  
  % loop over all unoccupied neighbors until we find the first boundary point.
  
  gridB = zeros(size(grid));    % the boundary we want to draw and to eliminate the repeatable counting
  finish = 0;
  bdP = [0 0; y_loc x_loc]; % a matrix to store all the boundary points
  
  while finish == 0
      finish = 1;
      
      for k=1:size(bdP,1)   % examine all the boundary points stored in bdP
          if ~(bdP(k,1) == 0 && bdP(k,2) == 0)  % skip all the [0 0] array
              finish = 0;
              neighb = [bdP(k,1)+1 bdP(k,2); bdP(k,1)-1 bdP(k,2); bdP(k,1) bdP(k,2)+1; bdP(k,1) bdP(k,2)-1; bdP(k,1)+1 bdP(k,2)+1; bdP(k,1)+1 bdP(k,2)-1; bdP(k,1)-1 bdP(k,2)+1; bdP(k,1)-1 bdP(k,2)-1];
              % all the neighbors of the boundary points listed in bdP
              
              for j=1:8
                  if (grid(neighb(j,1), neighb(j,2))==1 && gridB(neighb(j,1), neighb(j,2))==0)  
                      % if this neighbor is occupied and hasn't been stored in gridB
                      
                    if ~((grid(neighb(j,1)+1, neighb(j,2))==1) && (grid(neighb(j,1)-1, neighb(j,2))==1) && (grid(neighb(j,1), neighb(j,2)+1)==1) && (grid(neighb(j,1), neighb(j,2)-1)==1))
                        % and if this neighbor is not surrounded by 4 other neighbors 
                        
                        bdP = [bdP; neighb(j,1) neighb(j,2)];   % then this neighbor is a boundary point and add it to bdP
                        gridB(neighb(j,1), neighb(j,2))=1;  % mark it in gridB as 1
                    end
                  end
              end 
              
              gridB(bdP(k,1),bdP(k,2)) = 1; % the first boundary is already used and mark it as 1 in gridB
              bdP(k,1) = 0; 
              bdP(k,2) = 0; % set this point as [0 0] in bdP so it would be skipped
          end 
      end 
  end 
  
    if ~mod(i,100)
        fprintf('particle %g\n',i)
        plotgridB(x,gridB)
        drawnow
    end

end

plotgridB(x,gridB)


% algr 2: build the boundary incrementally 
function intdlaBuildBdryIncrem_suW17T(n_particle) 
%INTDLA   Simulate internal diffusion-limited aggregation.
%   INTDLA(N) simulates an internal DLA process with N particles.
%   (If unspecified N defaults to 10000.)
%
%   G = INTDLA(N) also returns the final occupation grid.

Npart = n_particle; 

drift = [0 0];
% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

rng('default')

v_dir = [1 0; 0 1; -1 0; 0 -1];
n_dir = 4;

for i = 1:Npart
  X = [0 0];
  while 1
    
    d = randi(n_dir);
    
    X(1) = X(1) + v_dir(d, 1) + drift(1); 
    X(2) = X(2) + v_dir(d, 2) + drift(2);
    
    % record the number of steps each particles takes before it reaches an empty grid, and store in a vector.  
    
    if X(2)+grid0 > 2*Ngrid+1
        xPos = 2*Ngrid+1;
    else
        xPos = X(2)+grid0; 
    end
     
    if X(1)+grid0 > 2*Ngrid+1
        yPos = 2*Ngrid+1;
    else
        yPos = X(1)+grid0;
    end
    
    if ~grid(yPos,xPos)
      grid(yPos,xPos) = 1;
      break
    end
    
  end  

    gridB = grid;
 
    for k = 2:(2*grid0-1)
        for l = 2:(2*grid0-1)
            if grid(k, l) == 1 && grid(k, l+1) == 1 && grid (k, l-1) == 1 && grid (k+1, l) == 1 && grid (k-1, l) == 1
                gridB(k, l) = 0; 

            end 
        end 
    end
    
    if ~mod(i,100)
        fprintf('particle %g\n',i)
        plotgridB(x,gridB)
        drawnow
    end

end 

plotgridB(x,gridB)


% algr 3: Discrete Laplace
function intdlaDiscreteLaplace_suW17W(n_particle)
%INTDLA   Simulate internal diffusion-limited aggregation.
%   INTDLA(N) simulates an internal DLA process with N particles.
%   (If unspecified N defaults to 10000.)
%
%   G = INTDLA(N) also returns the final occupation grid.

Npart = n_particle; 

drift = [0 0];
% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant
grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

rng('default')


v_dir = [1 0; 0 1; -1 0; 0 -1];
n_dir = 4;

for i = 1:Npart
  X = [0 0];
  while 1
    
    d = randi(n_dir);
    
    X(1) = X(1) + v_dir(d, 1) + drift(1);
    X(2) = X(2) + v_dir(d, 2) + drift(2);
    
    % record the number of steps each particles takes before it reaches an
    % empty sopt, and store in a vector.  
    
    
    if X(2)+grid0 > 2*Ngrid+1
        xPos = 2*Ngrid+1;
    else 
        xPos = X(2)+grid0; 
        if xPos <= 0 
            xPos = 1;
        end      
    end
    
    
    if X(1)+grid0 > 2*Ngrid+1
        yPos = 2*Ngrid+1;
    else
        yPos = X(1)+grid0;
    end
    
    if ~grid(yPos,xPos)
      grid(yPos,xPos) = 1;
      break
    end
    
  end
  
  % algr 3: Discrete Laplace 
    gridB = grid;

    for k = 2:(2*Ngrid)
        for l = 2:(2*Ngrid)
            delta_f = 4*grid(k, l) - grid(k, l+1) - grid(k, l-1) - grid(k+1, l) - grid(k-1, l);
            if delta_f == 0  
                gridB(k, l) = 0; 

            end 
        end 
    end
    
    if ~mod(i,100)
        fprintf('particle %g\n',i)
        plotgridB(x,gridB)
        drawnow
    end

end 

plotgridB(x,gridB)

%==========================================================================
function plotgridB(x,gridB)

imagesc(x,x,gridB.'), axis square, axis xy
xlabel('x'), ylabel('y')


% Compare the boundaries obtained from different algorithms

% problems: the cavity shown might cause difference between the boundaries
% obtained from the above 3 algorithms. 

function compareBdry(n_particle) 
Npart = n_particle; 

drift = [0 0];
% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

rng('default')


v_dir = [1 0; 0 1; -1 0; 0 -1];
n_dir = 4;


for i = 1:Npart
      X = [0 0];
      while 1

        d = randi(n_dir);

        X(1) = X(1) + v_dir(d, 1) + drift(1); 
        X(2) = X(2) + v_dir(d, 2) + drift(2);

        % record the number of steps each particles takes before it reaches an
        % empty sopt, and store in a vector.     
        if X(2)+grid0 > 2*Ngrid+1
            xPos = 2*Ngrid+1;
        else
            xPos = X(2)+grid0; 
        end

        if X(1)+grid0 > 2*Ngrid+1
            yPos = 2*Ngrid+1;
        else
            yPos = X(1)+grid0;
        end

        if ~grid(yPos,xPos)
          grid(yPos,xPos) = 1;
          break
        end

      end

      % algr 1
      
      % find the first boundary point
      y_loc = grid0;
      x_loc = 0; 

      isEnd = 0;

      for kk = 0:(grid0-1)
          x_loc2 = grid0 + kk; % start from the center of the circle and examine the grids outwards

          if (grid(y_loc, x_loc2) == 0 && isEnd == 0) % find the first grid which is not occupied
              x_loc = x_loc2 - 1;
              isEnd = 1;
          end 
      end

      % loop over all unoccupied neighbors until we find the first boundary point.

      gridB1 = zeros(size(grid));    % the boundary we want to draw and to eliminate the repeatable counting
      finish = 0;
      bdP = [0 0; y_loc x_loc]; % a matrix to store all the boundary points

      while finish == 0
          finish = 1;

          for k=1:size(bdP,1)   % examine all the boundary points stored in bdP
              if ~(bdP(k,1) == 0 && bdP(k,2) == 0)  % skip all the [0 0] array
                  finish = 0;
                  neighb = [bdP(k,1)+1 bdP(k,2); bdP(k,1)-1 bdP(k,2); bdP(k,1) bdP(k,2)+1; bdP(k,1) bdP(k,2)-1; bdP(k,1)+1 bdP(k,2)+1; bdP(k,1)+1 bdP(k,2)-1; bdP(k,1)-1 bdP(k,2)+1; bdP(k,1)-1 bdP(k,2)-1];
                  % all the neighbors of the boundary points listed in bdP

                  for j=1:8
                      if (grid(neighb(j,1), neighb(j,2))==1 && gridB1(neighb(j,1), neighb(j,2))==0)  
                          % if this neighbor is occupied and hasn't been stored in gridB

                        if ~((grid(neighb(j,1)+1, neighb(j,2))==1) && (grid(neighb(j,1)-1, neighb(j,2))==1) && (grid(neighb(j,1), neighb(j,2)+1)==1) && (grid(neighb(j,1), neighb(j,2)-1)==1))
                            % and if this neighbor is not surrounded by 4 other neighbors 

                            bdP = [bdP; neighb(j,1) neighb(j,2)];   % then this neighbor is a boundary point and add it to bdP
                            gridB1(neighb(j,1), neighb(j,2))=1;  % mark it in gridB as 1
                        end
                      end
                  end 

                  gridB1(bdP(k,1),bdP(k,2)) = 1; % the first boundary is already used and mark it as 1 in gridB
                  bdP(k,1) = 0; 
                  bdP(k,2) = 0; % set this point as [0 0] in bdP so it would be skipped
              end 
          end 
      end
      
     

      % algr 2
      gridB2 = grid;
        for k = 2:(2*Ngrid-1)
            for l = 2:(2*Ngrid-1)
                if grid(k, l) == 1 && grid(k, l+1) == 1 && grid (k, l-1) == 1 && grid (k+1, l) == 1 && grid (k-1, l) == 1
                    gridB2(k, l) = 0; 

                end 
            end 
        end

      % algr 3
      gridB3 = grid;
        for k = 2:(2*Ngrid-1)
            for l = 2:(2*Ngrid-1)
                delta_f = 4*grid(k, l) - grid(k, l+1) - grid(k, l-1) - grid(k+1, l) - grid(k-1, l);
                if delta_f == 0  
                    gridB3(k, l) = 0; 
                end 
            end 
        end
        
      
    % compare 3 boundaries 
   
    theyAgree = 0; 
    for m = 1:(2*Ngrid)
        for n = 1:(2*Ngrid)

          if ((gridB1(m,n)==gridB2(m,n)) && (gridB2(m,n)==gridB3(m,n)) && (gridB1(m,n)==gridB3(m,n)))
              theyAgree = 1; 

          else 
              theyAgree = 0;
              disp('They do not agree with each other. The step number i = ')
              disp(i)
          end  

        end
    end
    

end

   
if (theyAgree == 1)
        fprintf('\n They agree with each other in every step. \n \n ')
end 
   


% 5. plot radius versus particle number 
function PlotRadiusVersusNumber(n_particle)

Npart = n_particle; 

drift = [0 0];
% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

rng('default')

v_dir = [1 0; 0 1; -1 0; 0 -1];
n_dir = 4;



allAverage = [];
allStd = [];

for i = 1:Npart

  radius_collection = [];
  
  X = [0 0];
  while 1
    
    d = randi(n_dir);
    
    X(1) = X(1) + v_dir(d, 1) + drift(1); 
    X(2) = X(2) + v_dir(d, 2) + drift(2);
    
    % record the number of steps each particles takes before it reaches an empty grid, and store in a vector.  
    
    if X(2)+grid0 > 2*Ngrid+1
        xPos = 2*Ngrid+1;
    else
        xPos = X(2)+grid0; 
    end
     
    if X(1)+grid0 > 2*Ngrid+1
        yPos = 2*Ngrid+1;
    else
        yPos = X(1)+grid0;
    end
    
    if ~grid(yPos,xPos)
      grid(yPos,xPos) = 1;
      break
    end
    
  end  

  
  % algr 2
    gridB = grid;
 
    for k = 2:(2*Ngrid)
        for l = 2:(2*Ngrid)
            if grid(k, l) == 1 && grid(k, l+1) == 1 && grid (k, l-1) == 1 && grid (k+1, l) == 1 && grid (k-1, l) == 1
                gridB(k, l) = 0; 
            end 
        end 
    end
   
    % calculate the radius in this step: r = ((x^2)+y^2)^(1/2)
    for m = 1:(2*Ngrid+1)
        for n = 1:(2*Ngrid+1)
            if (gridB(m,n) == 1)
                radius = sqrt((m - grid0).^2 + (n - grid0).^2); 

                radius_collection = [radius_collection, radius];
                
            end
        end 
    end 
    
    radius_average = mean(radius_collection);
    standardDev = std(radius_collection);
 
    allAverage = [allAverage, radius_average];
    allStd = [allStd, standardDev];
    

end 

disp("std= "+allStd(2*Ngrid+1))


f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = ''; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';

subplot(1,2,1,'Parent',p) 
loglog(allAverage)
title('Average Distance versus Particle Numbers') 
xlabel('particle number') % x-axis label
ylabel('average distance') % y-axis label

subplot(1,2,2,'Parent',p) 
loglog(allStd)
title('Standard Deviation versus Particle Numbers')
xlabel('particle number') % x-axis label
ylabel('standard deviation') % y-axis label



%IDLA in four directions
function IDLA_fourDt(n_particle) % 8/28 _3suW11T
%INTDLA   Simulate internal diffusion-limited aggregation.
%   INTDLA(N) simulates an internal DLA process with N particles.
%   (If unspecified N defaults to 10000.)
%
%   G = INTDLA(N) also returns the final occupation grid.

Npart = n_particle; 

drift = [0 0];
% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

rng('default')


v_dir = [1 0; 0 1; -1 0; 0 -1];
n_dir = 4;

for i = 1:Npart
  X = [0 0];
  while 1

    d = randi(n_dir);

    X(1) = X(1) + v_dir(d, 1) + drift(1); 
    X(2) = X(2) + v_dir(d, 2) + drift(2);

    % record the number of steps each particles takes before it reaches an
    % empty sopt, and store in a vector.  


    if X(2)+grid0 > 2*Ngrid+1
        xPos = 2*Ngrid+1;
    else
        xPos = X(2)+grid0; 
    end

    if X(1)+grid0 > 2*Ngrid+1
        yPos = 2*Ngrid+1;
    else
        yPos = X(1)+grid0;
    end

    if ~grid(yPos,xPos)
      grid(yPos,xPos) = 1;
      break
    end

   if ~mod(i,100)
    fprintf('particle %g\n',i)
    plotgrid(x,grid)
    drawnow
   end
  end
end

plotgrid(x,grid)
    

function plotgrid(x,grid)

imagesc(x,x,grid.'), axis square, axis xy
xlabel('x'), ylabel('y')

