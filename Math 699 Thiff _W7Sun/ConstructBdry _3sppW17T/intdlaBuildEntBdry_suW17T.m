% algr 1: Build Entire Boundary 

function [varargout] = intdlaBuildEntBdry_suW17T(Npart)
%INTDLA   Simulate internal diffusion-limited aggregation.
%   INTDLA(N) simulates an internal DLA process with N particles.
%   (If unspecified N defaults to 10000.)
%
%   G = INTDLA(N) also returns the final occupation grid.

if nargin < 1
  Npart = 10000;                % number of particles
end
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
  
  
  % algr 1: Build Entire Boundary  

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
        plotgrid(x,gridB)
        drawnow
    end

end


if nargout > 0
  varargout{1} = grid;
end  


    
            


%==========================================================================
function plotgrid(x,gridB)

imagesc(x,x,gridB.'), axis square, axis xy
xlabel('x'), ylabel('y')
