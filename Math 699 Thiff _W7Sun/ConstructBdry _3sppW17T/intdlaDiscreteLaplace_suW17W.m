% algr 3: Discrete Laplace


function [varargout] = intdlaDiscreteLaplace_suW17W(Npart)
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
    
    
    if X(1)+grid0 > 2*Ngrid+1
        xPos = 2*Ngrid+1;
    else 
        xPos = X(1)+grid0; 
        if xPos <= 0 
            xPos = 1;
        end      
    end
    
    
    if X(2)+grid0 > 2*Ngrid+1
        yPos = 2*Ngrid+1;
    else
        yPos = X(2)+grid0;
    end
    
    if ~grid(yPos,xPos)
      grid(yPos,xPos) = 1;
      break
    end
    
  end
  
  % algr 3: Discrete Laplace 
    grid_bdy = grid;

    for k = 2:(2*Ngrid)
        for l = 2:(2*Ngrid)
            delta_f = 4*grid(k, l) - grid(k, l+1) - grid(k, l-1) - grid(k+1, l) - grid(k-1, l);
            if delta_f == 0  
                grid_bdy(k, l) = 0; 

            end 
        end 
    end
    
    if ~mod(i,100)
        fprintf('particle %g\n',i)
        plotgrid(x,grid_bdy)
        drawnow
    end

end 

plotgrid(x,grid_bdy)

if nargout > 0
  varargout{1} = grid;
end
            


%==========================================================================
function plotgrid(x,grid_bdy)

imagesc(x,x,grid_bdy.'), axis square, axis xy
xlabel('x'), ylabel('y')
