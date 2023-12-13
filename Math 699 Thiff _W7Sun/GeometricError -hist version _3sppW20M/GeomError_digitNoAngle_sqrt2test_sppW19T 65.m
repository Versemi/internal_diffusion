% calculate the geometric error _sppW19T 6/5
function [varargout] = GeomError_digitPt2_sqrt2test_sppW19T(radius_average) 

radius_average = input('\n Input the radius (The average radius is 55.9415 for 10000 particles in total): \n');
inputRadius(radius_average)


function inputRadius(radius_input)

Npart = 20000; 

% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

radius_average = radius_input; 

%n = 1000;
%theta = linspace(0,2*pi,n);
%x_drawn = radius_average*cos(theta)+grid0;
%y_drawn = radius_average*sin(theta)+grid0;


geomBdry = zeros(2*Ngrid+1); 
for ii= 1:(2*Ngrid+1)
    for jj = 1:(2*Ngrid+1)
        %for kk=1:n
            Radius_geomBdry = sqrt((ii - grid0).^2 + (jj - grid0).^2);
            %if (abs(ii - x_drawn(kk)) <= 1/2) && (abs(jj - y_drawn(kk)) <= 1/2) && (0 <= abs(ii - x_drawn(kk))) && (0 <= abs(jj - y_drawn(kk)))
            if (abs(Radius_geomBdry - radius_average) <= sqrt(2)/2) && (0 <= abs(Radius_geomBdry - radius_average))
                geomBdry(ii, jj) = 1; 
            end 
        %end 
    end
end 

plotGeomBdry(x,geomBdry)



geomRadiusCollect = [];

% calculate the radius: r = ((x^2)+y^2)^(1/2)
for mm = 1:(2*Ngrid+1)
    for nn = 1:(2*Ngrid+1)
        if (geomBdry(mm,nn) == 1)
            Radius_geomBdry = sqrt((mm - grid0).^2 + (nn - grid0).^2); 

            geomRadiusCollect = [geomRadiusCollect, Radius_geomBdry];
        end
    end 
end 

geomRadius_average = mean(geomRadiusCollect);
geomRadius_standardDev = std(geomRadiusCollect);

disp("The average radius of the grids on the geometric boundary is " + geomRadius_average)
disp("The standard deviation of the grids on the geometric boundary is " + geomRadius_standardDev)


number_geomBdry = 0; 
for ii= 1:(2*Ngrid+1)
    for jj = 1:(2*Ngrid+1)
        if (geomBdry(ii,jj) == 1)
            number_geomBdry = number_geomBdry + 1; 
        end 
    end
end 

disp("The number of the grids on the geometric boundary is " + number_geomBdry)


function inputParticleNum(n_particle)

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

end 

% algr 2 to find the boundary in last step
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

disp("The average radius is " + radius_average + " for " + Npart + " particles in total")

geomBdry = zeros(2*Ngrid+1); 
for ii= 1:(2*Ngrid+1)
    for jj = 1:(2*Ngrid+1)
        Radius_geomBdry = sqrt((ii - grid0).^2 + (jj - grid0).^2);
        if (abs(Radius_geomBdry - radius_average) <= sqrt(2)/2) && (0 <= abs(Radius_geomBdry - radius_average))
            geomBdry(ii, jj) = 1; 
        end 
    end
end 

plotGeomBdry(x,geomBdry)


geomRadiusCollect = [];

% calculate the radius: r = ((x^2)+y^2)^(1/2)
for mm = 1:(2*Ngrid+1)
    for nn = 1:(2*Ngrid+1)
        if (geomBdry(mm,nn) == 1)
            Radius_geomBdry = sqrt((mm - grid0).^2 + (nn - grid0).^2); 

            geomRadiusCollect = [geomRadiusCollect, Radius_geomBdry];
        end
    end 
end 

geomRadius_average = mean(geomRadiusCollect);
geomRadius_standardDev = std(geomRadiusCollect);

disp("The average radius of the grids on the geometric boundary is " + geomRadius_average)
disp("The standard deviation of the grids on the geometric boundary is " + geomRadius_standardDev)


number_geomBdry = 0; 
for ii= 1:(2*Ngrid+1)
    for jj = 1:(2*Ngrid+1)
        if (geomBdry(ii,jj) == 1)
            number_geomBdry = number_geomBdry + 1; 
        end 
    end
end 

disp("The number of the grids on the geometric boundary is " + number_geomBdry)



function plotGeomBdry(x,GeomBdry)

imagesc(x,x,GeomBdry.'), axis square, axis xy
xlabel('x'), ylabel('y')

