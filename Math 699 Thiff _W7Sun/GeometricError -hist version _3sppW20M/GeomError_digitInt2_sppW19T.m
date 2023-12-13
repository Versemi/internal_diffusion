% calculate the geometric error _sppW19T 6/5
function [Rstd,Ravg,Nint] = GeometricError0(R)

if nargin < 1, R = 100; end

if length(R) > 1
  Rstd = []; Ravg = []; Nint = [];
  for i = 1:length(R)
    [Rstd1,Ravg1,Nint1] = GeometricError(R(i));
    Rstd = [Rstd Rstd1];
    Ravg = [Ravg Ravg1];
    Nint = [Nint Nint1];
  end
  return
end

% Input the radius (The average radius is 55.9415 for 10000 particles in total)

% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(R)+1;  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

dtheta = .01/R;
n = ceil(2*pi/dtheta);
theta = linspace(0,2*pi,n);

x_drawn = R*cos(theta)+grid0-1;
y_drawn = R*sin(theta)+grid0-1;

geomBdry = zeros(2*Ngrid+1);

geomRadiusCollect = [];
number_geomBdry = 0; 

for ii= 1:n
    x_temp = int32(x_drawn(ii))+1; 
    y_temp = int32(y_drawn(ii))+1; 
    
    if (geomBdry(y_temp, x_temp)==0)
        geomBdry(y_temp, x_temp)=1;

        Radius_geomBdry = sqrt((double(y_temp) - grid0).^2 + (double(x_temp) - grid0).^2); 
        geomRadiusCollect = [geomRadiusCollect, Radius_geomBdry];
        
        number_geomBdry = number_geomBdry + 1; 
    end 
end 

number_insideCircle = 0; 

for jj=1:(2*Ngrid+1)
    for kk=1:(2*Ngrid+1)
        if (jj <= y_temp) && (kk <= x_temp)
            number_insideCircle = number_insideCircle +1;
        end 
    end 
end 

plotGeomBdry(x,geomBdry)


Ravg = mean(geomRadiusCollect);
Rstd = std(geomRadiusCollect);
Nint = number_insideCircle;

if 0
  disp("The average radius of the grids on the geometric boundary is " + geomRadius_average)
  disp("The standard deviation of the grids on the geometric boundary is " + geomRadius_standardDev)
  disp("The number of the grids inside the cirlce is " + number_insideCircle)
end

function plotGeomBdry(x,GeomBdry)

imagesc(x,x,GeomBdry.'), axis square, axis xy
xlabel('x'), ylabel('y')
