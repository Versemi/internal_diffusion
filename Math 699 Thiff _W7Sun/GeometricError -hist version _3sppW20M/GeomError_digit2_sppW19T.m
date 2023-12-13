% calculate the geometric error _sppW19T 6/5
function [varargout] = GeomError_digitPt_sppW19T(radius_average) 

radius_average = input('\n Input the radius (The average radius is 55.9415 for 10000 particles in total): \n');
inputRadius(radius_average)


function inputRadius(radius_input)

Npart = 2000000; 


% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

radius_average = radius_input; 

n = 1000;
theta = linspace(0,2*pi,n);
x_drawn = radius_average*cos(theta)+grid0;
y_drawn = radius_average*sin(theta)+grid0;


geomBdry = zeros(2*Ngrid+1); 
for ii= 1:(2*Ngrid+1)
    for jj = 1:(2*Ngrid+1)
        for kk = 1:n
            % Radius_geomBdry = sqrt((ii - grid0).^2 + (jj - grid0).^2);
            if (abs(ii - x_drawn(kk)) <= 1/2) && (abs(jj - y_drawn(kk)) <= 1/2) && (0 <= abs(ii - x_drawn(kk))) && (0 <= abs(jj - y_drawn(kk)))
                geomBdry(ii, jj) = 1; 
            end 
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



