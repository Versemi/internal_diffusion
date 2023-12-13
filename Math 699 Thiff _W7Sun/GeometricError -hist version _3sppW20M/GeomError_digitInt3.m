% calculate the geometric error _sppW19T 6/5
function [varargout] = GeometricError2(r_input) 

notFinish = 1;  

while (notFinish == 1) 
    x = input('What do we do? \n 1. Input the radius and plot a discretized circle; \n 2. See how the average radius of the discretized circle and the standard deviation vary with respect to the given radius; \n 3. Do nothing. \n \n')

    if (x==1)
        r_input = input('\n Input the radius (The average radius is 55.9415 for 10000 particles in total): \n');
        inputRadius(r_input)
    end 

    if (x==2)
        r_input = 5000;
        AverageStdVersusR(r_input)
    end 
    
    if (x==3)
        notFinish = 0;  
    end 

end 




function inputRadius(radius_input)

Npart = 20000000; 

% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;

r_input = radius_input; 

n = 1000000;
theta = linspace(0,2*pi,n);

x_drawn = r_input*cos(theta)+grid0-0.5;
y_drawn = r_input*sin(theta)+grid0-0.5;


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


geomRadius_average = mean(geomRadiusCollect);
geomRadius_standardDev = std(geomRadiusCollect);

disp("The average radius of the grids on the geometric boundary is " + geomRadius_average)
disp("The standard deviation of the grids on the geometric boundary is " + geomRadius_standardDev)
disp("The number of the grids inside the cirlce is " + number_insideCircle)


function plotGeomBdry(x,GeomBdry)

imagesc(x,x,GeomBdry.'), axis square, axis xy
xlabel('x'), ylabel('y')



function AverageStdVersusR(radius_max)

r_input = radius_max; 

Npart = 20000000; 

% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;


allAverage = [];
allStd = [];

for m=1:r_input 

    n = 1000000;
    theta = linspace(0,2*pi,n);

    x_drawn = m*cos(theta)+grid0-0.5;
    y_drawn = m*sin(theta)+grid0-0.5;


    geomBdry = zeros(2*Ngrid+1);

    geomRadiusCollect = [];

    for ii= 1:n
        x_temp = int32(x_drawn(ii))+1; 
        y_temp = int32(y_drawn(ii))+1; 

        if (geomBdry(y_temp, x_temp)==0)
            geomBdry(y_temp, x_temp)=1;

            Radius_geomBdry = sqrt((double(y_temp) - grid0).^2 + (double(x_temp) - grid0).^2); 
            geomRadiusCollect = [geomRadiusCollect, Radius_geomBdry];
        end 
    end 

    geomRadius_average = mean(geomRadiusCollect);
    geomRadius_standardDev = std(geomRadiusCollect);
    
    allAverage = [allAverage, geomRadius_average];
    allStd = [allStd, geomRadius_standardDev];

end 

f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = ''; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';

subplot(1,2,1,'Parent',p) 
loglog(allAverage)
title('Average Radius versus Given Radius') 
xlabel('Given Radius') % x-axis label
ylabel('average radius R_n') % y-axis label

subplot(1,2,2,'Parent',p) 
loglog(allStd)
title('Standard Deviation versus Given Radius')
xlabel('Given Radius') % x-axis label
ylabel('standard deviation') % y-axis label



