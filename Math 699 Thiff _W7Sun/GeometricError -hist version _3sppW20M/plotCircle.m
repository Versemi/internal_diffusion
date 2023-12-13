function [varargout] = plotCircle(r_input) 

Npart = 200000; 

% Need to change grid boundaries if there's a drift.
% Find an empirical formula first.
Ngrid = ceil(1.2*sqrt(Npart));  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
grid = zeros(2*Ngrid+1);
x = -Ngrid:Ngrid;



r_input = 1000; 

n = 100000*2;
theta = linspace(0,2*pi,n);

x_drawn = r_input*cos(theta)+grid0-0.5;
y_drawn = r_input*sin(theta)+grid0-0.5;


plot(x_drawn, y_drawn)

axis equal