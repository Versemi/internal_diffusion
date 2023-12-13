% calculate the geometric error _sppW20W 6/13

%{ 
Task: Draw a cirlce on the grids and mark all the grids which are passed by 
the boundary of the circle. Calculate the distance between the center of
each colored grids and the origin. Evaluate the average distance and the
standard deviation. Compute the number of the grids inside the circle
drawn.
%}

function [Rstd,Ravg,Nint] = GeometricError02(R)

if nargin < 1, R = 10; end

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


m1 = ceil(R - 1/2); n1 = 0; 

% find the boundary points on the upperRightQuarter
m_upperRightQuarter = []; % horizontal value of the coordinate
n_upperRightQuarter = []; % vertical value of the coordinate  

this_m = m1; this_n = n1; 
notEnd = 1; 

while (notEnd == 1)
    m_upper = this_m; n_upper = this_n + 1;
    m_left = this_m - 1; n_left = this_n; 
 
    m_neighbor = [m_upper m_left];
    n_neighbor = [n_upper n_left]; 

    for j=1:size(m_neighbor)
        ratio = n_neighbor(j)./(m_neighbor(j)+0.5);
        th = atan(ratio);
        x_intersect = R*cos(th);
        y_intersect = R*sin(th);
        
        if (x_intersect <= m_neighbor(j)+0.5) && ((m_neighbor(j)-0.5) <= x_intersect) && (y_intersect <= n_neighbor(j)+0.5) && ((n_neighbor(j)-0.5) <= y_intersect)
            this_m = m_neighbor(j);  
            this_n = n_neighbor(j);sw
            
            m_upperRightQuarter = [m_upperRightQuarter, this_m]; 
            n_upperRightQuarter = [n_upperRightQuarter, this_n];
            
        end

        if (this_m == 1) && (this_n == ceil(R-1/2))
            notEnd = 0; 
        end 

    end 
       
end 


all_radius = [];

for k=1:size(m_upperRightQuarter)
    radius = sqrt(m_upperRightQuarter(k).^2 + n_upperRightQuarter(k).^2); 
    all_radius = [all_radius, radius];
end 

radius_m1 = R; 
all_radius = [all_radius, radius_m1];
all_radius = [all_radius, radius_m1];
all_radius = [all_radius, radius_m1];
all_radius = [all_radius, radius_m1];

num_inside = m_upperRightQuarter(1)-1;

for kk=2:size(m_upperRightQuarter)
    if ~(n_upperRightQuarter(kk)==n_upperRightQuarter(kk-1))
        num_inside = num_inside + m(kk) -1; 
    end 
end 

% count the rest grids in the center
num_inside = num_inside + 4*(m1-1) + 1; 

Ravg = mean(all_radius);
Rstd = std(all_radius);
Nint = num_inside;
