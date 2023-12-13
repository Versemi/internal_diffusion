% calculate the geometric error _sppW19T 6/5
function [Rstd,Ravg,Nint] = GeometricError01(R)

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

Ngrid = ceil(R)+1;  % size of grid quadrant

grid0 = Ngrid+1;                % center of grid
x = -Ngrid:Ngrid;


m1 = ceil(R - 1/2); n1 = 1; 
m2 = -m1; n2 = -n1; 
m3 = 0; n3 = ceil(R - 1/2);
m4 = -m3; n4 = -n3; 


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

    % draw a circle
    dtheta = .01/R;
    num = ceil(2*pi/dtheta);
    theta = linspace(0,2*pi,num);

    x_drawn = R*cos(theta)+grid0-1;
    y_drawn = R*sin(theta)+grid0-1;

    for j=1:size(m_neighbor)

        %{
        ratio = n_neighbor(j)./m_neighbor(j);
        th = atan(ratio);
        m_intersect = R*cos(th);
        n_intersect = R*sin(th);
        %}
        
        %disp(n_intersect)
        
        for j2=1:num
            if (x_drawn(j2) <= m_neighbor(j)-1/2) && ((m_neighbor(j)-1/2-1) <= x_drawn(j2)) && (y_drawn(j2) <= n_neighbor(j)-1/2) && ((n_neighbor(j)-1/2-1) <= y_drawn(j2))
                this_m = m_neighbor(j); 
                this_n = n_neighbor(j); 
            end
        end 

        m_upperRightQuarter = [m_upperRightQuarter, this_m]; 
        n_upperRightQuarter = [n_upperRightQuarter, this_n];


        disp(n_upperRightQuarter)

        %{
        if (m_intersect <= m_neighbor(j)) && ((m_neighbor(j)-1) <= m_intersect) && (n_intersect <= n_neighbor(j)) && ((n_neighbor(j)-1) <= n_intersect)
            this_m = m_neighbor(j); 
            this_n = n_neighbor(j);
            
            m_upperRightQuarter = [m_upperRightQuarter, this_m]; 
            n_upperRightQuarter = [n_upperRightQuarter, this_n];
            
            
            disp(m_upperRightQuarter)
            disp(n_upperRightQuarter)
            
        end
        %}

        if (this_m == 1) && (this_n == ceil(R-1/2))
            notEnd = 0; 
        end 

    end 
       
end 

geomBdry = zeros(2*Ngrid+1);
for i2 = 1:size(m_upperRightQuarter)
    if (geomBdry(n(i2), m(i2))==0)
        geomBdry(n(i2), m(i2))=1;
    end
end 

plotGeomBdry(x,geomBdry)



%{
m = m_upperRightQuarter; n = n_upperRightQuarter; all_number_Bdry = number_upperRightQuarter; 

for k = 1:(number_upperRightQuarter-1)
    m_upperLeftQuarter = -m_upperRightQuarter(k); n_upperLeftQuarter = n_upperRightQuarter(k);

    all_number_Bdry = all_number_Bdry + 1; 

    m = [m, m_upperLeftQuarter]; 
    n = [n, n_upperLeftQuarter];

end

for l = 1:(number_upperRightQuarter-1)
    m_lowerRightQuarter = m_upperRightQuarter(k); n_lowerRightQuarter = -n_upperRightQuarter(k);

    all_number_Bdry = all_number_Bdry + 1; 

    m = [m, m_lowerRightQuarter]; 
    n = [n, n_lowerRightQuarter];
end

for ii = 1:(number_upperRightQuarter-1)
    m_lowerLeftQuarter = -m_upperRightQuarter(k); n_lowerLeftQuarter = -n_upperRightQuarter(k);

    all_number_Bdry = all_number_Bdry + 1; 

    m = [m, m_lowerLeftQuarter]; 
    n = [n, n_lowerLeftQuarter];
end

m = [m, m1, m2, m3, m4];
n = [n, n1, n2, n3, n4];
all_number_Bdry = all_number_Bdry + 3; 

geomBdry = zeros(2*Ngrid+1);
geomRadiusCollect = [];
number_insideCircle = 0; 

for jj = 1:all_number_Bdry
    if (geomBdry(n(jj), m(jj))==0)
        geomBdry(n(jj), m(jj))=1;
    end 

    Radius_geomBdry = sqrt((double(n(jj)) - grid0).^2 + (double(m(jj)) - grid0).^2); 
    geomRadiusCollect = [geomRadiusCollect, Radius_geomBdry];

    for kk=1:(2*Ngrid+1)
        for ll=1:(2*Ngrid+1)
            if (kk <= n(jj)) && (ll <= m(jj))
                number_insideCircle = number_insideCircle +1;
            end 
        end 
    end 
end


Ravg = mean(geomRadiusCollect);
Rstd = std(geomRadiusCollect);
Nint = number_insideCircle;

plotGeomBdry(x,geomBdry)
%}


function plotGeomBdry(x,GeomBdry)

imagesc(x,x,GeomBdry.'), axis square, axis xy
xlabel('x'), ylabel('y')
