% another plot tried _3suW2R 6/28
 
function [Rstd,Ravg,Nint] = GeometricError_remove(R,flag)
 
if nargin < 1, R = 10; end
if nargin < 2, flag = 'noremove'; end
 
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
 

if strcmpi(flag,'center0')
  % if the center is at lower left corner of the first grid 
  center0 = true;
else
  center0 = false;
end

if center0 
    m0 = ceil(R)-0.5; n0 = 0.5; 
else 
    m0 = ceil(R-1/2); n0 = 0; 
end 
 
% find the boundary points on the upperRightQuarter
m_urq = []; % horizontal value of the coordinate
n_urq = []; % vertical value of the coordinate  
num_inside = 0; 

all_d = [];
 
m3 = m0; n3 = n0; 
notEnd = 1; 
 
while (notEnd == 1) 
    mU = m3; nU = n3+1; 
    mDiag = m3 - 1; nDiag = n3 + 1;
    mL = m3 - 1; nL = n3;  
    
    mNeib = [mU mDiag mL]; nNeib = [nU nDiag nL];
    d_temp = [];
    
    % compute the distances between the boundary and the grid center to decide which one to delete
    for i=1:3
        x_gridCenter = mNeib(i);
        y_gridCenter = nNeib(i);

        x_in = R./sqrt(1+ (y_gridCenter./x_gridCenter).^2);    % the intersection between the line y=(n/m)x and the circle 
        y_in = R./sqrt(1+ (x_gridCenter./y_gridCenter).^2);
        d_sq = (x_gridCenter - x_in).^2 + (y_gridCenter - y_in).^2;
        d = sqrt(d_sq);

        d_temp = [d_temp d];
    end 
    
    [M, I] = min(d_temp);
    m3 = mNeib(I); n3 = nNeib(I);
    
    if center0 
        if (m3 == 0.5) %&& (n3 == ceil(R))
            notEnd = 0; 
        end
    else
        if (m3 == 1) %&& (n3 == ceil(R-1/2))
            notEnd = 0; 
        end
    end 

    m_urq = [m_urq m3]; n_urq = [n_urq n3];
end 

if center0
    m_urq = [m0 m_urq];
    n_urq = [n0 n_urq];
end 


sz_urq = size(m_urq);

for i=1:sz_urq(2)
    x_gridCenter = m_urq(i);
    y_gridCenter = n_urq(i);
    
    x_in = R./sqrt(1+ (y_gridCenter./x_gridCenter).^2);    % the intersection between the line y=(n/m)x and the circle 
    y_in = R./sqrt(1+ (x_gridCenter./y_gridCenter).^2);
    d_sq = (x_gridCenter - x_in).^2 + (y_gridCenter - y_in).^2;
    d = sqrt(d_sq);
    
    all_d = [all_d d];
end 

if center0 
    all_d = [all_d all_d all_d all_d];
else 
    d0 = ceil(R-1/2)-R;
    all_d = [all_d all_d all_d all_d d0 d0 d0 d0];    % all distances between grid center and the boundary 
end 

d_avg = mean(all_d);
d_sq_avg = d_avg.^2;    % square of average of d
disp("The square of the average of all the distances is " + d_sq_avg)
    
 
num_inside = 4*num_inside + 4*(ceil(R-1/2)-1) + 1; % grids in one quarter*4 + grids in the center + origin
 
 
Ngrid = ceil(R-1/2)+1;  % size of grid quadrant
grid0 = Ngrid;     % center of grid
x = -Ngrid:Ngrid;
geomBdry = sparse(2*Ngrid);

sz = size(m_urq); 
 
if center0
    for jj = 1:sz(2)
        %if (geomBdry(n_upperRightQuarter(jj), m_upperRightQuarter(jj))==0)
            geomBdry(n_urq(jj)+0.5+grid0, m_urq(jj)+0.5+grid0)=1;
            geomBdry(-n_urq(jj)+0.5+grid0, m_urq(jj)+0.5+grid0)=1;
            geomBdry(n_urq(jj)+0.5+grid0, -m_urq(jj)+0.5+grid0)=1;
            geomBdry(-n_urq(jj)+0.5+grid0, -m_urq(jj)+0.5+grid0)=1;
        %end 
    end 

else 
    for jj = 1:sz(2)
        %if (geomBdry(n_upperRightQuarter(jj), m_upperRightQuarter(jj))==0)
            geomBdry(n_urq(jj)+grid0, m_urq(jj)+grid0)=1;
            geomBdry(-n_urq(jj)+grid0, m_urq(jj)+grid0)=1;
            geomBdry(n_urq(jj)+grid0, -m_urq(jj)+grid0)=1;
            geomBdry(-n_urq(jj)+grid0, -m_urq(jj)+grid0)=1;
        %end 
    end       
            geomBdry(grid0 + ceil(R-1/2), grid0)=1;
            geomBdry(grid0 - ceil(R-1/2), grid0)=1;
            geomBdry(grid0, grid0 + ceil(R-1/2))=1;
            geomBdry(grid0, grid0 - ceil(R-1/2))=1;
end 
 
plotGeomBdry(x,geomBdry)
 
 
all_radius = [];
 
for k=1:sz(2)
    radius = sqrt(m_urq(k).^2 + n_urq(k).^2); 
    all_radius = [all_radius, radius];
end 

if center0 
    all_radius = [all_radius all_radius all_radius all_radius];
else
    r0 = m0;
    all_radius = [all_radius r0 r0 r0 r0];
end 
  
 
Ravg = mean(all_radius);
Rstd = std(all_radius);
Nint = num_inside;
 
 
function plotGeomBdry(x,GeomBdry)
 
% imagesc(x,x,GeomBdry.')
% n = (size(GeomBdry))/2;
% x = [-n/2:n/2];
% axis equal

spy(GeomBdry), axis square, axis xy
xlabel('x'), ylabel('y')