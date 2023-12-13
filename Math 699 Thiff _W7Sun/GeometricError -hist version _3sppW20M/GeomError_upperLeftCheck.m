% calculate the geometric error _sppW20W 6/13
 
function [Rstd,Ravg,Nint] = GeometricError(R,flag)
 
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
 
if strcmpi(flag,'remove')
  % Remove extra parts of boundary.
  remove = true;
else
  remove = false;
end

if strcmpi(flag,'center0')
  % if the center is at lower left corner of the first grid 
  center0 = true;
else
  center0 = false;
end

%switch lower(flag)
%  case 'remove'
%    % Do something
%  case 'noremove'
%    % Do something else
%  case 'gaussian'
%    % Do something else
%  otherwise
%    error('Unknown flag %s',flag);
%end
 
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
    m3_ul = m3 - 0.5; n3_ul = n3 + 0.5; % m_upperLeft & n_upperLeft
    m3_lr = m3_ul + 1; n3_lr = n3_ul - 1;   % m_lowerRight & n_lowerRight
    
    dist3_ul = sqrt(m3_ul.^2 + n3_ul.^2);
     
    if (dist3_ul < R)
        n3 = n3 + 1;
        num_inside = num_inside + m3 - 1;
     
    elseif (dist3_ul == R)
        m3 = m3 - 1; n3 = n3 + 1; 
        num_inside = num_inside + m3 - 2;
     
    else
        m3 = m3 - 1;
        num_inside = num_inside - 1;
    end 
    
    if center0 
        if (m3 == 0.5) %&& (n3 == ceil(R))
            notEnd = 0; 
        end
    else
        if (m3 == 1) && (n3 == ceil(R-1/2))
            notEnd = 0; 
        end
    end 

      
    m_urq = [m_urq m3]; % all m in upper right corner of the circle 
    n_urq = [n_urq n3]; % all n in upper right corner of the circle 
end 


if center0
    m_urq = [m0 m_urq];
    n_urq = [n0 n_urq];
end 

if remove
    sz_temp = size(m_urq);
    all_d_temp = [];

    for i=1:sz_temp(2)
        % compute the distances between the boundary and the grid center to decide which one to delete 
        x_gridCenter = m_urq(i);
        y_gridCenter = n_urq(i);

        x_in = R./sqrt(1+ (y_gridCenter./x_gridCenter).^2);    % the intersection between the line y=(n/m)x and the circle 
        y_in = R./sqrt(1+ (x_gridCenter./y_gridCenter).^2);
        d_sq = (x_gridCenter - x_in).^2 + (y_gridCenter - y_in).^2;
        d = sqrt(d_sq);

        all_d_temp = [all_d_temp d];
    end 

    m_urqR = m_urq; % all m after being removed
    n_urqR = n_urq;

    chose = false;
    for i=1:sz_temp(2)
        
        mE_ul = m_urq(i)-0.5; nE_ul = n_urq(i)-0.5;
        mE_lr = mE_ul+1; nE_lr = nE_ul-1;
        
        distE_ul = sqrt(mE_ul.^2 + nE_ul.^2);
        distE_lr = sqrt(mE_lr.^2 + nE_lr.^2);
        
        outside = (distE_ul>=R) && (distE_lr>=R);
        inside = (distE_ul<R) && (distE_lr<R);
        
        if (i>=2) && (i<=(sz_temp(2)-1))
            diff_m = abs(m_urq(i-1) - m_urq(i+1));
            diff_n = abs(n_urq(i-1) - n_urq(i+1));   % ?
       
          if (outside || inside)  
              if (chose)
                  if all_d_temp(i-1) < all_d_temp(i)
                      if (m_urqR(i-1) ~= 0)
                        m_urqR(i) = 0;
                        n_urqR(i) = 0;
                      end
                  else
                      if (m_urqR(i-2) ~= 0)
                        m_urqR(i-1) = 0;
                        n_urqR(i-1) = 0;
                      end
                  end 
              end
              chose = true;
          else
              if (chose)
                  if (m_urqR(i-2) ~= 0)
                        m_urqR(i-1) = 0;
                        n_urqR(i-1) = 0;
                  end
              end 
              chose = false;
          
%                 if (m_urqR(i-1) == 0) && (n_urqR(i-1) == 0)         
%                   if all_d_temp(i-1) < all_d_temp(i)
%                       m_urqR(i) = 0;
%                       n_urqR(i) = 0;
% 
%                       m_urqR(i-1) = m_urq(i-1);
%                       n_urqR(i-1) = n_urq(i-1);
%                   end 
%                 end 
%               
%               if (m_urqR(i-1) ~= 0) && (n_urqR(i-1) ~= 0) && (diff_m <= 1) && (diff_n <= 1)
%                   m_urqR(i) = 0;
%                   n_urqR(i) = 0;
%               end
          end
        end 

    end 

    m_urq = nonzeros(m_urqR); n_urq = nonzeros(n_urqR);
    m_urq = m_urq'; n_urq = n_urq';
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