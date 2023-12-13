% calculate the geometric error _sppW20W 6/13

% GeometricError(R,'nonremove') gives a plot of a discretized 
% circle and the standard deviation of the grids on the boundary. The
% boundary of the drawn circle has to pass the grids consisting a discretized
% one.

% GeometricError(R,'remove') gives a plot of a discretized 
% circle and the standard deviation of the grids on the boundary. The
% boundary of the drawn circle does not have to pass the grids consisting a 
% discretized one.
 
function [Rstd,Ravg,Nint] = GeometricError(R,flag)
 
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

if strcmpi(flag,'nonremove') %_3suW1
  % Remove extra parts of boundary.
  nonremove = true;
else
  nonremove = false;
end
 
if strcmpi(flag,'remove') %_3suW2R
  % Remove extra parts of boundary.
  remove = true;
else
  remove = false;
end

if strcmpi(flag,'remove_previous') %_3suW2
  % Remove extra parts of boundary.
  remove_previous = true;
else
  remove_previous = false;
end

if strcmpi(flag,'center0') %_3suW2
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
 

if nonremove
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
end

if remove  % _3suW2R 6/28
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
end 


if center0
    m_urq = [m0 m_urq];
    n_urq = [n0 n_urq];
end 

if remove_previous  % _3suW2
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
all_radius = [];
% disp(sz_urq(2))

% find an upper bound for d^2 _3suW3R
m_urq1 = []; n_urq1 = []; 
m_urq2 = []; n_urq2 = []; 
m_urq3 = []; n_urq3 = []; 
m_urq4 = []; n_urq4 = []; % put m, n into 4 cases

m_oct = []; n_oct = [];
m_oct1 = []; n_oct1 = [];
m_oct2 = []; n_oct2 = []; 
m_oct4 = []; n_oct4 = []; % in the first octant _3suW6F(Sh) 8/3

m_1a = []; n_1a = []; m_1b = []; n_1b = [];
m_2a = []; n_2a = []; m_2b = []; n_2b = [];
m_4a = []; n_4a = []; m_4b = []; n_4b = []; % 1a and 4b _3suW8F 8/10(sh)

all_ub1 = []; all_lb1 = []; 
all_ub2 = []; all_lb2 = []; 
all_ub3 = []; all_lb3 = []; 
all_ub4 = []; all_lb4 = []; % arrays which include all lower and upper bounds for 4 cases

if mod(sz_urq(2),2)==0  % _3suW6F
    sz_oct = sz_urq(2)/2;
else
    sz_oct = (sz_urq(2)-1)/2;
    md = sz_oct + 1; % add a test for m_urq(md)
    
    x_Amd = sqrt(R^2 - (n_urq(md)-1/2).^2); y_Bmd = sqrt(R^2 - (m_urq(md)-1/2).^2); 
    y_Cmd = sqrt(R^2 - (m_urq(md)+1/2).^2); x_Dmd = sqrt(R^2 - (n_urq(md)+1/2).^2);
    % y_Emd = sqrt(R^2 - (m_urq(md)+1/2).^2); y_Fmd = sqrt(R^2 - (m_urq(md)-1/2).^2);
    x_Gmd = sqrt(R^2 - (n_urq(md)-1/2).^2); x_Hmd = sqrt(R^2 - (n_urq(md)+1/2).^2);
    
    if ((m_urq(md)-1/2) <= x_Amd) && ((m_urq(md)+1/2) >= x_Amd) && ((n_urq(md)-1/2) <= y_Bmd) && ((n_urq(md)+1/2) >= y_Bmd)
        type_md = 1;
    end 
    
    if (((m_urq(md)-1/2) <= x_Dmd) && ((m_urq(md)+1/2) > x_Dmd) && ((n_urq(md)-1/2) < y_Cmd) && ((n_urq(md)+1/2) > y_Cmd)) || (((m_urq(md)-1/2) < x_Dmd) && ((m_urq(md)+1/2) > x_Dmd) && ((n_urq(md)-1/2) <= y_Cmd) && ((n_urq(md)+1/2) > y_Cmd))
        type_md = 2;
    end 
    
    if ((m_urq(md)-1/2) <= x_Gmd) && ((m_urq(md)+1/2) >= x_Gmd) && ((m_urq(md)-1/2) <= x_Hmd) && ((m_urq(md)+1/2) >= x_Hmd)
        type_md = 4;
    end
end

for i=1:sz_urq(2)
    % compute the average radius and the standard deviation _3suW1?
    radius = sqrt(m_urq(i).^2 + n_urq(i).^2); 
    all_radius = [all_radius, radius];

    d = abs(R - radius);
    all_d = [all_d d];

    % put m, n into 4 cases  
    x_A = sqrt(R^2 - (n_urq(i)-1/2)^2); y_B = sqrt(R^2 - (m_urq(i)-1/2)^2);
    y_C = sqrt(R^2 - (m_urq(i)+1/2)^2); x_D = sqrt(R^2 - (n_urq(i)+1/2)^2);
    y_E = sqrt(R^2 - (m_urq(i)+1/2)^2); y_F = sqrt(R^2 - (m_urq(i)-1/2)^2);
    x_G = sqrt(R^2 - (n_urq(i)-1/2)^2); x_H = sqrt(R^2 - (n_urq(i)+1/2)^2);
    
    % for type 1 in the upper right quarter 
    if ((m_urq(i)-1/2) <= x_A) && ((m_urq(i)+1/2) >= x_A) && ((n_urq(i)-1/2) <= y_B) && ((n_urq(i)+1/2) >= y_B)
        m_urq1 = [m_urq1 m_urq(i)];
        n_urq1 = [n_urq1 n_urq(i)];
        
        % type 1 in the 1st octant
        if i<=sz_oct
            m_oct1 = [m_oct1 m_urq(i)];
            n_oct1 = [n_oct1 n_urq(i)];
            
            % type 1a & 1b for remove case
            if remove
                if i>1
                   if (m_urq(i-1)==m_urq(i)) && (n_urq(i-1)+1==n_urq(i)) % 1a and 4b _3suW8F 8/10(sh)
                        m_1a = [m_1a m_urq(i)];
                        n_1a = [n_1a n_urq(i)];
                   end

                   if (m_urq(i-1)-1==m_urq(i)) && (n_urq(i-1)+1==n_urq(i))
                        m_1b = [m_1b m_urq(i)];
                        n_1b = [n_1b n_urq(i)];
                   end
                end
            end

        end
    end 
    
    % type 2 in the upper right quarter
    if (((m_urq(i)-1/2) <= x_D) && ((m_urq(i)+1/2) > x_D) && ((n_urq(i)-1/2) < y_C) && ((n_urq(i)+1/2) > y_C)) || (((m_urq(i)-1/2) < x_D) && ((m_urq(i)+1/2) > x_D) && ((n_urq(i)-1/2) <= y_C) && ((n_urq(i)+1/2) > y_C))
        m_urq2 = [m_urq2 m_urq(i)];
        n_urq2 = [n_urq2 n_urq(i)];
        
        % type 2 in the 1st octant        
        if i<=sz_oct
            m_oct2 = [m_oct2 m_urq(i)];
            n_oct4 = [n_oct4 n_urq(i)];
           
           % type 2b for remove case
           if remove
               if i>1 
                   if (m_urq(i-1)-1==m_urq(i)) && (n_urq(i-1)==n_urq(i))
                        m_2a = [m_2a m_urq(i)];
                        n_2a = [n_2a n_urq(i)];
                   end
                   if (m_urq(i-1)-1==m_urq(i)) && (n_urq(i-1)+1==n_urq(i))
                        m_2b = [m_2b m_urq(i)];
                        n_2b = [n_2b n_urq(i)];
                   end
               end
           end
           
        end
    end 
    
    % type 3 in the upper right quarter
    if ((n_urq(i)-1/2) <= y_E) && ((n_urq(i)+1/2) >= y_E) && ((n_urq(i)-1/2) <= y_F) && ((n_urq(i)+1/2) >= y_F)
        m_urq3 = [m_urq3 m_urq(i)];
        n_urq3 = [n_urq3 n_urq(i)];
    end 
    
    % type 4 in the upper right quarter
    if ((m_urq(i)-1/2) <= x_G) && ((m_urq(i)+1/2) >= x_G) && ((m_urq(i)-1/2) <= x_H) && ((m_urq(i)+1/2) >= x_H)
        m_urq4 = [m_urq4 m_urq(i)];
        n_urq4 = [n_urq4 n_urq(i)];
        
        % type 4 in the 1st octant
        if i<=sz_oct
            m_oct4 = [m_oct4 m_urq(i)];
            n_oct4 = [n_oct4 n_urq(i)];
            
            % type 4a & 4b for remove case
            if remove
                if i>1
                    if (m_urq(i-1)==m_urq(i)) && (n_urq(i-1)+1==n_urq(i))
                        m_4a = [m_4a m_urq(i)];
                        n_4a = [n_4a n_urq(i)];
                    end

                    if (m_urq(i-1)-1==m_urq(i)) && (n_urq(i-1)+1==n_urq(i)) % 1a and 4b _3suW8F 8/10(sh)
                        m_4b = [m_4b m_urq(i)];
                        n_4b = [n_4b n_urq(i)];
                    end
                end
            end 
            
        end
        
    end 
    
end 

if mod(sz_urq(2),2)~=0
    if type_md==1
        m_oct1 = [m_oct1 m_urq(md)];
        n_oct1 = [n_oct1 n_urq(md)];
    elseif type_md==2
        m_oct2 = [m_oct2 m_urq(md)];
        n_oct2 = [n_oct2 n_urq(md)];
    else
        m_oct4 = [m_oct2 m_urq(md)];
        n_oct4 = [n_oct4 n_urq(md)];
    end 

end

sz1 = size(m_urq1); sz2 = size(m_urq2); sz3 = size(m_urq3); sz4 = size(m_urq4); 
sz_tot = sz1(2)+sz2(2)+sz3(2)+sz4(2); 
frac1 = sz1(2)/sz_tot;
frac2 = sz2(2)/sz_tot;
frac3 = sz3(2)/sz_tot;
frac4 = sz4(2)/sz_tot;

disp('in the upper right quarter:')
disp("total in the upper right quarter = " +sz_tot+";")
disp("#1 = "+sz1(2)+"; fraction 1 = "+frac1+";")
disp("#2 = "+sz2(2)+"; fraction 2 = "+frac2+";")
disp("#3 = "+sz3(2)+"; fraction 3 = "+frac3+";")
disp("#4 = "+sz4(2)+"; fraction 4 = "+frac4+";")

sz1o = size(m_oct1); sz2o = size(m_oct2); sz4o = size(m_oct4); 
sz_tot_oct = sz1o(2)+sz2o(2)+sz4o(2); 
frac1o = sz1o(2)/sz_tot_oct;
frac2o = sz2o(2)/sz_tot_oct;
frac4o = sz4o(2)/sz_tot_oct;

disp('in the first octant:')
disp("total_octant = " + sz_tot_oct)
disp("#1 = "+sz1o(2)+"; fraction 1 = "+frac1o+";")
disp("#2 = "+sz2o(2)+"; fraction 2 = "+frac2o+";")
disp("#4 = "+sz4o(2)+"; fraction 4 = "+frac4o+";")

sz1a = size(m_1a); sz1b = size(m_1b);
disp("#1a = "+sz1a(2)+"; #1b = "+sz1b(2)+"; ")

sz2a = size(m_2a); sz2b = size(m_2b);  
disp("#2a = "+sz2a(2)+";   #2b = "+sz2b(2)+"; ")

sz4a_temp = size(m_4a); sz4a = sz4a_temp(2) + 1; sz4b = size(m_4b); % the 1st one is not counted in the loop
disp("#4a = "+sz4a+"; #4b = "+sz4b(2)+"; ") 


if center0 
    all_d = [all_d all_d all_d all_d];
else 
    d0 = ceil(R-1/2)-R;
    all_d = [all_d all_d all_d all_d d0 d0 d0 d0];    % all distances between grid center and the boundary 
end 

d_avg = mean(all_d);
d_sq_avg = d_avg.^2;    % square of average of d
% disp("d_avg^2 = " + d_sq_avg)
    
 
num_inside = 4*num_inside + 4*(ceil(R-1/2)-1) + 1; % grids in one quarter*4 + grids in the center + origin
 

if center0 
    all_radius = [all_radius all_radius all_radius all_radius];
else
    r0 = m0;
    all_radius = [all_radius r0 r0 r0 r0];
end 
  
 
Ravg = mean(all_radius);
Rstd = std(all_radius);
Nint = num_inside;







% plot _3suW1
Ngrid = ceil(R-1/2)+1;  % size of grid quadrant
grid0 = Ngrid;     % center of grid
x = -Ngrid:Ngrid;
geomBdry = sparse(2*Ngrid);
 
if center0
    for jj = 1:sz_urq(2)
        %if (geomBdry(n_upperRightQuarter(jj), m_upperRightQuarter(jj))==0)
            geomBdry(n_urq(jj)+0.5+grid0, m_urq(jj)+0.5+grid0)=1;
            geomBdry(-n_urq(jj)+0.5+grid0, m_urq(jj)+0.5+grid0)=1;
            geomBdry(n_urq(jj)+0.5+grid0, -m_urq(jj)+0.5+grid0)=1;
            geomBdry(-n_urq(jj)+0.5+grid0, -m_urq(jj)+0.5+grid0)=1;
        %end 
    end 

else 
    for jj = 1:sz_urq(2)
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
 
 
function plotGeomBdry(x,GeomBdry)
 
% imagesc(x,x,GeomBdry.')
% n = (size(GeomBdry))/2;
% x = [-n/2:n/2];
% axis equal

spy(GeomBdry), axis square, axis xy
xlabel('x'), ylabel('y')