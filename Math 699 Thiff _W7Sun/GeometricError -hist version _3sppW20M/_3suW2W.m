%     m_urqR = m_urq; % all m in upper right corner of the circle with several grids removed
%     n_urqR = n_urq;



%       m1_ul = m_urq(i-2)-0.5; n1_ul = n_urq(i-2)-0.5;
%       m2_ul = m_urq(i-1)-0.5; n2_ul = n_urq(i-1)-0.5;
% 
%       m1_lr = m1_ul+1; n1_lr = n1_ul-1; 
%       m2_lr = m2_ul+1; n2_lr = n2_ul-1;
% 
%       m_test = [m1_ul m1_lr m2_ul m2_lr m3_ul m3_lr];  % [m1_ul m1_lr m2_ul m2_lr m3_ul m3_lr]
%       n_test = [n1_ul n1_lr n2_ul n2_lr n3_ul n3_lr];
% 
%       sz_test = size(m_test); 
%       dist_test = [];   
%       
% 
%       for j=1:sz_test(2)
%           dist = sqrt(m_test(j).^2 + n_test(j).^2);
%           dist_test = [dist_test dist]; % [dist1_ul dist1_lr dist2_ul dist2_lr dist3_ul dist3_lr]
%       end 



    
%     outside = dist3_ul > R && dist3_lr > R;
%     inside = dist3_ul < R && dist3_lr < R;
    
    % decide which one to delete
    
%     if (sz_urq(2)>=3)
%         last_m = m_upperRightQuarter(sz_urq(2)-1);
%         last_n = n_upperRightQuarter(sz_urq(2)-1);
% 
%         m_ul_last = last_m - 0.5; n_ul_last = last_n + 0.5;  
%         m_lr_last = m_ul_last + 1; n_lr_last = n_ul_last - 1;
% 
%         dist_ul_last = sqrt(m_ul_last.^2 + n_ul_last.^2);
%         dist_lr_last = sqrt(m_lr_last.^2 + n_lr_last.^2); 
% 
%         outside_last = dist_ul_last > R && dist_lr_last > R;
%         inside_last = dist_ul_last < R && dist_lr_last < R;
%         
%         two_m = [m3 last_m];
%         two_n = [n3 last_n];
%         two_d = [];
% 
%         if (outside || inside) && (outside_last || inside_last)
%             for i=1:2
%                 x_gridCenter = two_m(i);
%                 y_gridCenter = two_n(i);
% 
%                 x_in = R./sqrt(1+ (y_gridCenter./x_gridCenter).^2);    % the intersection between the line y=(n/m)x and the circle 
%                 y_in = R./sqrt(1+ (x_gridCenter./y_gridCenter).^2);
%                 d_sq = (x_gridCenter - x_in).^2 + (y_gridCenter - y_in).^2;
%                 d = sqrt(d_sq);
%                 
%                 two_d = [two_d d];
%             end 
% 
%         if remove 
%             if (two_d(1)>=two_d(2))
%                 m_upperRightQuarter(sz_urq(2)) = []; 
%                 n_upperRightQuarter(sz_urq(2)) = [];
%                 % num_insdie doesn't change 
%             else 
%                 m_upperRightQuarter(sz_urq(2)-1) = []; 
%                 n_upperRightQuarter(sz_urq(2)-1) = [];
%                 num_inside = num_inside + 1; 
%             end 
%         end 


    %     
    %     if ~donotadd || ~remove

    %     end
    
% end 
% end 
           

  



     
     tf_m = ismember(m_urq, m_urqR);  % [for R=3, [1 1 1 1 1]]
     tf_n = ismember(m_urq, m_urqR);  % [for R=3]
     
     
     
         %        dist_ul = [all_ul(i) all_ul(i-1) all_ul(i-2)];
    %        dist_lr = [all_lr(i) all_lr(i-1) all_lr(i-2)];
    % 
    %       grid1A = (dist_ul(1)>R) && ((dist_lr(1)>R)||(dist_lr(1)<R));
    %       grid2A = (dist_ul(2)>R) && (dist_lr(2)>R);
    %       grid3A = (dist_lr(3)>R) && ((dist_ul(3)>R)||(dist_ul(3)<R));
    %       shapeA = grid1A && grid2A && grid3A;
    %       
    %       shapeA = grid2A; 
    % 
    %       grid1B = (dist_ul(1)<R) && ((dist_lr(1)>R)||(dist_lr(1)<R));
    %       grid2B = (dist_ul(2)<R) && (dist_lr(2)<R);
    %       grid3B = (dist_lr(3)<R) && ((dist_ul(3)>R)||(dist_ul(3)<R));
    %       shapeB = grid1B && grid2B && grid3B;
    % 
    %       shapeB = grid2B;