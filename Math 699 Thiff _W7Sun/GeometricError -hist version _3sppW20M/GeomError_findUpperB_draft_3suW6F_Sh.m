    
    % compute the average d^2 _3suW2
%     x_gridCenter = m_urq(i);
%     y_gridCenter = n_urq(i);
%     
%     x_in = R./sqrt(1+ (y_gridCenter./x_gridCenter).^2);    % the intersection between the line y=(n/m)x and the circle 
%     y_in = R./sqrt(1+ (x_gridCenter./y_gridCenter).^2);
%     d_sq = (x_gridCenter - x_in).^2 + (y_gridCenter - y_in).^2;



% find the upper bounds for d: lower bound for R - sqrt(m_urq(i)^2+n_urq(i)^2) <= d <= upper bound for R - sqrt(m_urq(i)^2+n_urq(i)^2)
for i=1:sz1(2)
    ub1 = abs(min(sqrt((m_urq(i)+1/2)^2+(n_urq(i)-1/2)^2),sqrt((m_urq(i)-1/2)^2+(n_urq(i)+1/2)^2))-sqrt(m_urq(i)^2+n_urq(i)^2)); 
    lb1 = abs(sqrt((m_urq(i)-1/2)^2+(n_urq(i)-1/2)^2)-sqrt(m_urq(i)^2+n_urq(i)^2));
    
    ub1 = max(ub1, lb1); lb1 = min(ub1, lb1);
    all_ub1 = [all_ub1 ub1]; all_lb1 = [all_lb1 lb1];
end

for i=1:sz2(2)
    ub2 = abs(sqrt((m_urq(i)+1/2)^2+(n_urq(i)+1/2)^2)-sqrt(m_urq(i)^2+n_urq(i)^2)); 
    lb2 = abs(max(sqrt((m_urq(i)+1/2)^2+(n_urq(i)-1/2)^2),sqrt((m_urq(i)-1/2)^2+(n_urq(i)+1/2)^2))-sqrt(m_urq(i)^2+n_urq(i)^2));
    
    ub2 = max(ub2, lb2); lb2 = min(ub2, lb2);
    all_ub2 = [all_ub2 ub2]; all_lb2 = [all_lb2 lb2];
end

for i=1:sz3(2)
    ub3 = abs(sqrt((m_urq(i)-1/2)^2+(n_urq(i)+1/2)^2)-sqrt(m_urq(i)^2+n_urq(i)^2));
    lb3 = abs(sqrt((m_urq(i)+1/2)^2+(n_urq(i)-1/2)^2)-sqrt(m_urq(i)^2+n_urq(i)^2));
    
    ub3 = max(ub3, lb3); lb3 = min(ub3, lb3);
    all_ub3 = [all_ub3 ub3]; all_lb3 = [all_lb3 lb3];
end

for i=1:sz4(2)
    ub4 = abs(sqrt((m_urq(i)+1/2)^2+(n_urq(i)-1/2)^2)-sqrt(m_urq(i)^2+n_urq(i)^2));
    lb4 = abs(sqrt((m_urq(i)-1/2)^2+(n_urq(i)+1/2)^2)-sqrt(m_urq(i)^2+n_urq(i)^2));
    
    ub4 = max(ub4, lb4); lb4 = min(ub4, lb4);
    all_ub4 = [all_ub4 ub4]; all_lb4 = [all_lb4 lb4];
end

avg_ub1 = mean(all_ub1); avg_lb1 = mean(all_lb1);
avg_ub2 = mean(all_ub2); avg_lb2 = mean(all_lb2);
avg_ub3 = mean(all_ub3); avg_lb3 = mean(all_lb3);
avg_ub4 = mean(all_ub4); avg_lb4 = mean(all_lb4);

all_avg_ub = [all_ub1 all_ub2 all_ub3 all_ub4];
all_avg_lb = [all_lb1 all_lb2 all_lb3 all_lb4];

ub = max(all_avg_ub);
lb = min(all_avg_lb);
ub_sq = ub^2;
lb_sq = lb^2;
%disp(lb_sq + " <= d_avg^2 <= " + ub_sq)


% for i=1:sz_oct(2)
%     if ((m_oct(i)-1/2) <= x_A) && ((m_oct(i)+1/2) >= x_A) && ((n_oct(i)-1/2) <= y_B) && ((n_oct(i)+1/2) >= y_B)
%         m_urq1 = [m_urq1 m_urq(i)];
%         n_urq1 = [n_urq1 n_urq(i)];
%     end 
%     
%     if (((m_oct(i)-1/2) <= x_D) && ((m_oct(i)+1/2) > x_D) && ((n_oct(i)-1/2) < y_C) && ((n_oct(i)+1/2) > y_C)) || (((m_oct(i)-1/2) < x_D) && ((m_urq(i)+1/2) > x_D) && ((n_urq(i)-1/2) <= y_C) && ((n_urq(i)+1/2) > y_C))
%         m_urq2 = [m_urq2 m_urq(i)];
%         n_urq2 = [n_urq2 n_urq(i)];
%     end 
%     
%     if ((m_urq(i)-1/2) <= x_G) && ((m_urq(i)+1/2) >= x_G) && ((m_urq(i)-1/2) <= x_H) && ((m_urq(i)+1/2) >= x_H)
%         m_urq4 = [m_urq4 m_urq(i)];
%         n_urq4 = [n_urq4 n_urq(i)];
%     end 
% end
