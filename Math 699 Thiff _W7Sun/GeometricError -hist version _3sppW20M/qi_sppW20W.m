            number_temp = number_temp + 1;
            m_upperRightQuarter = [m_upperRightQuarter, this_m]; 
            n_upperRightQuarter = [n_upperRightQuarter, this_n];



        if (size(all_this_m)==1)
            this_m = all_this_m(1); this_n = all_this_n(1);
        end 
        
        if (size(all_this_m)>1)
            if ((all_this_m(1)>all_this_m(2))&&(all_this_n(1)==all_this_n(2)))||((all_this_m(1)==all_this_m(2))&&(all_this_n(1)<all_this_n(2)))
                this_m = all_this_m(1); this_n = all_this_n(1);
            else 
                this_m = all_this_m(1); this_n = all_this_n(1);
            end 
        end 