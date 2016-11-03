function [ext_B] = ABCR_extrapolation(B, points, extra_p, k)

ext_B = zeros(size(B,1),size(B,2)+1);
for j=1:size(B,1)
    w = B(j,:);   
    ind = find(w);
    flag = 0;
    if (size(ind,2) < size(w,2))
        if (extra_p < points(1) && abs(w(1)) > 0)
            flag = 1;
        elseif (extra_p > points(end)  && abs(w(end)) > 0)
            flag = 1;
        else                
            ind_out = ~w;
            dis_in =min(abs(points(ind)-extra_p));
            dis_out =min(abs(points(ind_out)-extra_p));
             if dis_in<dis_out   
                 flag =1;
             end       
        end
    else
        flag = 1;
    end
    if flag > 0
       mat=[];
       improve_vdm_mat=[];
       half_size = floor((ind(end)-ind(1)+1)/2);
       part1 = ind(1):(ind(1)+half_size -1);
       part2 = (ind(1)+half_size):ind(end);
       
       dis_l =min(abs(points(part1)-extra_p));
       dis_r =min(abs(points(part2)-extra_p));
       mid=abs(points(part1(end))+points(part2(1)))/2;
       if (abs(extra_p-mid) < 0.001)
            w_local =  w(ind);
            part_p = points(ind);
       elseif (dis_l < dis_r)           
            w_local =  w(part1);
            part_p = points(part1);
       else
             w_local =  w(part2);
             part_p = points(part2);
       end;
         
       m = abs(part_p(1) + part_p(end))/2;
       s  = abs(part_p(1) - part_p(end))/2;
      
       for i=1:k
           mat(:,i) = part_p.^(i-1);
           improve_vdm_mat(:,i) = ((part_p-m)/s).^(i-1);
       end
       cond(mat);
       cond(improve_vdm_mat);

       [v ,z, u] = svd(improve_vdm_mat,0);
       err =  norm(w_local'-v*v'*w_local');

       r= v'* w_local';  % coeffiecients of the orthogonal representation

       err = norm(improve_vdm_mat*r - w_local');
     
       for i=1:k
           mat(size(part_p,2)+1,i) = extra_p.^(i-1);
           improve_vdm_mat(size(part_p,2)+1,i) = ((extra_p-m)/s).^(i-1);
       end

       % improve_vdm_mat = v*z*u'

         
       v = improve_vdm_mat*u*inv(z);

       w_new = v*r;

       ext_B(j,end) = w_new(end);
    end
    ext_B(j,1:end-1) = B(j,:);   
end
end

