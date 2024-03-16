function p_m=tar_det(ind_BS)
% this is to find all mappings satisfying the contraints (21) and (22)
global D_sum dist_m_0
p_m=[];
if length(ind_BS)==2
  p_m=D_sum{ind_BS(1),ind_BS(2)}.';
  return
end
for i=1:length(ind_BS)
    dist_m_0_i=dist_m_0{ind_BS(i)};
    ind_BS_i=setdiff(ind_BS,ind_BS(i));
    p_m_i=tar_det(ind_BS_i);
    for j=1:size(p_m_i,2)
        for l=1:length(dist_m_0_i)
            flag=1;
            for n=1:length(ind_BS_i)
                D_12=D_sum{ind_BS(i),ind_BS_i(n)};
                if length(D_12)==0
                    ind_12=[];
                else
                    ind_12=intersect(find(D_12(:,1)==dist_m_0_i(l)),find(D_12(:,2)==p_m_i(n,j)));
                end
                flag = flag*sign(length(ind_12));
                if flag==0
                    break
                end
            end
            if flag>0
                p_m=[p_m,[p_m_i(1:i-1,j);dist_m_0_i(l);p_m_i(i:length(ind_BS_i),j)]];
            end
        end
    end
end
end