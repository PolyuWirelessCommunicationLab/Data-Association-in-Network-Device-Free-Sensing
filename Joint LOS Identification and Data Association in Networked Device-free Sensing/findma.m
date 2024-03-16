function q_ma=findma(q_possible,itera)
%select as many elements from q_possible as possible with the first itera
%elements will not share the same element
q_ma={};
if length(q_possible)==0
    return
end
q_matrix=(sortrows(q_possible.',1)).';
d_unique=unique(q_possible(1,:));
length_1_row=length(d_unique);
p_ma={};q_ma={};qn=0;
rep_num=zeros(1,length_1_row);
for k=1:length_1_row
    rep_num(k)=length(find(q_matrix(1,:)==d_unique(k)));
end
for i=1:rep_num(1)
    qn=qn+1;
    q_ma{qn}=q_matrix(:,i);
end 
p_matrix=q_matrix;

for k=2:length_1_row
    p_matrix=p_matrix(:,rep_num(k-1)+1:length(p_matrix(1,:)));
    num=length(q_ma);
    p_ma=q_ma;
    q_ma={};
    qn=0;
    for i=1:num
        dist_i=p_ma{i};
        for j=1:rep_num(k)
                flag=1;
                dist_j=p_matrix(:,j);
                for m=1:itera
                   flag=flag*(length(intersect(dist_j(m),dist_i(m,:)))==0||dist_j(m)==0);
                end
                if flag==1
                     qn=qn+1; 
                     q_ma{qn}=[dist_i,p_matrix(:,j)];
                end  
        end
    end
    if qn==0
      q_ma=p_ma;  
    end
end
end
