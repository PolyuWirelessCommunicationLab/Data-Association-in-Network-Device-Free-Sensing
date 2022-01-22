global coordinate_bs dist_m M K R dist_bs delta
M=4; % the number of base stations
R=240; % the region of interest,a square
gt_coor={};
error_pro=zeros(3,6);
error_pro_da=zeros(1,7);
% error_dist_da=zeros(1,7);
delta_d=1.5;
sigma_d=delta_d/2;
 for K=2:7% the number of target
iteration_max=100000;
error_num_da=0;
% error_dist=0;
error_num1=0;
error_num2=0;
error_num3=0;
gt_prob=0;
iteration=1; 
diff_tg=0;
real_itera=0;
while iteration<iteration_max

%produce uniformly coordinates of base stations and targets
coordinate_bs=R*(rand(M,2)-0.5);
coordinate_tg=R*(rand(K,2)-0.5);

%calculate distances from base stations to targets
dist_r=zeros(M,K); 
d=zeros(M,K);
for i=1:M
    for j=1:K
       dist_r(i,j)=sqrt((coordinate_bs(i,1)-coordinate_tg(j,1))^2+(coordinate_bs(i,2)-coordinate_tg(j,2))^2);
    end
end

% calculate the distance estimated in Phase 1 
dist_m=zeros(M,K);
for m=1:M
   for k=1:K
       dist_m(m,k)=ceil(dist_r(m,k)/delta_d)*delta_d-delta_d/2; 
   end
end
% add noise to the accurate distances for true range model
% dist_m=abs(dist_r+sigma*randn(M,K)); 

% judge whether there are more one target in the same range
same_range=0;
for m=1:M
    if length(unique(dist_m(m,:)))<K
        same_range=1;
    end
end
if same_range==1
   continue
end
% caculate the distances between any two base stations
dist_bs=zeros(M,M);
for i=1:M
    for j=1:M
        dist_bs(i,j)=sqrt((coordinate_bs(i,1)-coordinate_bs(j,1))^2+(coordinate_bs(i,2)-coordinate_bs(j,2))^2);
    end

end


% select the distances satisfying triangular inequality
delta=3*sigma_d;
p={};
t=0;
for i=1:K
     for j=1:K
         for k=1:K
             if abs(dist_m(1,i)-dist_m(2,j))<dist_bs(1,2)+delta&&abs(dist_m(1,i)+dist_m(2,j))>dist_bs(1,2)-delta
                  if  abs(dist_m(2,j)-dist_m(3,k))<dist_bs(2,3)+delta&&abs(dist_m(2,j)+dist_m(3,k))>dist_bs(2,3)-delta
                      if  abs(dist_m(1,i)-dist_m(3,k))<dist_bs(1,3)+delta&&abs(dist_m(1,i)+dist_m(3,k))>dist_bs(1,3)-delta
                          t=t+1;
                          p(t)={[dist_m(1,i);dist_m(2,j);dist_m(3,k)]};
                      end
                  end
             end
         end
     end
end

i=3;
T=t;
q3={};
t=1;
q3_num=0;
%evaluate whether the mathcing is feasible 
while t<=T
         initial_dist=p{t};
         t=t+1;
         [dist_error,initial_tg]=coorest(initial_dist);
         if sum(dist_error>1.5*sigma_d+0.01)>0
            continue 
         end
         q3_num=q3_num+1;
         q3{q3_num}=initial_dist;
end
%convert cell to matrix
l_q3=q3_num;
q3_matrix=zeros(3,l_q3);
for i=1:l_q3
    q3_matrix(:,i)=q3{i};
end    

%get the consistent matrix for the three BSs
q_ma={}; 
qn=0;
rep_num=zeros(1,K);
for k=1:K
    rep_num(k)=length(find(q3_matrix(1,:)==dist_m(1,k)));
end
for i=1:rep_num(1)
    qn=qn+1;
    q_ma{qn}=q3_matrix(:,i);
end 
p3_matrix=q3_matrix;
dist3_msrt=sort(dist_m(1:3,:),2);
%find the mathcing statisfying (17) in the set (36)
for k=2:K
    p3_matrix=p3_matrix(:,rep_num(k-1)+1:length(p3_matrix(1,:)));
    num=length(q_ma);
    p_ma=q_ma;
    q_ma={};
    qn=0;
    for i=1:num
        for j=1:rep_num(k)
                flag=1;
                 dist_tem=[p_ma{i},p3_matrix(:,j)];
                 dist_temsrt=sort(dist_tem,2);
                for m=1:3
                   flag=flag*ismember(dist_temsrt(m,:),nchoosek(dist3_msrt(m,:),length(dist_temsrt(m,:))),'rows');
                end
                if flag==1
                   qn=qn+1; 
                   q_ma{qn}=[p_ma{i},p3_matrix(:,j)];
                end  
        end
    end
end
p_ma=q_ma;
pn=qn;
qn=0;
q_ma={};
for i=1:pn
    pmatrix=p_ma{i};
    pmatrix=sort(pmatrix,2);
    if pmatrix==dist3_msrt
        qn=qn+1;
        q_ma{qn}=p_ma{i};
    end
end
%if the mathcing is not feasible, then it is regarded as an error
if qn==0
    iteration=iteration+1;
    error_num1=error_num1+K;
    error_num2=error_num3+K;
    error_num3=error_num3+K;
    %error_num_da=error_num_da+1;
    [error_num1,error_num2,error_num3,K,iteration]
    continue
end
l=1;
sumdiffall=zeros(1,qn);
qma={};
coor_tgall={};

while   l<=qn
q_right=q_ma{l};         
i=4;

% for each BS m apply Hungarian Method
while i<=M
    i_tg=zeros(K,2);
% get the initial lcoation  and corresponding distance     
 for j=1:K
    [diff,corest_tg(j,:)]=coorest(q_right(:,j));
 end
 m_cost=zeros(K,K);
 % calculate the cost of match solution
for k=1:K
    for r=1:K
         m_cost(k,r)=abs(dist_m(i,k)-norm(corest_tg(r,:)-coordinate_bs(i,:)));
    end
end
% get the optimal assignment solution using Hungarian Method
[rowindex,cost] = munkres(m_cost);
 q_i=zeros(1,K);
 %allocate the distance to each target
for k=1:K
  q_i(rowindex(k))=dist_m(i,k);
end
q_right=[q_right;q_i];
i=i+1;
end

coor_tg=zeros(K,2);
sumdiffsum=0;
% get the final target location
for j=1:K
   [diff1,coor_tg(j,:)]=coorest(q_right(:,j));
   sumdiffsum=sumdiffsum+sum(diff1);
end
% save target and the correponding value of objetive function
sumdiffall(l)= sumdiffsum;
coor_tgall{l}=coor_tg;
qma{l}=q_right;
l=l+1;
end
%find the optimal matching
l_ind=find(sumdiffall==min(sumdiffall));
dist_opt=qma{l_ind(1)};
final_tg=coor_tgall{l_ind};
[diff_est,rowindex]=diffdist(final_tg,coordinate_tg);
diff_dist=abs(dist_opt-dist_m(:,rowindex));

% calculate the error probability
% error_ind=0;
% for k=1:K
%    error_column=length(find(diff_dist(:,k)>0));
%    if error_column>=M-2
%        error_ind=1;
%    end
% end
% if error_ind==1
%     error_num_da=error_num_da+1;
% end
 error_num1=error_num1+length(find(diff_est>1.5));
 error_num2=error_num2+length(find(diff_est>2));
 error_num3=error_num3+length(find(diff_est>2.5));
[error_num1,error_num2,error_num3,K,iteration]
   iteration=iteration+1;
 end
error_pro(1,K)=error_num1/(K*iteration_max);
error_pro(2,K)=error_num2/(K*iteration_max);
error_pro(3,K)=error_num3/(K*iteration_max);
% error_pro_da(K)=error_num_da/iteration_max;
% error_dist_da(K)=error_dist/iteration_max;
 end

