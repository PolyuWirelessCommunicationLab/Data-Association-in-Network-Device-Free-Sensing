function [diff,coordinate_estimation]=coorest(dist)
global coordinate_bs sigma ini_cor 
l_d=length(dist);
coorini=initial_est(dist(1:3)); 
 x_est=coorini(1);
 y_est=coorini(2);
beta_current=[x_est;y_est];
beta_next=zeros(2,1);
J=zeros(l_d,2);
sdiff=10;
itera=1;
sdiff_next=0;
sdiff_current=0;
while sdiff>0.1&&itera<50
    x_est=beta_current(1);
    y_est=beta_current(2);
    r=zeros(l_d,1);
for i=1:l_d
    J(i,1)=(x_est-coordinate_bs(i,1))/sqrt((x_est-coordinate_bs(i,1))^2+(y_est-coordinate_bs(i,2))^2);
    J(i,2)=(y_est-coordinate_bs(i,2))/sqrt((x_est-coordinate_bs(i,1))^2+(y_est-coordinate_bs(i,2))^2);
    r(i,1)=sqrt((x_est-coordinate_bs(i,1))^2+(y_est-coordinate_bs(i,2))^2)-dist(i);
end 
   lamada=0.2;
    beta_next=beta_current-lamada*pinv(J.'*J)*J.'*r;
    sdiff=norm(beta_next-beta_current);
%     if abs(sdiff_next-sdiff_current)<0.01
%         sdiff=1;
%      end
%       sdiff_next-sdiff_current
     itera=itera+1;
    beta_current=beta_next;
end
coordinate_estimation=beta_current.';
% coordinate_estimation=initial_est(dist(1:3)).';
diff=[];
for n=1:l_d
    diff(n)=abs(norm(coordinate_estimation-coordinate_bs(n,:))-dist(n));
end

end
