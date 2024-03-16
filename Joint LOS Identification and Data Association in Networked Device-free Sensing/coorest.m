function [diff,coordinate_estimation]=coorest(dist)
global coordinate_bs m_nonzero M
m_nonzero=find(dist>0);
dist=dist(dist>0);
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
while sdiff>0.001&&itera<30
    x_est=beta_current(1);
    y_est=beta_current(2);
    r=zeros(l_d,1);
for i=1:length(m_nonzero)
    J(i,1)=(x_est-coordinate_bs(m_nonzero(i),1))/sqrt((x_est-coordinate_bs(m_nonzero(i),1))^2+(y_est-coordinate_bs(m_nonzero(i),2))^2);
    J(i,2)=(y_est-coordinate_bs(m_nonzero(i),2))/sqrt((x_est-coordinate_bs(m_nonzero(i),1))^2+(y_est-coordinate_bs(m_nonzero(i),2))^2);
    r(i,1)=sqrt((x_est-coordinate_bs(m_nonzero(i),1))^2+(y_est-coordinate_bs(m_nonzero(i),2))^2)-dist(i);
end 
% for i=1:length(m_nonzero)-1
%     J(i+length(m_nonzero),1)=0.5*(x_est-coordinate_bs(m_nonzero(i),1))/sqrt((x_est-coordinate_bs(m_nonzero(i),1))^2+(y_est-coordinate_bs(m_nonzero(i),2))^2)+0.5*(x_est-coordinate_bs(m_nonzero(i+1),1))/sqrt((x_est-coordinate_bs(m_nonzero(i+1),1))^2+(y_est-coordinate_bs(m_nonzero(i+1),2))^2);
%     J(i+length(m_nonzero),2)=0.5*(y_est-coordinate_bs(m_nonzero(i),2))/sqrt((x_est-coordinate_bs(m_nonzero(i),1))^2+(y_est-coordinate_bs(m_nonzero(i),2))^2)+0.5*(y_est-coordinate_bs(m_nonzero(i+1),2))/sqrt((x_est-coordinate_bs(m_nonzero(i+1),1))^2+(y_est-coordinate_bs(m_nonzero(i+1),2))^2);
%     r(i+length(m_nonzero),1)=0.5*sqrt((x_est-coordinate_bs(m_nonzero(i),1))^2+(y_est-coordinate_bs(m_nonzero(i),2))^2)+0.5*sqrt((x_est-coordinate_bs(m_nonzero(i+1),1))^2+(y_est-coordinate_bs(m_nonzero(i+1),2))^2)-0.5*dist(i+length(m_nonzero));
% end
   lamada=0.2;
    beta_next=beta_current-lamada*pinv(J)*r;
    sdiff=norm(beta_next-beta_current);
    beta_current=beta_next;
    itera=itera+1;
end
coordinate_estimation=beta_current.';
% coordinate_estimation=initial_est(dist(1:3)).';
diff=[];
for n=1:length(m_nonzero)
    diff(n)=(norm(coordinate_estimation-coordinate_bs(m_nonzero(n),:))-dist(n))^2;
end
end
