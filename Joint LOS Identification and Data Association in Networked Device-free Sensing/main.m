global coordinate_bs dist_m M K R dist_bs sigma_d dist_sum_m D_sum  dist_m_0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%parameter zone
M=4;  % the number of BSs
R=80; % the size of region of intrest
P_nl=0.5;  % the probability that there is one NLOS path for each target
P_b=0.01;  % the probability that the LOS path is blocked
N=3300;     % the number of subcarriers   
power=20;   % the transmit power
B=4*10^8;   % the bandwidth
delta_d=3*10^8/B;  % the range resolution
L0=200;     % the maximum number of resolvable path  
Rd=75;  
L=220;     % L+tau_max
N0=10^(-16.9)*10^(-3);  % the noise power spectrum
snr=power/(N*B*N0);     %the SNR
W=zeros(N,N);       %DFT matrix
for i=1:N
   for j=1:N
       W(i,j)=exp(-1i*2*pi*(i-1)*(j-1)/N);
   end
end
G_t=W(:,1:M*L);   
inv_G_t=pinv(G_t);
error_pro=zeros(6,6);
sigma_d=delta_d/4;
delta_0=4*sigma_d;
dist_m_0={};
num_of_G1=zeros(2,7);
num_of_G2=zeros(2,7);
num_of_G3=zeros(2,7);
for K=2:7% the number of target
    tic
    iteration_max=10000; % the maximum number of iterations
    iteration=1; 
    error_fa_0=0;
    error_fa_1=0;
    error_fa_2=0;
    error_fa_3=0;
    error_md_0=0;
    error_md_1=0;
    error_md_2=0;
    error_md_3=0;
    num_error=0;
    K_fa=0;
    tau_error=0;
    num_error_tau=0;
    while iteration<=iteration_max
        coordinate_bs=R*(rand(M,2)-0.5);  % randomly generate the coordinats of BSs
        coordinate_tg=R*(rand(K,2)-0.5);  % randomly generate the coordinats of TGs
        dist_m=zeros(M,K);                % store the estimated direct ranges between BSs and targets
        dist_sum_m=zeros(M,M,K);         % store the estimated sum ranges between BSs and targets
        dist=zeros(M,K);                % store the true distances between BSs and targets
        same_range=0;
        % the true ranges between any two BSs
        dist_bs=zeros(M,M);
        for m1=1:M
            for m2=1:M
                dist_bs(m1,m2)=norm(coordinate_bs(m1,:)-coordinate_bs(m2,:));
            end
        end
        % the true ranges between targets and BSs
        dist_r=zeros(M,M);
         for i=1:M
            for j=1:K
                dist_r(i,j)=norm(coordinate_bs(i,:)-coordinate_tg(j,:));
            end
         end
        % judge whether there are more one target in the same range
        for m=1:M
             for u=1:M
                dist_m_u=[];
                    for k=1:K
                        dist_m_u=[dist_m_u,ceil((dist_r(m,k)+dist_r(u,k))/delta_d)];
                    end
                    if u~=m
                        dist_m_u=[dist_m_u,ceil(dist_bs(m,u)/delta_d)];
                    end
                    if length(unique(dist_m_u))<length(dist_m_u)
                        same_range=1;
                        break
                    end
             end
        end
        if same_range==1
           continue
        end
        % Consider detection coverage and LOS blockage 
        for i=1:M
            for j=1:K
               if norm(coordinate_bs(i,:)-coordinate_tg(j,:))<=Rd && rand(1)>P_b
                    dist(i,j)=dist_r(i,j);
               end
            end
        end
         %judge whether each target is detected by at least three BSs
        for k=1:K
            if length(find(dist(:,k)>0))<3
                same_range=1;
                break
            end
        end
        if same_range==1
           continue
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % generate channels
        path_loss=zeros(M,M,K+1);
        H_all=zeros(M,M,K+1);
        for m=1:M
            for u=1:M
                for k=1:K
                     if dist(u,k)>0&&dist(m,k)>0
                         path_loss(m,u,k) = 10^((-143.3- 20*log10(dist(u,k)*1e-3)-20*log10(dist(m,k)*1e-3))/10);
                         theta=2*pi*rand(1);
                         H_all(m,u,k)=exp(1i*theta)*sqrt(path_loss(m,u,k));
                     end
                end
                if u~=m
                     path_loss(m,u,K+1) = 10^((-103.3- 20*log10(norm(coordinate_bs(m,:)-coordinate_bs(u,:))*1e-3))/10);
                     theta=2*pi*rand(1);
                     H_all(m,u,K+1)=exp(1i*theta)*sqrt(path_loss(m,u,K+1));
                end
            end
        end
        H=zeros(M*L,M);
        flag=0;
        rand_nl=rand(K,1);
        ind={};
        % introduce STO
        tau_bs=zeros(M,M);
        for m1=1:M
            for m2=m1+1:M
                tau_bs(m1,m2)=ceil(min(floor(7.5/delta_d),ceil(dist_bs(m1,m2)/delta_d)-1)*rand(1));
                tau_bs(m2,m1)=-tau_bs(m1,m2);
            end
        end
        % obtain the real channels with STO
        for m=1:M
            ind_m=zeros(M,K+1);
            ind_nl_m=zeros(M,K+1);
            for u=1:M
                for k=1:K
                    ind_m(u,k)=ceil((norm(coordinate_bs(m,:)-coordinate_tg(k,:))+norm(coordinate_bs(u,:)-coordinate_tg(k,:)))/delta_d);
                    if abs(H_all(m,u,k))>0
                        H(ind_m(u,k)+tau_bs(m,u)+(u-1)*L,m)= H(ind_m(u,k)+tau_bs(m,u)+(u-1)*L,m)+H_all(m,u,k);
                    end
                end
                if u~=m
                    ind_m(u,K+1)=ceil(norm(coordinate_bs(m,:)-coordinate_bs(u,:))/delta_d);
                    H(ind_m(u,K+1)+tau_bs(m,u)+(u-1)*L,m)= H(ind_m(u,K+1)+tau_bs(m,u)+(u-1)*L,m)+H_all(m,u,K+1);
                end
                for n=1:K
                    if rand_nl(n)<P_nl && ind_m(u,n)<L-floor(7.5/delta_d)
                        ind_umn=[ind_m(u,n)+1:L-floor(7.5/delta_d)];
                        nz_m=ind_umn(randperm(numel(ind_umn),1));
                        H(nz_m+tau_bs(m,u)+(u-1)*L,m)=H(nz_m+tau_bs(m,u)+(u-1)*L,m)+(2.5+10*rand(1))*10^(-6)*exp(1i*2*pi*rand(1));
                        ind_nl_m(u,n)=nz_m;
                    end
                end
            end
            ind{m}=ind_m;
            ind_nl{m}=ind_nl_m;
        end
        flag_error=0;
        tau_bs_est=zeros(M,M);
        threshold=0.065;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Phase I
        % estimate the channels
        for m=1:M
            x_m=sqrt(snr)*H(:,m);
            x_ind=find(abs(x_m)>0);
            y_m=G_t*x_m+sqrt(0.5)*(randn(N,1)+1j*randn(N,1));
            x_est=inv_G_t*y_m;
            for u=m:M
                if u==m
                    ind_est_mm=find(abs(x_est(1+(u-1)*L:L+(u-1)*L))> threshold);
                    ind_est_mm=reshape(ind_est_mm,[1,length(ind_est_mm)]);
                    dist_m(m,[1:length(ind_est_mm)])=(ind_est_mm*delta_d-delta_d/2)/2;
                    ind_true_mm=find(abs(x_m(1+(u-1)*L:L+(u-1)*L))>0);
                    if length(setdiff(ind_true_mm,ind_est_mm))>0||length(setdiff(ind_est_mm,ind_true_mm))>0
                        flag_error=1;
                    end
                else
                    ind_est_um=find(abs(x_est(1+(u-1)*L:L+(u-1)*L))>threshold);
                    ind_est_um=reshape(ind_est_um,[1,length(ind_est_um)]);
                    tau_bs_est(m,u)=min(ind_est_um)-ceil(dist_bs(u,m)/delta_d);
                    ind_est_um=ind_est_um-tau_bs_est(m,u);
                    ind_true_um=find(abs(x_m(1+(u-1)*L:L+(u-1)*L))>0);
                    ind_true_um= ind_true_um-tau_bs(m,u);
                    if length(setdiff(ind_true_um,ind_est_um))>0||length(setdiff(ind_est_um,ind_true_um))>0
                        flag_error=1;
                    end
                    ind_est_um(find(ind_est_um==min(ind_est_um)))=[];
                    for i=1:length(ind_est_um)
                        dist_sum_m(m,u,i)=(ind_est_um(i)*delta_d-delta_d/2);
                        dist_sum_m(u,m,i)=(ind_est_um(i)*delta_d-delta_d/2);
                    end
                end
            end   
        end
        if flag_error==1
            num_error=num_error+1;
        end
        % obtain the estimated direct ranges for each BS m
        for m=1:M
            dist_m_0{m}=unique(dist_m(m,find(dist_m(m,:)>0)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Phase II
        %for any two BSs m and u to find the possible tuples (d_{m,k},d_{u,k},d_{m,u,k})
        D_sum={};
        ind_2bs=nchoosek([1:M],2);
        ind_2bs_length=length(ind_2bs(:,1));
        for ind_num=1:ind_2bs_length
          m1=ind_2bs(ind_num,1);
          m2=ind_2bs(ind_num,2);
          d1=dist_m(m1,:);
          d2=dist_m(m2,:);
          d12=reshape(dist_sum_m(m1,m2,:),[1,length(dist_sum_m(m1,m2,:))]);
          d1=unique(d1);d1=d1(find(d1>0));
          d2=unique(d2);d2=d2(find(d2>0));
          d12=unique(d12);d12=d12(find(d12>0));
          D12=[];
         for i=1:length(d1)
             for j=1:length(d2)
                for l=1:length(d12)
                    if abs(d12(l)-d1(i)-d2(j))<=delta_0
                        D12=[D12;[d1(i),d2(j)]];
                        break
                    end
                end
            end
         end
          D_sum{m1,m2}=D12;
          D_sum{m2,m1}=D12;
        end
        % iterativly localizing the targets that are exactly detected by M,
        % M-1,....,3 BSs
        D_sum_origin=D_sum;
        dist_opt=[];
        for m=0:M-3
            %find the targets detected by M-m BSs
            num_det_BS=M-m;
            ind_com=nchoosek([1:M],num_det_BS);
            q_m=[];
            for n=1:size(ind_com,1)
                num_g1=1;
                for l=1:length(ind_com(n,:))
                    num_g1=num_g1*length(dist_m_0{ind_com(n,l)});
                end
                % find {G}^l
                if m==0
                    num_of_G1(1,K)=num_of_G1(1,K)+ num_g1;
                else
                    num_of_G1(2,K)=num_of_G1(2,K)+num_g1;
                end
                q_m_n=tar_det(ind_com(n,:));
                q_m_all=zeros(M,size(q_m_n,2));
                q_m_all(ind_com(n,:),:)=q_m_n;
                % find bar{G}^l
                if m==0
                    num_of_G2(1,K)=num_of_G2(1,K)+size(q_m_n,2);
                else
                    num_of_G2(2,K)=num_of_G2(2,K)+size(q_m_n,2);
                end
                q_m_m=[];
                for i=1:size(q_m_n,2)
                    [diff,coor]=coorest(q_m_all(:,i));
                    if sum(diff)<=2*num_det_BS*sigma_d^2
                        q_m=[q_m,q_m_all(:,i)];
                        q_m_m=[q_m_m,q_m_all(:,i)];
                    end
                end
                % find tilde{G}^l
                if m==0
                    num_of_G3(1,K)=num_of_G3(1,K)+size(q_m,2);
                else
                    num_of_G3(2,K)=num_of_G3(2,K)+size(q_m_m,2);
                end
            end
            %select the targets detected by M-m BSs under the constraint (20)
            if size(q_m,2)>1
               if sum(q_m(1,:))==0 
                   q_m=q_m(2:M,:);
                   q_ma=findma(q_m,M-1);
                   for n=1:length(q_ma)
                        q_ma{n}=[zeros(1,size(q_ma{n},2));q_ma{n}];
                   end
               else 
                   q_ma=findma(q_m,M);
               end
            else
               q_ma=findma(q_m,M);
            end
            if length(q_ma)==0
                continue
            end
            % find the optimal one with the minimum localization residual
            if length(q_ma)>1
                cost_all=[];
                for l=1:length(q_ma)
                    dist_m_l=q_ma{l};
                    cost_l=0;
                    for k=1:size(dist_m_l,2)
                        [diff_k,cootg]=coorest(dist_m_l(:,k));
                        cost_l=cost_l+sum(diff_k);
                    end
                    cost_all=[cost_all,cost_l];
                end
                ind_opt=find(cost_all==min(cost_all));
                q_ma_opt=q_ma{ind_opt(1)};
            else
                q_ma_opt=q_ma{1};
            end
            dist_opt=[dist_opt,q_ma_opt];
            %delete those assigned measurements
            ind_j=[];
            for j=1:size(q_ma_opt,2)
                for m=1:M
                    if q_ma_opt(m,j)>0
                        dist_m_0{m}=setdiff(dist_m_0{m},q_ma_opt(m,j));
                    end
                end
                for l=1:length(ind_2bs)
                    ind_2bs_l=ind_2bs(l,:);
                    D_12=D_sum{ind_2bs_l(1),ind_2bs_l(2)};
                    if length(D_12)>0
                        ind_j=union(find(D_12(:,1)==q_ma_opt(ind_2bs_l(1),j)),find(D_12(:,2)==q_ma_opt(ind_2bs_l(2),j)));
                        D_12(ind_j,:)=[];
                        D_sum{ind_2bs_l(1),ind_2bs_l(2)}=D_12;
                        D_sum{ind_2bs_l(2),ind_2bs_l(1)}=D_12;
                    end
                end
            end
        end
        % delete those ineffective targets
        K_present=size(dist_opt,2);
        coor_tg=zeros(K_present,2);
        ind_coortg=[];
        for k=1:K_present
            [diff_k,coor_tg(k,:)]=coorest(dist_opt(:,k));
                if sum(diff_k)>2*M*sigma_d^2||abs(coor_tg(k,1))>R/2+delta_d||abs(coor_tg(k,2))>R/2+delta_d
                ind_coortg=[ind_coortg,k];
                end
        end
        coor_tg(ind_coortg,:)=[];
        K_present=size(coor_tg,1);
        %calculate MD and FA ratio
        diff_mat=[];
        for k1=1:K
            for k2=1:K_present
              diff_mat(k1,k2)=norm(coordinate_tg(k1,:)-coor_tg(k2,:));  
            end
        end
        correct_num_0=0;
        correct_num_1=0;
        correct_num_2=0;
        correct_num_3=0;
        diff_mat_origin=diff_mat;
        ind_1=[];
        ind_2=[];
            while min(min(diff_mat))<=1
                [k1,k2]=find(diff_mat==min(min(diff_mat)));
                diff_mat_k1k2=diff_mat(k1,k2);
                if diff_mat_k1k2<0.375
                    diff_mat(k1,:)=10;
                    diff_mat(:,k2)=10;
                    correct_num_0=correct_num_0+1;
                end
                if diff_mat_k1k2<0.5
                    diff_mat(k1,:)=10;
                    diff_mat(:,k2)=10;
                    correct_num_1=correct_num_1+1;
                end
                if diff_mat_k1k2<0.75
                    diff_mat(k1,:)=10;
                    diff_mat(:,k2)=10;
                    correct_num_2=correct_num_2+1;
                end
                if diff_mat_k1k2<1
                    diff_mat(k1,:)=10;
                    diff_mat(:,k2)=10;
                    correct_num_3=correct_num_3+1;
                end
            end
        md_num_0=K-correct_num_0;md_num_1=K-correct_num_1;md_num_2=K-correct_num_2;md_num_3=K-correct_num_3;
        fa_num_0=K_present-correct_num_0;fa_num_1=K_present-correct_num_1;fa_num_2=K_present-correct_num_2;fa_num_3=K_present-correct_num_3;
         error_fa_0=error_fa_0+fa_num_0;
         error_fa_1=error_fa_1+fa_num_1;
         error_fa_2=error_fa_2+fa_num_2;
         error_fa_3=error_fa_3+fa_num_3;
         error_md_0=error_md_0+md_num_0;
         error_md_1=error_md_1+md_num_1;
         error_md_2=error_md_2+md_num_2;
         error_md_3=error_md_3+md_num_3;
        K_fa=K_fa+K_present;
        [error_md_1,error_fa_1,K,iteration]
        iteration=iteration+1;
        tau_bs_est=tau_bs_est-tau_bs_est';
        tau_error=tau_error+sum(sum(abs(tau_bs-tau_bs_est)));
        if sum(sum(abs(tau_bs-tau_bs_est)))~=0
           num_error_tau= num_error_tau+1;
        end
    end
    error_pro(1,K)=error_md_0/(K*iteration_max);
    error_pro(2,K)=error_md_1/(K*iteration_max);
    error_pro(3,K)=error_md_2/(K*iteration_max);
    error_pro(4,K)=error_md_3/(K*iteration_max);
    error_pro(5,K)=error_fa_0/(K*iteration_max);
    error_pro(6,K)=error_fa_1/(K*iteration_max);
    error_pro(7,K)=error_fa_2/(K*iteration_max);
    error_pro(8,K)=error_fa_3/(K*iteration_max);
    error_pro(9,K)=num_error/(iteration_max); 
    error_pro(10,K)=num_error_tau/(iteration_max); 
    error_pro(11,K)=K_fa/(iteration_max); 
    toc
end
    num_of_G1=num_of_G1/iteration_max;
    num_of_G2=num_of_G2/iteration_max;
    num_of_G3=num_of_G3/iteration_max;
