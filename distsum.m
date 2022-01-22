function dist_sum=distsum(cortg,L,dist)
global coordinate_bs
dist_sum=0;
 for l=1:L
     dist_sum=dist_sum+abs(norm(cortg-coordinate_bs(l,:))-dist(l));
 end
end