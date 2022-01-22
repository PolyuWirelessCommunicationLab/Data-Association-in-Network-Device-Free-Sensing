function [diff,rowindex]=diffdist(A,B)
global K R
diff=zeros(1,K);
diffdistant=zeros(K,K);
for i=1:K
    for j=1:K
        ab=A(i,:)-B(j,:);
        diffdistant(i,j)=norm(ab);
        j=j+1;
    end
 i=i+1;
end
[rowindex,cost]=munkres(diffdistant);
cor_index=rowindex;
for i=1:K
  diff(i)=diffdistant(i,rowindex(i));  
end
end