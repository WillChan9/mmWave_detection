function [ index, XT ] = cfar_os( xc, N, k, pro_N, Pfa)
alpha=N.*(Pfa.^(-1./N)-1);

index=1+N/2+pro_N/2:length(xc)-N/2-pro_N/2;
XT=zeros(1,length(xc));

for i=index
    cell_left=xc(1,i-N/2-pro_N/2:i-pro_N/2-1);
    cell_right=xc(1,i+pro_N/2+1:i+N/2+pro_N/2);
    cell_all=cat(2,cell_left,cell_right);
    cell_sort=sort(cell_all);

    Z=cell_sort(1,k);

    XT(1,i)=Z.*alpha;
end
XT(1:N/2+pro_N/2)=XT(1+N/2+pro_N/2);
XT(length(xc)-N/2-pro_N/2+1:length(xc))=XT(length(xc)-N/2-pro_N/2);
end

