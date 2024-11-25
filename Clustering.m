function IDX=Clustering(X,epsilon,MinPts)
%input X£¬column value :x,y,amp,v
    %cluster number
    C=0;
    %data size
    n=size(X,1);
    %amp flag£¬if amp>mean then not take for noise
    amp_mean=mean(X(:,3));
    amp_flag=X(:,3)>(zeros(size(X,1),1)+amp_mean);
    %data clusters ID
    IDX=zeros(n,1);
    % n*n distance matrix
    D=pdist2(X(:,1:2),X(:,1:2));
    %mark all the plot as not visited
    visited=false(n,1);
    isnoise=false(n,1);
    %go through all the plots
    for i=1:n
        if ~visited(i)
            %if plots is not visited then 
            visited(i)=true; 
            % get epsilon neigbourhood index
            Neighbors=RegionQuery(i);
            % compare neighbourhood point quantity with MinPts
            if (numel(Neighbors)<MinPts && amp_flag(i)==0) 
                isnoise(i)=true;               
            elseif (X(i,1)<1 && X(i,1)>-1 && X(i,2)<3)%judge if it's ground noise
                isnoise(i)=true;  
            else
                %create cluster C using X(i,:) as core
                C=C+1;
                %find other plots in neighbourhood 
                vi=X(i,4);
                ai=X(i,3);
                ExpandCluster(i,Neighbors,C,vi,ai);
            end
            
        end
    end
    
    function ExpandCluster(i,Neighbors,C,vi,ai)
        %mark i-th plot in C
        IDX(i)=C;
        k = 1;
        while true
            j = Neighbors(k);
            if ~visited(j)
                visited(j)=true;
                Neighbors2=RegionQuery(j);
                if (numel(Neighbors2)>=MinPts && abs(X(j,3)-ai)<1900 && abs(X(j,4)-vi)<2.7)
                    %expand neighbourhood
                    Neighbors=[Neighbors Neighbors2];   
                end
            end
            %if not mark as cluster C
            if IDX(j)==0
                IDX(j)=C;
            end
            k = k + 1;
            if k > numel(Neighbors) %go through all the points in Neighbor
                break;
            end
        end
    end
    %get neighboring plots
    function Neighbors=RegionQuery(i)
        Neighbors=find(D(i,:)<=epsilon);
    end
end
