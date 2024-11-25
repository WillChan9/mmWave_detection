function PlotClusterinResult(X,cluster,IDX,xx)
    k=max(IDX);
    Colors=hsv(k);
    Legends = {''};
    for i=0:k
        Xi=X(IDX==i,:);
        if i~=0
            Style = 'x';
            MarkerSize = 5;
            Color = Colors(i,:);
            Legends{end+1} = ['Cluster #' num2str(i)];
        else
            Style = 'o';
            MarkerSize = 5;
            Color = [0 0 0];
            if ~isempty(Xi)
                Legends{end+1} = 'Noise';
            end
        end
        if ~isempty(Xi)
            plot(Xi(:,1),Xi(:,2),Style,'MarkerSize',MarkerSize,'Color',Color);
        end
        hold on;
    end
    for j=1:k
        x=cluster(j,6);% bug when only one cluster
        y=cluster(j,8);
        w=cluster(j,5)-cluster(j,6);
        h=cluster(j,7)-cluster(j,8);
        rectangle('Position',[x,y,w,h],'Curvature',[0,0],'LineWidth',2,'LineStyle','-','EdgeColor',Colors(j,:));
        text(x-0.4,y-0.4, num2str(j),'FontSize',12,'FontWeight','Bold','Color',Colors(j,:));
    end
    hold off;
    set(gca,'XLim',[-20 20]);
    set(gca,'YLim',[0 40]);
    %axis equal;
    grid on;
    legend(Legends);
    legend('Location', 'NorthEastOutside');
    title(['¾ÛÀà½á¹û frame No:' num2str(xx)])
end

