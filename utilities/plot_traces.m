%% compare whiten data

%%
%dX=1:24;
%dX=[(1:24) 60:120];
dX = 1:spl/2;
dtah = reshape(dtac,4,spl,[]);
dtah(:,dX,:)=[];
dtah = reshape(dtah,4,[]);
%%
nn = size(dtah,2);
C = dtah*dtah'/nn;
[V,D] = eig(C);
W = V/sqrt(D)*V';
[sum(W,2) W]
%dtaw = W*dtac;
%%
dtaw = readmda('/media/psf/PH/pre.mda');
dX = 1:spl/2;
dtah = reshape(dtaw,4,spl,[]);
dtah(:,dX,:)=[];
dtah = reshape(dtah,4,[]);
dtas = [dtaw;[dtah,dtah]];
%%
sp = 24;
spl = sl + sp;
events=12; morevent = 200;
X=(1:events*spl) + morevent*spl; x_lim=[X(1)-1,X(end)];
Nplot=size(dtas,1);

figure(1);clf
for i = 1:Nplot
    n = mod(2*(i-1),Nplot-1)+1;
    if i == Nplot; n = Nplot; end
    subplot(4,2,n)    
    hold on
    plot(X,dtas(i,X))
    if i>4
        plot(x_lim,[2.5 2.5],'m')
        ylim([-3 6])
    else
        plot(x_lim,[2.5 2.5],'m')
        ylim([-3 6])
    end
    xlim(x_lim)
    y_lim = get(gca,'Ylim');y1=y_lim(1);y2=y_lim(2); 
    for j=0:(events-1)
        x1 = x_lim(1)+j*spl+0.5; 
        x2 = x1 + spl; x3 = x1 + 15;
        if mod(j,2)==0
            fill([x1,x2,x2,x1],[y1,y1,y2,y2],'k','linestyle','none');
        end
%         fill([x1,x3,x3,x1],[y1,y1,y2,y2],'r','linestyle','none');
%         fill([x1,x3,x3,x1],[y1,y1,y2,y2],'r','linestyle','none');
        alpha(0.05);
    end
    ylim(y_lim)
end

%%

count_edges = 1:96:(size(T)+1);
t20 = firing(2,firing(3,:)==20);
t14 = firing(2,firing(3,:)==14);
[counts14,~] = histcounts(t14,count_edges);
