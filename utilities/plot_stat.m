%%
clear
noth = [3 5];


%%
% %% isoi distribution

for nIsoI = 1:length(noth) 
IsoI = IsoIs{nIsoI};

%IsoI = IsoI(IsoI(:,5)>5,:); % MS clusters
%IsoI = IsoI(IsoI(:,4)>2,:);

IsoI = IsoI(IsoI(:,1)>0,:);
IsoI = IsoI(IsoI(:,2)>0,:);
IsoIp{nIsoI} = IsoI(:,2);
end

figure;hold on
edges =0.25:0.5:14.25;

hist(IsoIp{1},edges); 
hist(IsoIp{2},edges); 
h = findobj(gca,'Type','patch');
h(1).FaceColor = 'y'; h(1).FaceAlpha = 0.4;
h(2).FaceColor = 'g'; h(2).FaceAlpha = 0.7;
xlim([edges(1) edges(end)]);
xticks(ceil(edges(1)):floor(edges(end)));
xlabel('IsoI')
legend({'noise\_overlap=0.03','noise\_overlap=0.05'})

%%
i14 = F(3,:)==14;
i20 = F(3,:)==20;

t14 = F(2,i14);
t20 = F(2,i20);
%%
figure(4);clf;hold on
ax = gca;
outerpos = ax.OuterPosition; ti = ax.TightInset; 
left = outerpos(1) + ti(1)*1.5;
bottom = outerpos(2) + ti(2)*0.5;
ax_width = outerpos(3) - ti(1)*1.5 - ti(3)*2.5;
ax_height = outerpos(4) - ti(2)*0.5 - ti(4);
ax.Position = [left bottom ax_width ax_height];

plot(t20(t20<200000),1,'b+')
plot(t14(t14<200000),0.8,'r+')

ylim([0,2])

%%
edges_end = ceil((max(F(2,:))/96))*96;
count_edges = 1:96:edges_end;

T(1,:) = F(2,F(3,:)==14);
T(2,:) = F(2,F(3,:)==20);

[counts14,~] = histcounts(t14,count_edges);





