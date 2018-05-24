clear
cd('/Users/hupeiyao/Movies/PH')
addpath(genpath(pwd))

%% !! Read Tetrode Data !!

doTetrodeData=1; doPosition=1; doTetrodeTS=1;
drIn = 'SpikeSort/'; 
%bpf_name = '0615R16BC-cl.bpf';
bpf_name = 'M16_2.bpf';

bpf2matVolts([drIn bpf_name],'',0,doTetrodeData,doTetrodeTS,doPosition,0,0);  
%%
bpf_name = 'M16_2.bpf';
load(['DATA/' bpf_name(1:end-4) '.mat']);
%%
iCh = 1;

nchs=4;
spl = tetrodeSamplesPerBlock;
Chs = unique(tetrodeChannel)' 
Data = reshape(tetrodeData,nchs,spl,[]);
idxs = tetrodeChannel==Chs(iCh);
T = tetrodeTimestamp(idxs);
U = tetrodeUnit(idxs);
min_diffT = min(diff(T)) 
unique(U')
dtac = reshape(Data(:,:,idxs),nchs,[]);


%% -- Plot Peak Distribution --

F = readmda('/Users/hupeiyao/Movies/PH/firing.mda');
lmap = readmda('/Users/hupeiyao/Movies/PH/lmap.mda');
% fets = readmda('/Users/hupeiyao/Movies/PH/fets.mda');

% resc = 2; lmap(resc,1)=resc;
% lmap(:,1)=lmap(:,2);
% disc=6; lmap(disc,1)=0;
% disc=7; lmap(disc,1)=0;
% disc=3; lmap(disc,1)=0;

sp=24;
spl=sl+sp;
figure;
edges = min(mod(F(2,:),spl)):max(mod(F(2,:),spl))+1;
%edges = 7:80;

Nplot=length(lmap);
for iu = 1:Nplot
    subplot(4,ceil(Nplot/4),iu); hold on
    idx = F(3,:)==iu;
    Tu = F(2,idx);
    Tr = mod(Tu,spl);
        if lmap(iu,1)==0
            histogram(Tr,edges,'FaceColor','r','FaceAlpha',0.25);title(['reject' num2str(iu)])
        elseif lmap(iu,1)~=lmap(iu,2)
            histogram(Tr,edges,'FaceColor','b','FaceAlpha',0.25);title(['merge' num2str(lmap(iu,1))])
        else
            histogram(Tr,edges); title(['keep' num2str(iu)])
        end
    y_lim = get(gca,'Ylim');
    plot([39 39], y_lim, 'm', 'LineWidth',1);
    ylim(y_lim)
    xlim([31 edges(end)])
end

%% -- Plot Place Fields --
%resc = [2 4 7];
%lmap(resc,1)=lmap(resc,2);

figure;
pl1 = 4; pl2 = 7;
masking = 1;

bS = 6;
sigma = 1;
colormap hot;
%clim_scale=0.9;

Sorter = {'M','W'};
Color  = {'g','w'};
for Wc = 0:1 % Wc = 0 to do MS
if Wc
    keepCluster = U;
else
    keepCluster = zeros(1,size(F,2));
    for i = 1:size(lmap,1)
        keepCluster(F(3,:)==lmap(i,2)) = lmap(i,1);
    end
end

units = unique(keepCluster);
if units(1)==0; units(1)=[]; end % exclude noise cluster

Ax = ceil(double(max(roomXY(:,1)))/bS);
Ay = ceil(double(max(roomXY(:,2)))/bS);
A = zeros(Ax,Ay,length(units)+1);

times = 3:length(roomTimeStamps);
%roomXY(times(1),1)*roomXY(times(1),2)

counts = zeros(length(units),size(times,2)-1);

for u = 1:length(units) 
    if Wc
        Tu = T(keepCluster==units(u));
    else
        f = F(2,keepCluster==units(u)); 
        Tr = round((mod(f,spl)-10-sp)/4.8);
        Tu = T(ceil(f/spl))+cast(Tr','uint32');
    end
    [counts(u,:),~] = histcounts(Tu,roomTimeStamps(times));
end

[~,sid] = sort(sum(counts,2),'descend');
counts = counts(sid,:);

for i = times(1:end-1)
    x = double(roomXY(i,1)); y = double(roomXY(i,2)); 
    if x+y == 0 
        continue 
    end
    A(ceil(x/bS),ceil(y/bS),1) = A(ceil(x/bS),ceil(y/bS),1) + roomTimeStamps(i)-roomTimeStamps(i-1);
    for u = 1:length(units)
        A(ceil(x/bS),ceil(y/bS),u+1) = A(ceil(x/bS),ceil(y/bS),u+1) + counts(u,i-times(1)+1);
    end
end


% -- plot place fields
rowNeed = ceil(sum(lmap(:,1)>0)/pl2);
DT = A(:,:,1); %imgaussfilt(A(:,:,1),sigma);
subplot(pl1,pl2,(Wc+1)*rowNeed*pl2);

imagesc(DT);ax=gca;ax.CLim = 0.8*ax.CLim;axis off;%axis equal;
title([num2str(round(roomTimeStamps(end-1)/10000)) 's'],'color','y');
imAlpha = ~A(:,:,1)==0;

%

kernel = [1 1 1; 1 0 1; 1 1 1]/8;

for nplot = 1:size(A,3)-1
    subplot(pl1,pl2,nplot+Wc*rowNeed*pl2)
    PF = A(:,:,nplot+1)./A(:,:,1);
    Pa = ~isnan(PF);
    PF(isnan(PF))=0; 
    F8 = conv2(PF,kernel,'same');
    Fn = conv2(double(Pa),kernel,'same');
    F_8 = F8./Fn; 
    coh = corr(PF(Pa),F_8(Pa));
    
%    imagesc(PF,'AlphaData',Pa); 
    Pn = imgaussfilt(double(Pa),sigma); 
    PF = imgaussfilt(PF,sigma);
    if masking
        imagesc(PF./Pn,'AlphaData',Pa); 
    else
        imagesc(PF./Pn);
    end
    
    %colorbar; 
    %axis equal
    axis off
    ax=gca;
    C_Lim = ax.CLim;
    Title = [sprintf('%.2f',coh) '<' num2str(units(sid((nplot)))) '>' sprintf('%.1f',C_Lim(2)*10000)];
    title([Title '-' num2str(sum(counts(nplot,:)))],'fontsize',12,'color',Color{Wc+1});
    ax.CLim = clim_scale*C_Lim;
end
set(gcf,'Color',[0 0 0])

end

%% MountainSort v.s. Wclust

kCluster = zeros(1,size(F,2));
for j = 1:size(lmap,1)
    kCluster(F(3,:)==lmap(j,2)) = lmap(j,1);
end

count_edges = 1:spl:(size(dtac,2)+1);
[ccs,~] = histcounts(F(2,kCluster~=0),count_edges);

Cms0=sum(kCluster==0); Cmn0=sum(U==0);
Cmn=[];Cms=[];
for i = 1:max(U)
    Cmn = [Cmn; sum(U==i)];
end
for i = unique(kCluster)
   if i==0; continue; end
   Cms =  [Cms; sum(kCluster==i)];
end

[Cms0,Cmn0;
[zeros(length(Cmn)-length(Cms),1);sort(Cms)],[zeros(length(Cms)-length(Cmn),1);sort(Cmn)];
 sum(Cms),sum(Cmn)]

%[mean(ccs(U==0)) mean(ccs(U~=0));mean(U(ccs==0)~=0) mean(U(ccs~=0)~=0)]
