clear
cd('/Users/hupeiyao/Movies/PH')
addpath(genpath(pwd))

%% !! Read Data to MAT !!
% drIn = '/Users/hupeiyao/Movies/PH/DATA/'; 
% files = dir('*.bpf');
% for i = 1:length(files)
%     bpf_name = files(i).name
%     bpf2matPF([drIn bpf_name]);
% end

%% load data

bpf_name = 'M16_2.bpf';
load(['DATA/' bpf_name(1:end-4) '.mat']);

nchs=4;
sl = tetrodeSamplesPerBlock;
Chs = unique(tetrodeChannel)' 
Data = reshape(tetrodeData,nchs,sl,[]);

%%
iCh = 1;
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
% lmap(:,1) = lmap(:,2);
% disc=6; lmap(disc,1)=0;
% disc=7; lmap(disc,1)=0;
% disc=3; lmap(disc,1)=0;

sp=24;
spl=sl+sp;
figure;
edges = min(mod(F(2,:),spl)):max(mod(F(2,:),spl))+1;

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
    plot([40 40], y_lim, 'm', 'LineWidth',1);
    ylim(y_lim)
    xlim([31 edges(end)])
end

%% -- Plot Place Fields --
resc = 4;

lmap(:,1) = lmap(:,2);

figure;
pl1 = 4; pl2 = 7;
sorting = 0;
masking = 1;

bS = 8;
sigma = 1; % for smooting place field

colormap hot;
Sorter = {'M','W'};
Color  = {'g','w'};

for Wc = 0:0 % Wc = 0 to do MS
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
A = zeros(Ax, Ay, length(units) + 1);

times = 3:length(roomTimeStamps);
%roomXY(times(1),1)*roomXY(times(1),2)

counts = zeros(length(units),size(times,2)-1);

for u = 1:length(units) 
    if Wc
        Tu = T(keepCluster==units(u));
    else
        f = F(2,keepCluster==units(u)); 
        Tr = round((mod(f,spl)-10-sp)/4.8);
        Tu = T(ceil(f/spl)) + cast(Tr','uint32');
    end
    [counts(u,:),~] = histcounts(Tu,roomTimeStamps(times));
end

if sorting
    [~,sid] = sort(sum(counts,2),'descend');
    counts = counts(sid,:);
else
    sid = 1:max(units);
end

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

% -- Plot Place Fields --
rowNeed = ceil(sum(lmap(:,1)>0)/pl2);
PT = A(:,:,1); %imgaussfilt(A(:,:,1),sigma);
subplot(pl1,pl2, 2*(Wc+1)*pl2);
imagesc(PT);ax=gca;ax.CLim = 0.8*ax.CLim;axis off;%axis equal;
title([num2str(round(roomTimeStamps(end-1)/10000)) 's'],'color','y');


P_valid = ~(PT==0);
kernel = [1 1 1; 1 0 1; 1 1 1]/8;

for nplot = 1:size(A,3)-1
    subplot(pl1,pl2,nplot)%+Wc*rowNeed*pl2)
    PS = A(:,:,nplot+1);
    PF = PS./PT;
    
    % surround
    PS_s = conv2(PS,kernel,'same'); 
    PT_s = conv2(PT,kernel,'same');
    PF_s = PS_s./PT_s;
    to_use = and(PT>2000,~(PT_s==0));
    coh = corr(PF(to_use),PF_s(to_use));
    
    % filtered
    PS_f = imgaussfilt(PS,sigma);
    PT_f = imgaussfilt(PT,sigma);
    PF_f = PS_f./PT_f;   
    %P_n = imgaussfilt(P_valid,sigma);

    Pnan = ~isnan(PF);
    PF(isnan(PF))=0; 
    F8 = conv2(PF,kernel,'same');
    Fn = conv2(double(Pnan),kernel,'same');
    F_8 = F8./Fn; 
    
    if masking
        imagesc(PF_f,'AlphaData',P_valid); 
    else
        imagesc(PF_f);
    end
    
    %colorbar;
    %axis equal
    axis off
    ax=gca;
    C_Lim = ax.CLim;
    Title = [sprintf('%.2f',coh) '<' num2str(units(sid((nplot)))) '>' sprintf('%.1f',C_Lim(2)*10000)];
    title([Title '-' num2str(sum(counts(nplot,:)))],'fontsize',13,'color',Color{Wc+1});
    if coh>0.49
        ax.Title.Color = 'y';
    end
end
set(gcf,'Color',[0 0 0])

end


%% Test bursting timing 

units = unique(F(3,:));

Pair = [1 4];

f = F(2,F(3,:)==Pair(1,1)); 
Tr = round((mod(f,spl)-10-sp)/4.8);
T1 = T(ceil(f/spl)) + cast(Tr','uint32');

f = F(2,F(3,:)==Pair(1,2)); 
Tr = round((mod(f,spl)-10-sp)/4.8);
T2 = T(ceil(f/spl)) + cast(Tr','uint32');
 



%% MountainSort v.s. Wclust (confusion matrix?)

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
