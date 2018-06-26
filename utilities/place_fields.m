clear
cd('/Users/hupeiyao/Movies/PH')
addpath(genpath(pwd))

%% !! Read Data to MAT !!
drIn = '/Users/hupeiyao/Movies/PH/DATA/'; 
files = dir([drIn '*.bpf']);
for i = 1:length(files)
    bpf_name = files(i).name
    bpf2matPF([drIn bpf_name]);
end

%bpf2matPF('/Users/hupeiyao/Movies/PH/DATA/M28_111914_blackbox-cl.bpf')

%% load data
%mat_name='M9_052214_blackbox-cl.mat';
%mat_name = '0615R16BC-cl.mat';

mat_name = 'M16_2.mat';
load(['DATA/' mat_name]);

nchs=4;
sl = tetrodeSamplesPerBlock;
Chs = unique(tetrodeChannel)' 
Data = reshape(tetrodeData,nchs,sl,[]);
%%

iCh = 1; idxs = tetrodeChannel==Chs(iCh); T = tetrodeTimestamp(idxs); U = tetrodeUnit(idxs);
round(roomTimeStamps(end-1)/10000)
%% -- Plot Peak Distribution --
F = readmda('/Users/hupeiyao/Movies/PH/firing.mda');
%lmap = readmda('/Users/hupeiyao/Movies/PH/lmap.mda');
% F = readmda('Fs/F_M9_052214_blackbox-cl_C3_11.mda');
unitss = unique(F(3,:)); lmap = [unitss' unitss'];


% resc = 2; lmap(resc,1)=resc;
% lmap(:,1) = lmap(:,2);
% disc=[7]; lmap(disc,1)=0;



figure; 
Nplot=length(lmap); sp=24; spl=sl+sp; clip_shift = 8; trigger_point = 10;
edges = min(mod(F(2,:),spl)):max(mod(F(2,:),spl))+1; edges = edges - sp - clip_shift;
for iu = 1:Nplot
    subplot(4,ceil(Nplot/4),iu); hold on
    idx = F(3,:)==iu;
    Tu = F(2,idx);
    Tr = mod(Tu,spl) - sp - clip_shift;
        if lmap(iu,1)==0
            histogram(Tr,edges,'FaceColor','r','FaceAlpha',0.25);title(['('  num2str(size(Tu,2)) ') ' 'reject' num2str(iu)])
        % elseif lmap(iu,1)~=lmap(iu,2)
            % histogram(Tr,edges,'FaceColor','b','FaceAlpha',0.25);title(['('  num2str(size(Tu,1)) ') ' 'merge' num2str(lmap(iu,1))])
        else
            histogram(Tr,edges); title(['('  num2str(size(Tu,2)) ') ' 'keep' num2str(iu)])
        end
    y_lim = get(gca,'Ylim');
    plot((trigger_point+1)*[1,1], y_lim, 'm', 'LineWidth',1);
    ylim(y_lim)
    xlim([0 edges(end)])
end

%% -- Plot Place Fields --
iCh = 1; idxs = tetrodeChannel==Chs(iCh); T = tetrodeTimestamp(idxs); U = tetrodeUnit(idxs);
%%
F = readmda('/Users/hupeiyao/Movies/PH/firing.mda');
lmap = readmda('/Users/hupeiyao/Movies/PH/lmap.mda');
lmap(:,1) = lmap(:,2);
%rescue = [16 17]; lmap(rescue,1) = rescue;
%my_map = [7 8]; lmap(my_map,1) = my_map(1);
%unitss = unique(F(3,:)); lmap = [unitss' unitss'];

s_title = [num2str(2.8) ' ' num2str(round(roomTimeStamps(end-1)/10000)) 's'];

figure; pl1 = 4; pl2 = 7;
masking = 1;
sorting = 0;
smooing = 1;

bS = 8; % bin size for average
sigma = 1; % for smooting place field
T_valid = 3000;

Color  = {'y','w'};
time_offset = trigger_point + clip_shift + sp - 1; % -1 since time start from zero.
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
A = zeros(Ax, Ay, length(units) + 1);

times = 3:length(roomTimeStamps);
%roomXY(times(1),1)*roomXY(times(1),2)

counts = zeros(length(units),size(times,2)-1);

for u = 1:length(units) 
    if Wc
        Tu = T(keepCluster==units(u));
    else
        f = F(2,keepCluster==units(u)); 
        Tr = round((mod(f,spl)-time_offset)/4.8);
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
    if x*y==0; continue; end
    A(ceil(x/bS),ceil(y/bS),1) = A(ceil(x/bS),ceil(y/bS),1) + roomTimeStamps(i)-roomTimeStamps(i-1);
    for u = 1:length(units)
        A(ceil(x/bS),ceil(y/bS),u+1) = A(ceil(x/bS),ceil(y/bS),u+1) + counts(u,i-times(1)+1);
    end
end

% -- Plot Place Fields --
PT = A(:,:,1);
P_valid = ~(PT==0);
kernel = [1 1 1; 1 0 1; 1 1 1]/8;
rowNeed = ceil(sum(lmap(:,1)>0)/pl2);

subplot(pl1,pl2, 2*(Wc+1)*pl2);
imagesc(PT,'AlphaData',P_valid); axis off; %axis equal;
title(s_title,'color','c','fontsize',13);

for nplot = 1:size(A,3)-1
    subplot(pl1,pl2,nplot+Wc*rowNeed*pl2)
    PS = A(:,:,nplot+1);
    PF = PS./PT;  % place field = spike counts / times
    
    % surround
    PS_s = conv2(PS,kernel,'same'); 
    PT_s = conv2(PT,kernel,'same');
    PF_s = PS_s./PT_s;
    to_use = and(PT>T_valid,~(PT_s==0));
    coh = corr(PF(to_use),PF_s(to_use));
    cor = corr(PS(to_use),PT(to_use));
    
    if smooing
        % filtered
        PS_f = imgaussfilt(PS,sigma);
        PT_f = imgaussfilt(PT,sigma);
        PF = PS_f./PT_f;   
        colormap hot
    else
        colormap default
    end
     
    if masking
        imagesc(PF,'AlphaData',P_valid); 
    else
        imagesc(PF);
    end
    
    %colorbar;
    %axis equal
    axis off
    ax=gca;
    C_Lim = ax.CLim;
    Title = [sprintf('%.2f',cor) '|' sprintf('%.2f',coh) '<' num2str(units(sid((nplot)))) '>' sprintf('%.1f',C_Lim(2)*10000)];
    title([Title '|' num2str(sum(counts(nplot,:)))],'fontsize',13,'color',Color{Wc+1});
    if coh-1.5*cor>0
        ax.Title.Color = 'g';
    end
    if C_Lim(2)<1e-4
        ax.CLim = [0 1e-4];
    end
end
set(gcf,'Color',[0 0 0])

end

%% Test timing of bursting pairs

Pair = [7 8];
bin_width = 2; % 1ms
num_show = 200;

f = F(2,F(3,:)==Pair(1,1)); 
Tr = round((mod(f,spl)-time_offset)/4.8);
T1 = T(ceil(f/spl)) + cast(Tr','uint32');
T1 = double(T1)/10;

f = F(2,F(3,:)==Pair(1,2)); 
Tr = round((mod(f,spl)-time_offset)/4.8);
T2 = T(ceil(f/spl)) + cast(Tr','uint32');
T2 = double(T2)/10; % in ms

[T_sorted,sort_id] = sort([T1;T2]);
T_labels = [zeros(size(T1));ones(size(T2))];
T_diff = [diff(T_sorted) diff(T_labels(sort_id))];
T2_1 = T_diff(T_diff(:,2)==-1); % 2 followed by 1
T1_2 = T_diff(T_diff(:,2)==1);  % 1 followed by 2

figure;hold on
histogram(diff(T1),'BinWidth',bin_width,'Normalization','probability'); 
histogram(diff(T2),'BinWidth',bin_width,'Normalization','probability'); 
xlim([-num_show*bin_width/100 num_show*bin_width]);legend()
%%
bin_width = 1; % 1ms
num_show = 400;
figure;hold on
histogram(T2_1,'BinWidth',bin_width,'Normalization','probability'); 
histogram(T1_2,'BinWidth',bin_width,'Normalization','probability'); 
xlim([-num_show*bin_width/100 num_show*bin_width])

%% MountainSort v.s. Wclust (confusion matrix?)

kCluster = zeros(1,size(F,2));
for j = 1:size(lmap,1)
    kCluster(F(3,:)==lmap(j,2)) = lmap(j,1);
end

count_edges = 0:spl:size(dtac,2);
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
