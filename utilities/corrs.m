clear
cd('/Users/hupeiyao/Movies/PH')
addpath(genpath(pwd))

drIn = '/Users/hupeiyao/Movies/PH/Fs/';
Files = dir([drIn '*.mda']);
mat_name = 'name0'; 

bS = 8; % bin size for average
T_valid = 3000;
cohs_m=[]; cors_m=[];
cohs_w=[]; cors_w=[];

for n_file = 1:length(Files)

firing_name = Files(n_file).name;
mat_name_new = [firing_name(3:end-10) '.mat'];
F = readmda([firing_name]);
    
if ~strcmp(mat_name_new, mat_name)
    mat_name = mat_name_new;
    load(['DATA/' mat_name]);
    sl = tetrodeSamplesPerBlock;
    nchs = 4; sp = 24; spl = sl+sp;
    Chs = unique(tetrodeChannel);
    Data = reshape(tetrodeData,nchs,sl,[]);
end

if str2double(firing_name(end-5:end-4))==0; continue; end

iCh = str2double(firing_name(end-7));
idxs = tetrodeChannel==Chs(iCh);
T = tetrodeTimestamp(idxs);
U = tetrodeUnit(idxs);
% dtac = reshape(Data(:,:,idxs),nchs,[]);

for Wc = 0:1 

if Wc % Wclust
    keepCluster = U;
else % Wc = 0 to do MS
    keepCluster = F(3,:);
end
    
units = unique(keepCluster);
if isempty(units); continue; end
if units(1)==0; units(1)=[]; end % exclude noise cluster
if isempty(units); continue; end

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
        Tr = round((mod(f,spl)-15-sp)/4.8);
        Tu = T(ceil(f/spl)) + cast(Tr','uint32');
    end
    [counts(u,:),~] = histcounts(Tu,roomTimeStamps(times));
end

for i = times(1:end-1)
    x = double(roomXY(i,1)); y = double(roomXY(i,2)); 
    if x*y==0; continue; end
    A(ceil(x/bS),ceil(y/bS),1) = A(ceil(x/bS),ceil(y/bS),1) + roomTimeStamps(i)-roomTimeStamps(i-1);
    for u = 1:length(units)
        A(ceil(x/bS),ceil(y/bS),u+1) = A(ceil(x/bS),ceil(y/bS),u+1) + counts(u,i-times(1)+1);
    end
end

PT = A(:,:,1);
P_valid = ~(PT==0);

for nplot = 1:size(A,3)-1
    PS = A(:,:,nplot+1);
    
    if or(sum(PS(:))<100, sum(PS(:))>10000); continue; end
    % only analyze clusters with 100~10000 events

    % place_field = place_spike_counts / time_spent 
    % about ~30x30 (already binned in 8x8 pixels)
    PF = PS./PT;
    % average over 8 surrounding point
    kernel = [1 1 1; 1 0 1; 1 1 1]/8; 
    PS_s = conv2(PS,kernel,'same'); 
    PT_s = conv2(PT,kernel,'same');
    % firing rate of surrounding
    PF_s = PS_s./PT_s; 
    % only use the data with at least 0.3s
    to_use = and(PT>T_valid,~(PT_s==0));
    % coherence
    coh = corr(PF(to_use),PF_s(to_use));
    % correlation of spike-count and time
    cor = corr(PS(PT>T_valid),PT(PT>T_valid));   
    
    if Wc==0
        cohs_m = [cohs_m, coh];
        cors_m = [cors_m, cor];
    else
        cohs_w = [cohs_w, coh];
        cors_w = [cors_w, cor];
    end
end

end
end
%%
bin_width=0.08;fac = 1.5;
figure;hold on
histogram(cohs_w-fac*cors_w,'BinWidth',bin_width,'Normalization','probability'); 
histogram(cohs_m-fac*cors_m,'BinWidth',bin_width,'Normalization','probability'); 
legend({'Wclust','MountainSort'})

%% this 
bin_width=0.05;
figure;hold on
histogram(cohs_w,'BinWidth',bin_width,'Normalization','probability'); 
histogram(cohs_m,'BinWidth',bin_width,'Normalization','probability'); 
legend({'Wclust','MountainSort'})

%%
 
bin_width=0.05;
figure;hold on
histogram(cors_w,'BinWidth',bin_width,'Normalization','probability'); 
histogram(cors_m,'BinWidth',bin_width,'Normalization','probability'); 
legend({'Wclust','MountainSort'})

%%
bin_width=0.04;
figure;hold on
histogram(cohs_w,'BinWidth',bin_width) 
histogram(cohs_m,'BinWidth',bin_width)
legend({'Wclust','MountainSort'})












