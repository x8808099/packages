clear
cd('~/ss')
addpath(genpath(pwd))

%mat_name = 'M9_052214_blackbox-cl.mat';
mat_name = 'M16_2.mat';
% mat_name = '0615R16BC-cl.mat';

load(['/media/psf/PH/DATA/' mat_name]);
Chs = unique(tetrodeChannel)

sl = tetrodeSamplesPerBlock; % sample length
nchs=4; % tetrode
Data = reshape(tetrodeData,nchs,sl,[]);

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

iCh = 1;
sp=24; % padding

idxs = tetrodeChannel==Chs(iCh);
dtac = reshape(Data(:,:,idxs),nchs,[]);
U = tetrodeUnit(idxs); unique(U') % T = tetrodeTimestamp(idxs); min(diff(T)) 
dtac = dtac - mean(dtac,2); % subtract mean before padding with zeros
dta = reshape(dtac,nchs,sl,[]); dta = cat(2,zeros(nchs,sp,size(dta,3)),dta);
dtac = reshape(dta,nchs,[]); dtac = [dtac,zeros(nchs,sp)];
%dtac = [dtac(2:4,:);dtac(1,:)]; % actually affect the result...
run_count = 1;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

tic
clip_size = 64; clip_shift = 8;
input_clip_size = sp + sl;
det_th = 2.8; % event detection threshold in stdev
det_itv = 30; % event detection blocking period

% For discarding noisy clusters
ndt = sp + 12 + det_itv; % noisy detection time
td_th = 0.5; % allowed ratio of timestamps after ndt
nd_th = 0.05; % noise overlap threshold
oe_th = 0.5; % overlapping events fraction threshold

discard_or_not = 'true';
merge_or_not = 'true';
fit_or_not = merge_or_not;
whitening = 'true'; %also need change th
extract_features = 0; 

name_ = [mat_name(1:end-4),'_s',num2str(clip_size),'i',num2str(det_itv),...
    't',num2str(det_th),'C',num2str(iCh),'_',num2str(max(U),'%02d')];
firing_name = '/media/psf/PH/firing.mda';
metric_name = '/media/psf/PH/metric.json';

% C = dtac*dtac'/size(dtac,2); [V,D] = eig(C); W = V/sqrt(D)*V'   
% pre-whitening
% dX = [1:24, 72:120];
% dta(:,dX,:)=[];
% dta = reshape(dta,4,[]);
% nn = size(dta,2);
% C = dta*dta'/nn;
% [V,D] = eig(C);
% W = V/sqrt(D)*V';
% dtac = W*dtac;

mda_name = 'ss8.mda'; 
writemda(dtac,mda_name,'float32');

command_sort = ['mlp-run'...
    ' ~/mountainlab/mountainsort3.mlp sort'...
    ' --clip_shift=' num2str(clip_shift)...
    ' --clip_size=' num2str(clip_size)...
    ' --clip_padding=' num2str(sp)...
    ' --detect_interval=' num2str(det_itv)...
    ' --detect_threshold=' num2str(det_th)...
    ' --noise_detect_time=' num2str(ndt)...
    ' --detect_time_discard_thresh=' num2str(td_th)...
    ' --noise_overlap_discard_thresh=' num2str(nd_th)...
    ' --input_clip_size=' num2str(input_clip_size)...
    ' --event_fraction_threshold=' num2str(oe_th)...
    ' --consolidation_factor=0.95'...
    ' --label_map_out=/media/psf/PH/lmap.mda'... 
    ' --whiten=' whitening...
    ' --discard_noisy_clusters=' discard_or_not...
    ' --merge_across_channels=' merge_or_not...
    ' --fit_stage=' fit_or_not...
    ' --cluster_metrics_out=' metric_name...
    ' --firings_out=' firing_name];

if run_count == 1
    command_sort = [command_sort ' --filt=' mda_name ' --pre_out=pre.mda'];
else
    command_sort = [command_sort ' --pre=pre.mda'];
end
run_count = run_count + 1;    

if (extract_features)
    command_sort = [command_sort ' --extract_fets=true --features_out=fets.mda'];
    system(command_sort);

lmap = readmda('/media/psf/PH/lmap.mda');
fets = readmda('fets.mda');
F = readmda(firing_name);
keepCluster = zeros(1,size(F,2));
for i = 1:size(lmap,1)
    keepCluster(F(3,:)==lmap(i,2)) = lmap(i,1);
end
Nkeep = sum(unique(keepCluster)>0);

[nj,ni] = size(fets);
fid = fopen('fets.txt','w');
fprintf(fid,'Cluster f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 \n');
for i = 1:ni
    fprintf(fid,'%d\t',keepCluster(i));
    for j = 1:nj
        fprintf(fid,'%f\t',fets(j,i));
    end
    fprintf(fid,'\n');
end
fclose(fid);
isoi_name = ['I_' name_ '_' num2str(Nkeep,'%02d') '.txt'];
system(['isoi fets.txt /media/psf/PH/Is/' isoi_name]);

else
    system(command_sort);
end

command_view = ['mountainview'...
    ' --raw=/home/parallels/ss/pre.mda'...
    ' --samplerate=48000'... 
    ' --cluster_metrics=' metric_name...
    ' --firings=' firing_name];
clipboard('copy',command_view)

toc