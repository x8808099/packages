cd('~/ss')
addpath(genpath(pwd))

doTetrodeData=1;
doEEG=0; doTetrodeTS=0; doPosition=0; doSync=0; doInput=0;  
drIn = 'SpikeSort/'; 
nchs=4;
%%
bpf_name = '0615R16BC-cl.bpf';

[tetrodeData,tetrodeTimestamp,tetrodeUnit,tetrodeChannel,sl,sample_rate] = ...
    bpf2mat_py([drIn bpf_name],doEEG,doTetrodeData,doTetrodeTS,doPosition,doSync,doInput);    
%%
bpf_name = 'M16_2.bpf';
load(['/media/psf/PH/DATA/' bpf_name(1:end-4) '.mat'],'');
Chs = unique(tetrodeChannel)
%%
sl = tetrodeSamplesPerBlock; sp=24;
iCh = 1;
det_itv = 30;
clip_size = 60;
noth = 10; % noise_overlap_threshold 
csof = 99;
merge_or_not = 'true';
whiten_or_not = 'true';
det_th = 2.5;

Chs = unique(tetrodeChannel)
Data = reshape(tetrodeData,nchs,sl,[]);
idxs = tetrodeChannel==Chs(iCh);
U = tetrodeUnit(idxs); 
T = tetrodeTimestamp(idxs); 
[min(diff(T)) unique(U')]
dtac = reshape(Data(:,:,idxs),nchs,[]);

%C = dtac*dtac'/size(dtac,2); [V,D] = eig(C); W = V/sqrt(D)*V'   

dtac = dtac-mean(dtac,2);
dta = reshape(dtac,nchs,sl,[]);
dta = cat(2,zeros(nchs,24,size(dta,3)),dta);
dtac = reshape(dta,nchs,[]);

name_ = [bpf_name(1:end-4),'_s',num2str(clip_size),'i',num2str(det_itv),...
    't',num2str(det_th),'C',num2str(iCh),'_',num2str(max(U),'%02d')];
firing_name = '/media/psf/PH/firing.mda';
metric_name = ['/media/psf/PH/M_' name_ '.json'];

% pre-whitening
% dtah = reshape(dtac,4,sl,[]);
% dtah(:,1:48,:)=[];
% dtah = reshape(dtah,4,[]);
% nn = size(dtah,2);
% C = dtah*dtah'/nn;
% [V,D] = eig(C);
% W = V/sqrt(D)*V';
% dtac = W*dtac;
%%
mda_name = 'ss.mda'; 
writemda(dtac,mda_name,'float32')

system(['mlp-run ~/mountainlab/mountainsort3.mlp sort'...
    ' --clip_size=' num2str(clip_size)...
    ' --detect_interval=' num2str(det_itv)...
    ' --noise_overlap_thresh=' num2str(noth/100)...
    ' --detect_threshold=' num2str(det_th)...
    ' --consolidation_factor=' num2str(csof/100)...
    ' --filt=' mda_name...
    ' --pre_out=/media/psf/PH/pre.mda'...
    ' --label_map_out=/media/psf/PH/lmap.mda'... % ' --features_out=fets.mda'...
    ' --cluster_metrics_out=' metric_name...
    ' --whiten=' whiten_or_not...
    ' --merge_across_channels=' merge_or_not...
    ' --firings_out=' firing_name]);
command = ['mountainview'...
    ' --raw=/media/psf/PH/pre.mda'...
    ' --samplerate=48000'... 
    ' --cluster_metrics=' metric_name...
    ' --firings=' firing_name];
clipboard('copy',command)
%%
% /home/parallels/ss/
command = ['mountainview'...
    ' --raw=/media/psf/PH/pre.mda'...
    ' --samplerate=48000'... 
    ' --cluster_metrics=' metric_name...
    ' --firings=' firing_name];
clipboard('copy',command)

%%
lmap=readmda('/media/psf/PH/lmap.mda');
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
