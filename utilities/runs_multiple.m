clear
cd('~/ss')
addpath(genpath(pwd))

%%
doTetrodeData=1;
doEEG=0; doTetrodeTS=0; doPosition=0; doSync=0; doInput=0;  
drIn = 'SpikeSort/'; List = dir([drIn '*.bpf']);
nchs=4;

noth = 3; % noise_overlap_threshold 
det_itv = 15;
clip_size = 60;
runcount = 0;

for nfile = 1:length(List)
    
bpf_name = List(nfile).name;
[tetrodeData,tetrodeTimestamp,tetrodeUnit,tetrodeChannel,sl,sample_rate] = ...
    bpf2mat_py([drIn bpf_name],doEEG,doTetrodeData,doTetrodeTS,doPosition,doSync,doInput);    

Chs = unique(tetrodeChannel); 
Data = reshape(tetrodeData,nchs,sl,[]);

for iCh = 1:length(Chs)
runcount=runcount+1;

idxs = tetrodeChannel==Chs(iCh);
U = tetrodeUnit(idxs); %T = tetrodeTimestamp(idxs); [min(diff(T)) unique(U')]
name_ = [bpf_name(1:end-4),'_s',num2str(clip_size),'_i',...
    num2str(det_itv),'_o',num2str(noth),'_C',num2str(iCh),'_',num2str(max(U),'%02d')];

dtac = reshape(Data(:,:,idxs),nchs,[]);
if size(dtac,2)<200000; continue; end

mda_name = ['ss' num2str(runcount) '.mda']; firing_name = ['Fs/F_' name_ '.mda'];
writemda(dtac,mda_name,'float32')

system(['mlp-run ~/mountainlab/mountainsort3.mlp sort'...
    ' --label_map_out=lmap.mda'...
    ' --features_out=fets.mda'...
    ' --clip_size=' num2str(clip_size)...
    ' --detect_interval=' num2str(det_itv)...
    ' --noise_overlap_thresh=' num2str(noth/100)...
    ' --filt=' mda_name...
    ' --firings_out=' firing_name]);
    
lmap = readmda('lmap.mda');
firing = readmda(firing_name);
fets = readmda('fets.mda');
delete(mda_name)

keepCluster = zeros(1,size(fets,2));
for j = 1:size(lmap,1)
    keepCluster(firing(3,:)==lmap(j,2)) = lmap(j,1);
end

Nkeep = length(unique(keepCluster))-1;
if Nkeep==0; continue; end

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
disp(runcount);
isoi_name = ['I_' name_ '_' num2str(Nkeep,'%02d') '.txt'];
system(['isoi fets.txt /media/psf/PH/Is/' isoi_name]);
% fid = fopen('~/isoi/iso.txt');isoi=reshape(fscanf(fid,'%f %f %d\n'),3,[]); isoi'

end

end

disp('Done!')