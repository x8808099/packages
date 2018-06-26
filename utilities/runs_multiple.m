clear
cd('~/ss')
addpath(genpath(pwd))

Files = dir('/media/psf/PH/DATA/*.mat');

% Pipeline options, curation = 'false'; extract_features = 'fasle';
whitening = 'true'; % need to lower threshold if not
merge_or_not = 'true';
fit_or_not = 'true';

clip_size = 60; clip_shift = 6; sp = 24;
det_th = 2.8; % event detection threshold in stdev
det_itv = clip_size/2; % event detection blocking period

% For discarding noisy clusters
noise_detect_time = sp + 12 + det_itv; % 12 is the expected peak time inside the clip
td_th = 0.5; % allowed ratio of timestamps after noise detection time
nd_th = 0.05; % noise overlap threshold
oe_th = 0.5; % overlapping events fraction threshold


%%
tic
run_count = 1;

for n_file = 1:length(Files)  

mat_name = Files(n_file).name;
load(['/media/psf/PH/DATA/' mat_name]);
sl = tetrodeSamplesPerBlock;
spl=sp+sl; nchs=4;
Chs = unique(tetrodeChannel);
Data = reshape(tetrodeData,nchs,sl,[]);

for iCh = 1:length(Chs)
run_count = run_count + 1;
mda_name = ['ss' num2str(mod(run_count,2)) '.mda'];

idxs = tetrodeChannel==Chs(iCh);
U = tetrodeUnit(idxs); % T = tetrodeTimestamp(idxs); min(diff(T))
firing_name = ['/media/psf/PH/Fs_2/F_' mat_name(1:end-4) '_C' num2str(iCh) '_' num2str(max(U),'%02d') '.mda'];
% if exist(firing_name,'file')==2; continue; end 

data = reshape(Data(:,:,idxs),nchs,[]); 
data = data-mean(data,2); % subtract mean before padding with zeros
dta = reshape(data,nchs,sl,[]); dta = cat(2,zeros(nchs,sp,size(dta,3)),dta);
data = reshape(dta,nchs,[]); data = [data,zeros(nchs,sp)];

if size(data,2)<1200000; continue; end
% only sort data with at least 10000 clips

writemda(data,mda_name,'float32');
    
command_sort = ['mlp-run'...
    ' ~/mountainlab/mountainsort3.mlp sort'...
    ' --clip_size=' num2str(clip_size)...
    ' --clip_shift=' num2str(clip_shift)...
    ' --detect_interval=' num2str(det_itv)...
    ' --detect_threshold=' num2str(det_th)...
    ' --detect_time_discard_thresh=' num2str(td_th)...
    ' --noise_overlap_discard_thresh=' num2str(nd_th)... 
    ' --noise_detect_time=' num2str(noise_detect_time)...
    ' --event_fraction_threshold=' num2str(oe_th)...
    ' --consolidation_factor=0.95'...
    ' --whiten=' whitening...
    ' --curation=false'...
    ' --extract_fets=false'...
    ' --merge_across_channels=' merge_or_not...
    ' --fit_stage=' fit_or_not...
    ' --filt=' mda_name...
    ' --firings_out=' firing_name];

system(command_sort);

end
end

disp('Done!')
toc