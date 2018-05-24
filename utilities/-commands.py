prlctl enter e542210f-e5dd-4a98-8a82-2e6fb455792c
prlctl start 
prlctl stop ..


mlp-run /home/parallels/ms/mountainsort3.mlp sort --_params=/media/psf/PH/params.json --pre_out=/media/psf/PH/pre.mda --firings_out=/media/psf/PH/firing.mda --label_map_out=/media/psf/PH/lmap.mda --features_out=/media/psf/PH/fets.mda --cluster_metrics_out=/media/psf/PH/metrics.json --filt=/media/psf/PH/a.mda




mountainview --raw=/media/psf/PH/pre.mda --samplerate=48000 --cluster_metrics=/media/psf/PH/metrics.json --firings=/media/psf/PH/firing.mda

mp-run-process pyms.extract_clips --timeseries=/media/psf/PH/pre.mda --firings=/media/psf/PH/firing.mda --clips_out=/media/psf/PH/clips.mda --clip_size=60

# caculating IsoI

/home/parallels/ms/isoitools/isoi /media/psf/PH/isoi/fets.txt /media/psf/PH/isoi/IsoI.txt

mountainview --samplerate=48000 --cluster_metrics=/media/psf/PH/M_0615R16BC-cl_s60_i15_o3_C1_06.json --raw=/media/psf/PH/pre21.mda --firings=/media/psf/PH/firing21.mda
mountainview --samplerate=48000 --cluster_metrics=/media/psf/PH/M_0615R16BC-cl_s60_i15_o3_C1_06.json --raw=/media/psf/PH/pre13.mda --firings=/media/psf/PH/firing13.mda


#

gedit ~/.bashrc
export PATH=$PATH:/home/parallels/mountainlab/packages/mountainview/bin
export PATH=$PATH:/home/parallels/mountainlab/packages/mlpipeline/bin
export PATH=$PATH:/home/parallels/mountainlab/bin


mlp-larinet
mp-spec
mp-list-processors

apt update
apt upgrade mountainlab mlpipeline mountainsort mountainview

