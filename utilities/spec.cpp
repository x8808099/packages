/* Define the spec */
inputs_opt('raw','filt','pre','geom');
outputs('firings_out');
outputs_opt('filt_out','pre_out');
outputs_opt('cluster_metrics_out','label_map_out','features_out');

param('samplerate',48000);
param('freq_min',300);
param('freq_max',6000);
param('freq_wid',1000);
param('whiten','true');
param('merge_across_channels','true');

param('adjacency_radius',-1);
param('detect_sign',1);
param('detect_threshold',3);
param('detect_interval',15);
param('clip_size',60);
param('consolidation_factor',0.9);

param('firing_rate_thresh',0.05);
param('isolation_thresh',0.95);
param('noise_overlap_thresh',0.03);
param('peak_snr_thresh',1.5);      

param('num_features',10);      

_Pipeline.run=function(X) {
  var pp=X.parameters;
  
  var pre='pre';
  if (!X.hasInput('pre')) {
    
    var filt='filt';
    if (!X.hasInput('filt')) {
      if (!X.hasInput('raw')) {
        console.error('Missing input: raw, filt or pre');
        return -1;
      }
      X.runProcess('ms3.bandpass_filter',
                   {timeseries:'raw'},
                   {timeseries_out:'filt_out'},
                   {samplerate:pp.samplerate,freq_min:pp.freq_min,freq_max:pp.freq_max,freq_wid:pp.freq_wid}
                  );
      filt='filt_out';
    }
  
  
    if (pp.whiten=='true') {
      X.runProcess('ms3.whiten',
                   {timeseries:filt},
                   {timeseries_out:'pre_out'},
                   {}
                  );
    }
    else {
      X.runProcess('pyms.normalize_channels',
                   {timeseries:filt},
                   {timeseries_out:'pre_out'},
                   {}
                  );
    }
    pre='pre_out';
  }

  var p={
    detect_threshold:pp.detect_threshold,
    detect_sign:pp.detect_sign,
    detect_interval:pp.detect_interval,
    adjacency_radius:pp.adjacency_radius,
    clip_size:pp.clip_size,
    merge_across_channels:pp.merge_across_channels,
    consolidation_factor:pp.consolidation_factor
  };
  
  var inputs={timeseries:pre};
  if (X.hasInput('geom')) {
    inputs.geom='geom';
  }
  X.runProcess('mountainsortalg.ms3alg',
               inputs,
               {firings_out:'firings_out'},
               p);
  
  var pc={
    clip_size:pp.clip_size,
    samplerate:pp.samplerate,
    isolation_thresh:pp.isolation_thresh,
    noise_overlap_thresh:pp.noise_overlap_thresh, 
    firing_rate_thresh:pp.firing_rate_thresh,
    peak_snr_thresh:pp.peak_snr_thresh,
    num_features:pp.num_features
  };

  X.runPipeline('curate',
               {pre:pre,firings:'firings_out'},
               {cluster_metrics_out:'cluster_metrics_out',
                label_map_out:'label_map_out',
                features_out:'features_out'},
                pc);

};

/////////////////////////////////////////////////////////////////////

function param(str,val) {
      if (val===undefined) {
        _Pipeline.spec.parameters.push({name:str});
      }
      else {
        _Pipeline.spec.parameters.push({name:str,optional:true,default_value:val});
      }
}
                
function inputs(str1,str2,str3,str4) {
  if (str1) _Pipeline.spec.inputs.push({name:str1});
  if (str2) _Pipeline.spec.inputs.push({name:str2});
  if (str3) _Pipeline.spec.inputs.push({name:str3});
  if (str4) _Pipeline.spec.inputs.push({name:str4});
}

function inputs_opt(str1,str2,str3,str4) {
  if (str1) _Pipeline.spec.inputs.push({name:str1,optional:true});
  if (str2) _Pipeline.spec.inputs.push({name:str2,optional:true});
  if (str3) _Pipeline.spec.inputs.push({name:str3,optional:true});
  if (str4) _Pipeline.spec.inputs.push({name:str4,optional:true});
}

function outputs(str1,str2,str3,str4) {
  if (str1) _Pipeline.spec.outputs.push({name:str1});
  if (str2) _Pipeline.spec.outputs.push({name:str2});
  if (str3) _Pipeline.spec.outputs.push({name:str3});
  if (str4) _Pipeline.spec.outputs.push({name:str4});
}

function outputs_opt(str1,str2,str3,str4) {
  if (str1) _Pipeline.spec.outputs.push({name:str1,optional:true});
  if (str2) _Pipeline.spec.outputs.push({name:str2,optional:true});
  if (str3) _Pipeline.spec.outputs.push({name:str3,optional:true});
  if (str4) _Pipeline.spec.outputs.push({name:str4,optional:true});
}