/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef P_MOUNTAINSORT3_H
#define P_MOUNTAINSORT3_H

#include <QString>
typedef int64_t bigint;

struct P_mountainsort3_opts {
    
    // ~ indicates parameters that does not depend on values here (set in pipeline) 
    
    double adjacency_radius = -1;// ~
    double detect_threshold = 3;// ~
    int detect_interval = 30;//~
    int detect_sign = 1;//~

    int clip_size = 60;//~
    int clip_shift = 6;//~

    int clip_padding = 24;

    int num_features = 10;//~
    int num_features_per_channel = 10;//~

    bigint max_pca_samples = 10000;

    bool consolidate_clusters = true;
    double consolidation_factor = 0.95;//~

    bool discard_noisy_clusters = true;
    bool merge_across_channels = true;//~
    double event_fraction_threshold = 0.5;//~

    int input_clip_size = 120;//~
    int noise_detect_time = 66;//~
    double detect_time_discard_thresh = 0.5;//~
    double noise_overlap_discard_thresh = 0.05;//~

    bool fit_stage = true;//~

    double t1 = -1;
    double t2 = -1;
};

bool p_mountainsort3(QString timeseries, QString geom, QString firings_out, QString temp_path, const P_mountainsort3_opts& opts);

#endif // P_MOUNTAINSORT3_H
