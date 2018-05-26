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
#ifndef DISCARD_NOISY_CLUSTERS_H
#define DISCARD_NOISY_CLUSTERS_H

#include <QVector>
#include "mda32.h"
#include "diskreadmda32.h"

struct Discard_noisy_clusters_opts {
    int clip_size = 60;
    int max_num_to_use = 500;
    int min_num_to_use = 100;
    int num_features = 10;
    int K_nearest = 6;
    int exhaustive_search_num = 100;
    int input_clip_size = 120; 
    int noise_detect_time = 48;
    double detect_time_discard_thresh = 0.5;
    double noise_overlap_discard_thresh = 0.1;
};

void discard_noisy_clusters(QVector<double>& times, QVector<int>& labels, QVector<int>& central_channels, Mda32& templates, const DiskReadMda32& X, Discard_noisy_clusters_opts opts);

#endif //DISCARD_NOISY_CLUSTERS_H

