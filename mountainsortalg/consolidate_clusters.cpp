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

#include "consolidate_clusters.h"
#include <cmath>

bool should_use_template(const Mda32& template0, Consolidate_clusters_opts opts);

QMap<int, int> consolidate_clusters(const Mda32& templates, Consolidate_clusters_opts opts)
{   
    int M = templates.N1();
    int T = templates.N2();
    int K = templates.N3();
    QVector<int> to_use(K + 1, 0);

    for (int k = 1; k <= K; k++) {
        Mda32 template0;
        templates.getChunk(template0, 0, 0, k - 1, M, T, 1);
        if (should_use_template(template0, opts))
            to_use[k] = 1;
    }

    QMap<int, int> label_map;
    label_map[0] = 0;
    int knext = 1;
    for (int k = 1; k <= K; k++) {
        if (to_use[k]) {
            label_map[k] = knext;        
            knext++;
        }
        else
            label_map[k] = 0; 
    }
    return label_map;
}

bool should_use_template(const Mda32& template0, Consolidate_clusters_opts opts)
{
    int central_channel = 0;

    int peak_location_tolerance = 15;

    bigint M = template0.N1();
    bigint T = template0.N2();
    bigint Tmid = (int)((T + 1) / 2) - 1 - opts.shifting;

    double abs_peak_on_central_channel = 0;
    bigint abs_peak_t_on_central_channel = Tmid;
    {
        for (bigint t = 0; t < T; t++) {
            double val = template0.value(central_channel, t);
            if (val > abs_peak_on_central_channel) {
                abs_peak_on_central_channel = val;
                abs_peak_t_on_central_channel = t;
            }
        }
    }

    double abs_peak = 0;
    bigint abs_peak_t = Tmid;
    for (bigint t = 0; t < T; t++) {
        for (bigint m = 0; m < M; m++) {
            double val = template0.value(m, t);
            if (val > abs_peak) {
                abs_peak = val;
                abs_peak_t = t;
            }
        }
    }
    {
        if (abs_peak_on_central_channel < abs_peak * opts.consolidation_factor)
            return false;
        // if (fabs(abs_peak_t_on_central_channel - Tmid) > peak_location_tolerance)
            // return false;
        if (fabs(abs_peak_t - Tmid) > peak_location_tolerance)
            return false;
    }
    return true;
}
