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
#include "get_sort_indices.h"
#include <cmath>
#include <QVector>
#include <mda.h>
#include <mda32.h>
#include "pca.h"
#include "kdtree.h"

using std::fabs;
using std::sqrt;

bool should_use_template(const Mda32& template0, Consolidate_clusters_opts opts);
Mda32 extract_clips(const Mda32& X, const QVector<double>& times, int clip_size);
Mda32 compute_noise_shape(const Mda32& noise_clips, const Mda32& template0);
void regress_out_noise_shape(Mda32& clips, const Mda32& shape);
Mda32 compute_mean_clip(const Mda32& clips);
double compute_noise_overlap(const Mda32& X, const QVector<double>& times, Consolidate_clusters_opts opts);

QMap<int, int> consolidate_clusters(QVector<double>& times, QVector<int>& labels, const Mda32& templates, const Mda32& X, Consolidate_clusters_opts opts)
{   
    bigint L = labels.count();
    int M = templates.N1();
    int T = templates.N2();
    int K = MLCompute::max(labels);
    QVector<int> to_use(K + 1, 0);

    for (int k = 1; k <= K; k++) {
        Mda32 template0;
        templates.getChunk(template0, 0, 0, k - 1, M, T, 1);
        if (should_use_template(template0, opts)) {
            to_use[k] = 1;
        }
    }

    QMap<int, int> label_map;
    label_map[0] = 0;
    int knext = 1;

    for (int k = 1; k <= K; k++) {
        label_map[k] = 0;
        if (to_use[k]) {
            QVector<double> timek;
            bigint det_counts = 0; 
            for (bigint i = 0; i < L; i++) {
                if (labels[i]==k) {
                    timek << times[i];
                    int det_time = (bigint)(times[i])%120;
                    if (det_time>47)
                        det_counts++; 
                    }
            }
            double noise_det_ratio = det_counts * 1.0 / timek.count();
            printf("Noise detection ratio is %f \n", noise_det_ratio);
            if (noise_det_ratio<0.5) {
                double noise_overlap = compute_noise_overlap(X, timek, opts);
                printf("Overlap is %f \n", noise_overlap);
                if (noise_overlap<0.3) {
                    label_map[k] = knext;
                    knext++;
                }
            }
        } 
    }
    return label_map;
}

bool should_use_template(const Mda32& template0, Consolidate_clusters_opts opts)
{
    int central_channel = 0;

    int peak_location_tolerance = 15;

    bigint M = template0.N1();
    bigint T = template0.N2();
    bigint Tmid = (int)((T + 1) / 2) - 6;

    double abs_peak_on_central_channel = 0;
    bigint abs_peak_t_on_central_channel = Tmid;
    {
        for (bigint t = 0; t < T; t++) {
            double val = template0.value(central_channel, t);
            if (val > abs_peak_on_central_channel) {            //  fabs(val)
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
        if (fabs(abs_peak_t_on_central_channel - Tmid) > peak_location_tolerance)
            return false;
        // if (fabs(abs_peak_t - Tmid) > peak_location_tolerance)
          //  return false;
    }
    return true;
}

bigint random_time(bigint N, bigint clip_size)
{
    if (N <= clip_size * 2)
        return N / 2;
    return clip_size + (qrand() % (N - clip_size * 2));
}

QVector<double> sample(const QVector<double>& times, bigint num)
{
    QVector<double> random_values(times.count());
    for (bigint i = 0; i < times.count(); i++) {
        random_values[i] = sin(i * 12 + i * i);
    }
    QList<bigint> inds = get_sort_indices_bigint(random_values);
    QVector<double> ret;
    for (bigint i = 0; (i < num) && (i < inds.count()); i++) {
        ret << times[inds[i]];
    }
    return ret;
}

double compute_noise_overlap(const Mda32& X, const QVector<double>& times, Consolidate_clusters_opts opts)
{

    bigint N = X.N2();
    int T = opts.clip_size;
    bigint num_to_use = qMin(opts.max_num_to_use, times.count());

    QVector<double> times_subset = sample(times, num_to_use);
    QVector<bigint> labels_subset;
    for (bigint i = 0; i < times_subset.count(); i++) {
        labels_subset << 1;
    }
    QVector<bigint> all_labels = labels_subset; //0 and 1
    QVector<double> all_times = times_subset;
    
    QVector<double> noise_times;
    QVector<bigint> noise_labels;
    for (bigint i = 0; i < times_subset.count(); i++) {
        noise_times << random_time(N, T);
        noise_labels << 0;
    }
    all_times.append(noise_times);
    all_labels.append(noise_labels);

 
    Mda32 clips = extract_clips(X, times_subset, T);
    Mda32 template1 = compute_mean_clip(clips);
    
    QVector<double> noise_times1;
    for (bigint i = 0; i < 1000; i++) {
        noise_times1 << random_time(N, T);
    }
    Mda32 noise_clips = extract_clips(X, noise_times1, T);
    Mda32 noise_shape = compute_noise_shape(noise_clips, template1); 
    
    Mda32 all_clips = extract_clips(X, all_times, T);
    regress_out_noise_shape(all_clips, noise_shape);

    Mda32 all_clips_reshaped(all_clips.N1() * all_clips.N2(), all_clips.N3());
    bigint NNN = all_clips.totalSize();
    for (bigint iii = 0; iii < NNN; iii++) {
        all_clips_reshaped.set(all_clips.get(iii), iii);
    }

    Mda32 FF;
    Mda32 CC, sigma;
    bool subtract_mean = false;
    pca(CC, FF, sigma, all_clips_reshaped, opts.num_features, subtract_mean);

    KdTree tree;
    tree.create(FF);
    double num_correct = 0;
    double num_total = 0;
    for (bigint i = 0; i < FF.N2(); i++) {
        QVector<float> p;
        for (bigint j = 0; j < FF.N1(); j++) {
            p << FF.value(j, i);
        }
        QList<int> indices = tree.findApproxKNearestNeighbors(FF, p, opts.K_nearest, opts.exhaustive_search_num);
        for (bigint a = 0; a < indices.count(); a++) {
            if (indices[a] != i) {
                if (all_labels[indices[a]] == all_labels[i])
                    num_correct++;
                num_total++;
            }
        }
    }

    if (!num_total)
        return 0;
    return 1 - (num_correct * 1.0 / num_total);
}


Mda32 extract_clips(const Mda32& X, const QVector<double>& times, int clip_size)
{
    bigint M = X.N1();
    bigint T = clip_size;
    bigint Tmid = (bigint)((T + 1) / 2) - 1;
    bigint L = times.count();
    
    Mda32 clips;
    clips.allocate(M, T, L);
    for (bigint i = 0; i < L; i++) {
        bigint t1 = times.value(i) - Tmid;
    for (int t = 0; t < T; t++) {
        for (int m = 0; m < M; m++) {
                clips.set(X.value(m, t1), m, t, i);
            }
        }
    }
    return clips;
}

Mda32 compute_mean_clip(const Mda32& clips)
{
    bigint M = clips.N1();
    bigint T = clips.N2();
    bigint L = clips.N3();
    Mda32 ret;
    ret.allocate(M, T);
    bigint aaa = 0;
    for (bigint i = 0; i < L; i++) {
        bigint bbb = 0;
        for (bigint t = 0; t < T; t++) {
            for (bigint m = 0; m < M; m++) {
                ret.set(ret.get(bbb) + clips.get(aaa), bbb);
                aaa++;
                bbb++;
            }
        }
    }
    if (L) {
        for (bigint t = 0; t < T; t++) {
            for (bigint m = 0; m < M; m++) {
                ret.set(ret.get(m, t) / L, m, t);
            }
        }
    }
    return ret;
}

Mda32 compute_noise_shape(const Mda32& noise_clips, const Mda32& template0)
{
    bigint peak_channel = 0;
    bigint peak_timepoint = 0;
    double best_val = 0;
    for (bigint t = 0; t < template0.N2(); t++) {
        for (bigint m = 0; m < template0.N1(); m++) {
            double val = template0.value(m, t);
            if (val > best_val) {
                best_val = val;
                peak_channel = m;
                peak_timepoint = t;
            }
        }
    }

    Mda ret(template0.N1(), template0.N2());
    for (bigint i = 0; i < noise_clips.N3(); i++) {
        double weight = noise_clips.value(peak_channel, peak_timepoint, i) / noise_clips.N3();
        for (bigint t = 0; t < template0.N2(); t++) {
            for (bigint m = 0; m < template0.N1(); m++) {
                ret.set(ret.get(m, t) + noise_clips.get(m, t, i) * weight, m, t);
            }
        }
    }
    
    Mda32 ret2(ret.N1(), ret.N2());
    for (bigint i = 0; i < ret.totalSize(); i++) {
        ret2.set(ret.get(i), i);
    }
    return ret2;
}


void regress_out_noise_shape(Mda32& clips, const Mda32& shape)
{
    double shape_norm = MLCompute::norm(shape.totalSize(), shape.constDataPtr());
    for (bigint i = 0; i < clips.N3(); i++) {
        Mda32 clip;
        clips.getChunk(clip, 0, 0, i, clips.N1(), clips.N2(), 1);
        double inner_product = MLCompute::dotProduct(clip.totalSize(), clip.constDataPtr(), shape.constDataPtr());
        if (shape_norm) {
            double coef = inner_product / (shape_norm * shape_norm);
            for (bigint j = 0; j < clip.totalSize(); j++) {
                clip.set(clip.get(j) - coef * shape.get(j), j);
            }
        }
        clips.setChunk(clip, 0, 0, i);
    }
}
