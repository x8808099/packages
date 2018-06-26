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
#include "discard_noisy_clusters.h"
#include "get_sort_indices.h"

#include <cmath>
#include <diskreadmda.h>
#include <diskreadmda32.h>
#include <mda.h>
#include <mda32.h>
#include <QVector>
#include "pca.h"
#include "kdtree.h"

using std::fabs;
using std::sqrt;

Mda32 compute_mean_clip(const Mda32& clips);
Mda32 extract_clips(const DiskReadMda32& X, const QVector<double>& times, int clip_size);
Mda32 compute_noise_shape(const Mda32& noise_clips, const Mda32& template0);
QVector<double> sample(const QVector<double>& times, bigint num);
void regress_out_noise_shape(Mda32& clips, const Mda32& shape);
double compute_noise_overlap(const DiskReadMda32& X, const QVector<double>& times, const Mda32& noise_clips0, Discard_noisy_clusters_opts opts);

void discard_noisy_clusters(QVector<double>& times, QVector<int>& labels, QVector<int>& central_channels, Mda32& templates, const DiskReadMda32& X, Discard_noisy_clusters_opts opts)
{   
    int M = templates.N1();
    int T = templates.N2();
    int K = templates.N3();
    bigint L = labels.count();
    bigint N_in =  X.N2() / opts.input_clip_size;
    bigint num_noise_sample = 1000;
    if (N_in <= num_noise_sample)
        return;
  
    QVector<double> noise_times0;
    double increment = N_in * 1.0 / num_noise_sample;
    double j = 0;
    for (int i = 0; i < num_noise_sample; i++) { 
            noise_times0 <<  ((bigint)j + 1) * opts.input_clip_size - opts.clip_size/2 - 1;
            j += increment;
    }
    Mda32 noise_clips0 = extract_clips(X, noise_times0, T);

    // bigint num_noise_sample = qMin(N_in, 1000);
    // qsrand(times.count()); 
    // for (bigint i = 0; i < num_noise_sample; i++) {
    //      noise_times0 << (qrand() % N_in + 1) * opts.input_clip_size - opt.clip_size/2 - 1;
    // }
    
    //#pragma omp parallel for
    // for (int jj = 0; jj < cluster_numbers.count(); jj++) {
    //     DiskReadMda32 X0;
    //     QVector<double> times_k;
    //     int k;
    //     P_isolation_metrics_opts opts0;
    // #pragma omp critical
    //     {
    //         X0 = X;
    //         k = cluster_numbers[jj];
    //         for (bigint i = 0; i < labels.count(); i++) {
    //             if (labels[i] == k)
    //                 times_k << times[i];
    //         }
    //         opts0 = opts;
    //     }

    //     get_cluster_metrics(X0, times_k, opts0);

    QMap<int, int> label_map;
    label_map[0] = 0;
    int knext = 1;

    for (int k = 1; k <= K; k++) {
        label_map[k] = 0;
        QVector<double> time_k;
        bigint det_counts = 0; 
	int central_channel = 0;
	int ii = 0;
        while (central_channel == 0) {
	        if (labels[ii] == k) { 
	            central_channel = central_channels[ii];
            }
            ii++;
        }
        for (bigint i = 0; i < L; i++) {
            if (labels[i]==k) {
                time_k << times[i];
                int det_time = (bigint)(times[i])%opts.input_clip_size;
                if (det_time >= opts.noise_detect_time) {
                    det_counts++; 
                }
            }
        }
        double noise_det_ratio = det_counts * 1.0 / time_k.count();
        if (noise_det_ratio < opts.detect_time_discard_thresh) {
            double noise_overlap = compute_noise_overlap(X, time_k, noise_clips0, opts);
            if (noise_overlap < opts.noise_overlap_discard_thresh) {
                qDebug().noquote() << QString("Ch%1 |%2 ratio:%3 % overlap:%4 % -> C%5").arg(central_channel).arg(time_k.count(),5).arg(noise_det_ratio*100,5,'f',1).arg(noise_overlap*100,6,'f',2).arg(knext);            
                label_map[k] = knext;
                knext++;
            }
       	    else {
                qDebug().noquote() << QString("Ch%1 |%2 ratio:%3 % overlap:%4 %").arg(central_channel).arg(time_k.count(),5).arg(noise_det_ratio*100,5,'f',1).arg(noise_overlap*100,6,'f',2);
            }
        }
        else {
            qDebug().noquote() << QString("Ch%1 |%2 ratio:%3 %").arg(central_channel).arg(time_k.count(),5).arg(noise_det_ratio*100,5,'f',1);
        }
    }

    // finalizing...

    QVector<bigint> inds_to_use;
    for (bigint ii = 0; ii < L; ii++) {
        int ij = labels[ii];
        if (label_map[ij] > 0)
            inds_to_use << ii;
    }

    QVector<double> times2(inds_to_use.count());
    QVector<int> labels2(inds_to_use.count());
    QVector<int> central_channels2(inds_to_use.count());
    for (bigint i = 0; i < inds_to_use.count(); i++) {
        int ik = inds_to_use[i];
        times2[i] = times[ik];
        labels2[i] = label_map[labels[ik]];
        central_channels2[i] = central_channels[ik];
    }
    
    int K2 = MLCompute::max(labels2);
    Mda32 templates2(M, T, K2);
    for (int k2 = 1; k2 <= K2; k2++) {
        int k1 = 0;
        for (int a = 1; a <= K; a++) {
            if (label_map[a] == k2)
                k1 = a;
        }
        if (k1 > 0) {
            Mda32 tmp;
            templates.getChunk(tmp, 0, 0, k1 - 1, M, T, 1);
            templates2.setChunk(tmp, 0, 0, k2 - 1);
        }
    }

    times = times2;
    labels = labels2;
    central_channels = central_channels2;
    templates = templates2;
}

double compute_noise_overlap(const DiskReadMda32& X, const QVector<double>& times, const Mda32& noise_clips0, Discard_noisy_clusters_opts opts)
{
    int T = opts.clip_size;
    bigint N_in =  X.N2() / opts.input_clip_size;
    bigint num_to_use = qMin(opts.max_num_to_use, times.count());

    QVector<double> times_subset = sample(times, num_to_use);
    QVector<bigint> labels_subset;
    for (bigint i = 0; i < times_subset.count(); i++)
        labels_subset << 1;
    
    QVector<double> noise_times;
    QVector<bigint> noise_labels;
     for (bigint i = 0; i < times_subset.count(); i++) {
         bigint ni = (bigint)(times_subset[i] + 0.5) / opts.input_clip_size + 100;
         noise_times <<  (ni % N_in + 1) * opts.input_clip_size - opts.clip_size/2 - 1 ;
         noise_labels << 0;
    }

    QVector<bigint> all_labels = labels_subset;
    QVector<double> all_times = times_subset;
    all_times.append(noise_times);
    all_labels.append(noise_labels);

    Mda32 clips = extract_clips(X, times_subset, T);
    Mda32 template1 = compute_mean_clip(clips);
    Mda32 noise_shape = compute_noise_shape(noise_clips0, template1); 
    
    Mda32 all_clips = extract_clips(X, all_times, T);
    regress_out_noise_shape(all_clips, noise_shape);

    Mda32 all_clips_reshaped(all_clips.N1() * all_clips.N2(), all_clips.N3());
    bigint NNN = all_clips.totalSize();
    for (bigint iii = 0; iii < NNN; iii++)
        all_clips_reshaped.set(all_clips.get(iii), iii);

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
        for (bigint j = 0; j < FF.N1(); j++)
            p << FF.value(j, i);
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

Mda32 extract_clips(const DiskReadMda32& X, const QVector<double>& times, int clip_size)
{   
    bigint M = X.N1();
    bigint T = clip_size;
    bigint Tmid = (bigint)((T + 1) / 2) - 1;
    bigint L = times.count();
    Mda32 clips(M, T, L);
    for (bigint i = 0; i < L; i++) {
        bigint t0 = (bigint)(times[i] - Tmid + 0.5);
        Mda32 tmp;
        if (!X.readChunk(tmp, 0, t0, M, T)) {
            qWarning() << "Problem reading chunk in extract_clips";
        }
        clips.setChunk(tmp, 0, 0, i);
    }
    return clips;
}

// bigint random_time(bigint N_in, bigint clip_size)
// {
//     if (N <= clip_size * 2)
//         return N / 2;
//     return clip_size + (qrand() % (N - clip_size * 2));
// }
