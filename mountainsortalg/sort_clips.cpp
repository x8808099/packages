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
#include "sort_clips.h"

#include <QDir>
#include <QTime>
#include "pca.h"
#include "isosplit5.h"

namespace P_sort_clips {
QVector<int> sort_clips_subset(const Mda32& clips, const QVector<bigint>& indices, Sort_clips_opts opts);
}

QVector<int> sort_clips(const Mda32& clips, const Sort_clips_opts& opts)
{
    bigint M = clips.N1();
    bigint T = clips.N2();
    bigint L = clips.N3();

    // qDebug().noquote() << QString("Sorting reduced clips %1x%2x%3").arg(M).arg(T).arg(L);

    QVector<bigint> indices(L);
    for (bigint i = 0; i < L; i++) {
        indices[i] = i;
    }

    QVector<int> labels = P_sort_clips::sort_clips_subset(clips, indices, opts);

    return labels;
}

namespace P_sort_clips {

QVector<int> sort_clips_subset(const Mda32& clips, const QVector<bigint>& indices, Sort_clips_opts opts)
{
    bigint M = clips.N1();  //4
    bigint T = clips.N2(); //10
    bigint L0 = indices.count();

    Mda32 FF;
    {
        // do this inside a code block so memory gets released
        Mda32 clips_reshaped(M * T, L0);
        bigint iii = 0;
        for (bigint ii = 0; ii < L0; ii++) {
            for (bigint t = 0; t < T; t++) {
                for (bigint m = 0; m < M; m++) {
                    clips_reshaped.set(clips.value(m, t, indices[ii]), iii);
                    iii++;
                }
            }
        }

        Mda32 CC, sigma;
        pca_subsampled(CC, FF, sigma, clips_reshaped, opts.num_features, false, opts.max_samples); //should we subtract the mean?
    }

    isosplit5_opts i5_opts;
    i5_opts.isocut_threshold = opts.isocut_threshold;
    i5_opts.K_init = opts.K_init;

    QVector<int> labels0(L0);
    {
        QTime timer;
        timer.start();
        if (!isosplit5(labels0.data(), opts.num_features, L0, FF.dataPtr(), i5_opts)) {

            qWarning() << "Problem in isosplit5 while processing clips from: "+opts.debug_description;
            QString debug_clips_file = QDir::tempPath()+"/debug_isosplit5_clips.mda";
            QString debug_features_file = QDir::tempPath()+"/debug_isosplit5_features.mda";

            qWarning() << "For debug purposes, writing clips to "+debug_clips_file+" and features to "+debug_features_file+" but note that the features are for a subset of the clips.";
            clips.write32(debug_clips_file);
            FF.write32(debug_features_file);

            qWarning() << "Isosplit5 returned with an error. See above for debug information. Aborting";
            abort();
        }
        //qDebug().noquote() << QString("Time elapsed for isosplit (%1x%2) - K=%3: %4 sec").arg(FF.N1()).arg(FF.N2()).arg(MLCompute::max(labels0)).arg(timer.elapsed() * 1.0 / 1000);
    }

    bigint K0 = MLCompute::max(labels0);
    if (K0 <= 1) {  //only 1 cluster
        //qDebug().noquote() << QString("samples %1,%2,%3,%4,%5,%6,%7,%8,%9").arg(labels0[0]).arg(labels0[2]).arg(labels0[4]).arg(labels0[9]).arg(labels0[11]).arg(labels0[22]).arg(labels0[23]).arg(labels0[25]).arg(labels0[29]);
        return labels0;
    }
    else {
        //branch method
        QVector<int> labels_new(L0);
        bigint k_offset = 0;
        for (bigint k = 1; k <= K0; k++) {
            QVector<bigint> indices2, indices2b;
            for (bigint a = 0; a < L0; a++) {
                if (labels0[a] == k) {
                    indices2 << indices[a];
                    indices2b << a;
                }
            }
            if (indices2.count() > 0) {
                QVector<int> labels2 = sort_clips_subset(clips, indices2, opts);
                bigint K2 = MLCompute::max(labels2);
                for (bigint bb = 0; bb < indices2.count(); bb++) {
                    labels_new[indices2b[bb]] = k_offset + labels2[bb];
                }
                k_offset += K2;
            }
        }

        return labels_new;
    }
}
Mda32 compute_templates(Mda32& clips, const QVector<int>& labels)
{
    int M = clips.N1();
    int T = clips.N2();
    bigint L = clips.N3();
    int Kmax = MLCompute::max(labels);
    Mda sum(M, T, Kmax);
    QVector<double> counts(Kmax);
    double* sumptr = sum.dataPtr();
    float* cptr = clips.dataPtr();
    bigint aa = 0;
    for (bigint i = 0; i < L; i++) {
        int k = labels[i];
        if (k > 0) {
            bigint bb = M * T * (k - 1);
            for (bigint jj = 0; jj < M * T; jj++) {
                sumptr[bb] += cptr[aa];
                aa++;
                bb++;
            }
            counts[k - 1]++;
        }
    }
    Mda32 ret(M, T, Kmax);
    for (int k = 0; k < Kmax; k++) {
        double factor = 0;
        if (counts[k])
            factor = 1.0 / counts[k];
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                ret.setValue(sum.value(m, t, k) * factor, m, t, k);
            }
        }
    }
    return ret;
}
}

