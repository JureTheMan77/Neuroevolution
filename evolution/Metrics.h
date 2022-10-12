//
// Created by jure on 10/11/22.
//

#ifndef NEUROEVOLUTION_METRICS_H
#define NEUROEVOLUTION_METRICS_H


#include <vector>

namespace evolution {
    class Metrics {
    private:
        std::vector<double> fitnessList;
        std::vector<double> cumulativeFitnessList;
        std::vector<double> accuracyList;
        std::vector<double> cumulativeAccuracyList;
        std::vector<double> mccList;
        std::vector<double> cumulativeMccList;
        bool locked;

    private:
        static double fixPercentile(double percentile);

    public:
        void addFitness(double fitness);
        void addAccuracy(double accuracy);
        void addMcc(double mcc);

        double getFitnessPercentile(double percentile);
        double getAccuracyPercentile(double percentile);
        double getMccPercentile(double percentile);

        double getBestFitness();
        double getBestAccuracy();
        double getBestMcc();

        void lock();
    };
}


#endif //NEUROEVOLUTION_METRICS_H
