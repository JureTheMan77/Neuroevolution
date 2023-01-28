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
        std::vector<double> bottomFitnessList;
        std::vector<double> topFitnessList;
        std::vector<double> accuracyList;
        std::vector<double> bottomAccuracyList;
        std::vector<double> topAccuracyList;
        std::vector<double> mccList;
        std::vector<double> bottomMccList;
        std::vector<double> topMccList;
        bool locked;

    private:
        static double fixPercentile(double percentile);

        [[nodiscard]] unsigned long percentageToIndex(double percentile) const;

        void checkNotLocked() const;

        void checkLocked() const;

    public:
        void addFitness(double fitness);

        void addAccuracy(double accuracy);

        void addMcc(double mcc);

        double getBottomFitnessPercentile(double percentile);

        double getBottomAccuracyPercentile(double percentile);

        double getBottomMccPercentile(double percentile);

        double getTopFitnessPercentile(double percentile);

        double getTopAccuracyPercentile(double percentile);

        double getTopMccPercentile(double percentile);

        double getBestFitness();

        double getBestAccuracy();

        double getBestMcc();

        double getWorstFitness();

        double getWorstAccuracy();

        double getWorstMcc();

        double getAverageSliceFitness(double bottomStartPercentage, double bottomEndPercentage);

        double getAverageSliceAccuracy(double bottomStartPercentage, double bottomEndPercentage);

        double getAverageSliceMcc(double bottomStartPercentage, double bottomEndPercentage);

        void lock();

        std::vector<double> getFitnessList();
    };
}


#endif //NEUROEVOLUTION_METRICS_H
