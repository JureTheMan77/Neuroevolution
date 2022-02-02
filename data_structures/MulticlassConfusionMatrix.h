//
// Created by jure on 1. 02. 22.
//

#ifndef NEUROEVOLUTION_MULTICLASSCONFUSIONMATRIX_H
#define NEUROEVOLUTION_MULTICLASSCONFUSIONMATRIX_H

#include <vector>
#include "../evolution/Agent.h"

namespace data_structures {
    class MulticlassConfusionMatrix {
    private:
        unsigned int numberOfClasses;
        std::vector<std::vector<unsigned int>> matrix;
        double accuracy;
        double matthewsCorrelationCoefficient;

    public:
        MulticlassConfusionMatrix(const std::shared_ptr<evolution::Agent> &agent,
                                  const std::vector<std::shared_ptr<data_structures::DataInstance>> &testingSet,
                                  unsigned int numOfClasses);

        double getAccuracy();

        double getMatthewsCorrelationCoefficient();

        std::string toString(std::vector<std::string> labels);

    private:
        void calculateAccuracy();

        /**
         * source: https://en.wikipedia.org/wiki/Phi_coefficient#cite_note-gorodkin2004comparing-32
         */
        void calculateMatthewsCorrelationCoefficient();


    };
}


#endif //NEUROEVOLUTION_MULTICLASSCONFUSIONMATRIX_H
