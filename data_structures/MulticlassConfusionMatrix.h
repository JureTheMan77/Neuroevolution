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
        unsigned int numberOfClasses = 0;
        std::vector<std::vector<unsigned int>> matrix;
        double accuracy = 0;
        double matthewsCorrelationCoefficient = 0;

    public:
        /**
         * Create a new multiclass confusion matrix.
         * @param agent agent to test
         * @param testingSet data to test with
         * @param numOfClasses number of classes in the testing set
         */
        MulticlassConfusionMatrix(const std::shared_ptr<evolution::Agent> &agent,
                                  const std::vector<std::shared_ptr<data_structures::DataInstance>> &testingSet,
                                  unsigned int numOfClasses);

        /**
         * Default destructor.
         */
        MulticlassConfusionMatrix() = default;

        /**
         * Get the calculated accuracy.
         * @return accuracy
         */
        [[nodiscard]] double getAccuracy() const;

        /**
         * Get the calculated mcc.
         * @return mcc
         */
        [[nodiscard]] double getMatthewsCorrelationCoefficient() const;

        /**
         * Get a pretty string representation of this object.
         * @param labels dataset labels
         * @return pretty string representation
         */
        std::string toString(std::vector<std::string> labels);

    private:
        /**
         * Calculate the accuracy.
         */
        void calculateAccuracy();

        /**
         * Calculate the mcc.
         * source: https://en.wikipedia.org/wiki/Phi_coefficient#cite_note-gorodkin2004comparing-32
         */
        void calculateMatthewsCorrelationCoefficient();
    };
}


#endif //NEUROEVOLUTION_MULTICLASSCONFUSIONMATRIX_H
