//
// Created by jure on 11. 02. 21.
//

#ifndef NEUROEVOLUTION_DATAINSTANCE_H
#define NEUROEVOLUTION_DATAINSTANCE_H

#include <memory>
#include <vector>

namespace data_structures {
    class DataInstance {
    private:
        std::vector<double> values;
        unsigned int correctIndex{};

    public:
        /**
         * Destroys the values vector.
         */
        ~DataInstance() = default;

        /**
         * Creates a new DataInstance object with values.
         * @param valuesAndLabelIndex all values are input values, the last one is the correct index
         */
        explicit DataInstance(std::vector<double> valuesAndLabelIndex) : values(std::move(valuesAndLabelIndex)) {};

        static std::shared_ptr<DataInstance>
        createDataInstance(const std::vector<double> &valuesAndLabelIndex);

        /**
         * Get the vector of input values.
         * @return vector of input values
         */
        std::vector<double> getValues();

        /**
         * Get the index of the correct answer.
         * @return the index of the correct answer
         */
        unsigned int getCorrectIndex();

    };
}


#endif //NEUROEVOLUTION_DATAINSTANCE_H
