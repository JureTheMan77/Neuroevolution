//
// Created by jure on 1. 02. 22.
//

#include "MulticlassConfusionMatrix.h"
#include <cmath>
#include <sstream>
#include <iomanip>

data_structures::MulticlassConfusionMatrix::MulticlassConfusionMatrix(const std::shared_ptr<evolution::Agent> &agent,
                                                                      const std::vector<std::shared_ptr<data_structures::DataInstance>> &testingSet,
                                                                      unsigned int numOfClasses) {
    this->numberOfClasses = numOfClasses;
    // create the confusion matrix based on the testing data and the provided agent
    // rows: correct classes
    // columns: predicted classes
    std::vector<std::vector<unsigned int>> confusionMatrix(this->numberOfClasses + 1);
    // initialize rows to a fixed size
    for (int i = 0; i < this->numberOfClasses + 1; i++) {
        confusionMatrix.at(i) = std::vector<unsigned int>(this->numberOfClasses + 1);
    }
    this->matrix = confusionMatrix;

    // run through all the testing values
    for (const std::shared_ptr<data_structures::DataInstance> &di: testingSet) {
        // reset the testing agent
        agent->getGraph()->reset();

        agent->getGraph()->traverse(di);

        // check if the prediction is correct
        unsigned int predictedIndex = agent->getGraph()->getLargestOutputValueIndex();
        // increment the correct cell
        this->matrix.at(di->getCorrectIndex()).at(predictedIndex) += 1;

        // reset the agent
        // agent->getGraph()->reset();
    }

    // sum the rows
    for (int i = 0; i < this->numberOfClasses; i++) {
        //sum the values of each row
        unsigned int sum = 0;
        for (int j = 0; j < this->numberOfClasses; j++) {
            sum += this->matrix.at(i).at(j);
        }
        // write the sum in the last column
        this->matrix.at(i).at(this->numberOfClasses) = sum;
    }

    // sum the columns
    for (int i = 0; i < this->numberOfClasses + 1; i++) {
        unsigned int sum = 0;
        for (int j = 0; j < this->numberOfClasses; j++) {
            sum += this->matrix.at(j).at(i);
        }
        // write the sum in the last row
        this->matrix.at(this->numberOfClasses).at(i) = sum;
    }

    // all output vertices have at least one input edge
    bool noInput = false;
    for(const auto &vertex : agent->getGraph()->getOutputVertices()){
        if(vertex->getInputEdges().empty()){
            noInput = true;
            break;
        }
    }

    if(noInput){
        this->accuracy = 0;
        this->matthewsCorrelationCoefficient = -1;
    } else {

        // calculate accuracy
        this->calculateAccuracy();

        // calculate the Matthews correlation coefficient
        this->calculateMatthewsCorrelationCoefficient();
    }
}

void data_structures::MulticlassConfusionMatrix::calculateAccuracy() {
// sum the diagonal
    double acc = 0;
    for (int i = 0; i < this->numberOfClasses; i++) {
        acc += this->matrix.at(i).at(i);
    }
    unsigned int total = this->matrix.at(this->numberOfClasses).at(this->numberOfClasses);
    acc = acc / total;
    this->accuracy = acc;
}

void data_structures::MulticlassConfusionMatrix::calculateMatthewsCorrelationCoefficient() {
    // the total number of samples correctly predicted
    // essentially, sum the diagonal, except the last value
    double c = 0;
    for (int i = 0; i < this->numberOfClasses; i++) {
        c += (double) this->matrix.at(i).at(i);
    }

    // the total number of samples
    auto s = (double) this->matrix.at(this->numberOfClasses).at(this->numberOfClasses);
    double s2 = s * s;

    double tp = 0;
    double tt = 0;
    double pp = 0;
    for (int i = (int) this->numberOfClasses - 1; i >= 0; i--) {
        tp = tp + (double) (this->matrix.at(this->numberOfClasses).at(i) *
                            this->matrix.at(i).at(this->numberOfClasses));

        pp = pp + (double) (this->matrix.at(this->numberOfClasses).at(i) *
                            this->matrix.at(this->numberOfClasses).at(i));

        tt = tt + (double) (this->matrix.at(i).at(this->numberOfClasses) *
                            this->matrix.at(i).at(this->numberOfClasses));
    }

    double numerator = c * s - tp;
    double denominator = std::sqrt(s2 - pp) * std::sqrt(s2 - tt);
    double mcc = -1;

    // denominator can be 0 if all predictions are for the same class
    if (denominator != 0) {
        mcc = numerator / denominator;
    }

    this->matthewsCorrelationCoefficient = mcc;

}

double data_structures::MulticlassConfusionMatrix::getAccuracy() const {
    return this->accuracy;
}

double data_structures::MulticlassConfusionMatrix::getMatthewsCorrelationCoefficient() const {
    return this->matthewsCorrelationCoefficient;
}

std::string data_structures::MulticlassConfusionMatrix::toString(std::vector<std::string> labels) {
    // print the matrix
    // get the longest output label
    int longestOutputLabel = 0;
    for (auto const &label: labels) {
        if (longestOutputLabel < label.length()) {
            longestOutputLabel = (int) label.length();
        }
    }
    // extra space
    longestOutputLabel += 1;

    std::ostringstream stream;
    stream << "\n";

    for (int i = 0; i < numberOfClasses + 1; i++) {
        if (i == 0) {
            // print labels
            for (int j = 0; j < numberOfClasses + 1; j++) {
                if (j == 0) {
                    // empty space
                    stream << std::setw(longestOutputLabel) << "" << " ";
                }
                if (j == numberOfClasses) {
                    // sum column
                    stream << std::setw(7) << "sum" << " ";
                } else {
                    stream << labels.at(j) << " ";
                }

            }
            stream << "\n";
        }
        for (int j = 0; j < numberOfClasses + 1; j++) {
            if (j == 0) {
                // print label
                if (i == numberOfClasses) {
                    stream << std::setw(longestOutputLabel) << "sum" << " ";
                } else {
                    stream << std::setw(longestOutputLabel) << labels.at(i) << " ";
                }
            }
            if (j == numberOfClasses) {
                stream << std::setw(7) << std::fixed << std::setprecision(3) << this->matrix.at(i).at(j);
            } else {
                stream << std::setw((int) labels.at(j).length()) << std::fixed << std::setprecision(3)
                       << this->matrix.at(i).at(j) << " ";
            }
        }
        stream << "\n";
    }
    stream << "Accuracy: " << std::fixed << std::setprecision(3) << this->accuracy << "\n";
    stream << "Matthews correlation coefficient: " << std::fixed << std::setprecision(3)
           << this->matthewsCorrelationCoefficient;

    return stream.str();
}
