//
// Created by jure on 10/11/22.
//

#include <algorithm>
#include <cmath>
#include "Metrics.h"

void evolution::Metrics::addFitness(double fitness) {
    this->checkLocked();
    this->fitnessList.push_back(fitness);
}

void evolution::Metrics::addAccuracy(double accuracy) {
    this->checkLocked();
    this->accuracyList.push_back(accuracy);
}

void evolution::Metrics::addMcc(double mcc) {
    this->checkLocked();
    this->mccList.push_back(mcc);
}

double evolution::Metrics::getBottomFitnessPercentile(double percentile) {
    this->checkNotLocked();
    percentile = evolution::Metrics::fixPercentile(percentile);
    auto index = this->percentageToIndex(percentile);
    double cumulative = this->bottomFitnessList.at(index);
    return cumulative;
}

double evolution::Metrics::getBottomAccuracyPercentile(double percentile) {
    this->checkNotLocked();
    percentile = evolution::Metrics::fixPercentile(percentile);
    auto index = this->percentageToIndex(percentile);
    double cumulative = this->bottomAccuracyList.at(index);
    return cumulative;
}

double evolution::Metrics::getBottomMccPercentile(double percentile) {
    this->checkNotLocked();
    percentile = evolution::Metrics::fixPercentile(percentile);
    auto index = this->percentageToIndex(percentile);
    double cumulative = this->bottomMccList.at(index);
    return cumulative;
}

void evolution::Metrics::lock() {
    std::sort(this->fitnessList.begin(), this->fitnessList.end());
    std::sort(this->accuracyList.begin(), this->accuracyList.end());
    std::sort(this->mccList.begin(), this->mccList.end());

    double bottomCumulative = 0;
    double topCumulative = 0;
    for (unsigned long i = 0; i < this->fitnessList.size(); i++) {
        bottomCumulative += this->fitnessList.at(i);
        topCumulative += this->fitnessList.at(this->fitnessList.size() - 1 - i);
        this->bottomFitnessList.push_back(bottomCumulative / (double) (i + 1));
        this->topFitnessList.push_back(topCumulative / (double) (i + 1));
    }
    bottomCumulative = 0;
    topCumulative = 0;
    for (unsigned long i = 0; i < this->accuracyList.size(); i++) {
        bottomCumulative += this->accuracyList.at(i);
        topCumulative += this->accuracyList.at(this->fitnessList.size() - 1 - i);
        this->bottomAccuracyList.push_back(bottomCumulative / (double) (i + 1));
        this->topAccuracyList.push_back(topCumulative / (double) (i + 1));
    }
    bottomCumulative = 0;
    topCumulative = 0;
    for (unsigned long i = 0; i < this->mccList.size(); i++) {
        bottomCumulative += this->mccList.at(i);
        topCumulative += this->mccList.at(this->fitnessList.size() - 1 - i);
        this->bottomMccList.push_back(bottomCumulative / (double) (i + 1));
        this->topMccList.push_back(topCumulative / (double) (i + 1));
    }
    this->locked = true;
}

double evolution::Metrics::getBestFitness() {
    this->checkNotLocked();
    return this->fitnessList.at(this->fitnessList.size() - 1);
}

double evolution::Metrics::getBestAccuracy() {
    this->checkNotLocked();
    return this->accuracyList.at(this->accuracyList.size() - 1);
}

double evolution::Metrics::getBestMcc() {
    this->checkNotLocked();
    return this->mccList.at(this->mccList.size() - 1);
}

double evolution::Metrics::getWorstFitness() {
    this->checkNotLocked();
    return this->fitnessList.at(0);
}

double evolution::Metrics::getWorstAccuracy() {
    this->checkNotLocked();
    return this->accuracyList.at(0);
}

double evolution::Metrics::getWorstMcc() {
    this->checkNotLocked();
    return this->mccList.at(0);
}

double evolution::Metrics::getTopFitnessPercentile(double percentile) {
    this->checkNotLocked();
    percentile = evolution::Metrics::fixPercentile(percentile);
    unsigned long index = this->percentageToIndex(percentile);
    double cumulative = this->topFitnessList.at(index);
    return cumulative;
}

double evolution::Metrics::getTopAccuracyPercentile(double percentile) {
    this->checkNotLocked();
    percentile = evolution::Metrics::fixPercentile(percentile);
    unsigned long index = this->percentageToIndex(percentile);
    double cumulative = this->topAccuracyList.at(index);
    return cumulative;
}

double evolution::Metrics::getTopMccPercentile(double percentile) {
    this->checkNotLocked();
    percentile = evolution::Metrics::fixPercentile(percentile);
    unsigned long index = this->percentageToIndex(percentile);
    double cumulative = this->topMccList.at(index);
    return cumulative;
}

double evolution::Metrics::getAverageSliceFitness(double bottomStartPercentage, double bottomEndPercentage) {
    this->checkNotLocked();
    bottomStartPercentage = evolution::Metrics::fixPercentile(bottomStartPercentage);
    bottomEndPercentage = evolution::Metrics::fixPercentile(bottomEndPercentage);

    auto startIndex = this->percentageToIndex(bottomStartPercentage);
    if (startIndex > 0) {
        startIndex += 1;
    }
    auto endIndex = this->percentageToIndex(bottomEndPercentage);

    double sum = 0;
    for (unsigned long i = startIndex; i <= endIndex; i++) {
        sum += this->fitnessList.at(i);
    }
    double average = sum / ((double) (endIndex - startIndex + 1));
    return average;
}

double evolution::Metrics::getAverageSliceAccuracy(double bottomStartPercentage, double bottomEndPercentage) {
    this->checkNotLocked();
    bottomStartPercentage = evolution::Metrics::fixPercentile(bottomStartPercentage);
    bottomEndPercentage = evolution::Metrics::fixPercentile(bottomEndPercentage);

    auto startIndex = this->percentageToIndex(bottomStartPercentage);
    if (startIndex > 0) {
        startIndex += 1;
    }
    auto endIndex = this->percentageToIndex(bottomEndPercentage);

    double sum = 0;
    for (unsigned long i = startIndex; i <= endIndex; i++) {
        sum += this->accuracyList.at(i);
    }
    double average = sum / ((double) (endIndex - startIndex + 1));
    return average;
}

double evolution::Metrics::getAverageSliceMcc(double bottomStartPercentage, double bottomEndPercentage) {
    this->checkNotLocked();
    bottomStartPercentage = evolution::Metrics::fixPercentile(bottomStartPercentage);
    bottomEndPercentage = evolution::Metrics::fixPercentile(bottomEndPercentage);

    auto startIndex = this->percentageToIndex(bottomStartPercentage);
    if (startIndex > 0) {
        startIndex += 1;
    }
    auto endIndex = this->percentageToIndex(bottomEndPercentage);

    double sum = 0;
    for (unsigned long i = startIndex; i <= endIndex; i++) {
        sum += this->mccList.at(i);
    }
    double average = sum / ((double) (endIndex - startIndex + 1));
    return average;
}

void evolution::Metrics::checkNotLocked() const {
    if (!locked) {
        throw std::domain_error("This object is not yet locked!");
    }
}

unsigned long evolution::Metrics::percentageToIndex(double percentile) const {
    unsigned long length = this->fitnessList.size();
    auto end = (unsigned long) round((double) length * (percentile / 100));
    if (end == 0) {
        return 0;
    } else {
        return end - 1;
    }
}

double evolution::Metrics::fixPercentile(double percentile) {
    if (percentile > 100) {
        return 100;
    } else if (percentile < 0) {
        return 0;
    } else {
        return percentile;
    }
}

void evolution::Metrics::checkLocked() const {
    if (this->locked) {
        throw std::domain_error("This object is locked!");
    }
}

std::vector<double> evolution::Metrics::getFitnessList()  {
    return this->fitnessList;
}





