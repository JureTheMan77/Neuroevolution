//
// Created by jure on 10/11/22.
//

#include <algorithm>
#include <cmath>
#include "Metrics.h"

void evolution::Metrics::addFitness(double fitness) {
    if (this->locked) {
        throw std::domain_error("This object is locked!");
    }
    this->fitnessList.push_back(fitness);
}

void evolution::Metrics::addAccuracy(double accuracy) {
    if (this->locked) {
        throw std::domain_error("This object is locked!");
    }
    this->accuracyList.push_back(accuracy);
}

void evolution::Metrics::addMcc(double mcc) {
    if (this->locked) {
        throw std::domain_error("This object is locked!");
    }
    this->mccList.push_back(mcc);
}

double evolution::Metrics::getFitnessPercentile(double percentile) {
    if (!this->locked) {
        throw std::domain_error("This object is not yet locked!");
    }

    percentile = fixPercentile(percentile);

    unsigned long length = this->fitnessList.size();
    auto end = (unsigned long) std::round((double) length * (percentile / 100));
    if (end == 0) {
        end = 1;
    }

    double cumulative = this->cumulativeFitnessList.at(end - 1);
    double avg = cumulative / (double) end;

    return avg;
}

double evolution::Metrics::getAccuracyPercentile(double percentile) {
    if (!this->locked) {
        throw std::domain_error("This object is not yet locked!");
    }

    percentile = fixPercentile(percentile);

    unsigned long length = this->accuracyList.size();
    auto end = (unsigned long) std::round((double) length * (percentile / 100));
    if (end == 0) {
        end = 1;
    }

    double cumulative = this->cumulativeAccuracyList.at(end - 1);
    double avg = cumulative / (double) end;

    return avg;
}

double evolution::Metrics::getMccPercentile(double percentile) {
    if (!this->locked) {
        throw std::domain_error("This object is not yet locked!");
    }
    percentile = fixPercentile(percentile);

    unsigned long length = this->mccList.size();
    auto end = (unsigned long) std::round((double) length * (percentile / 100));
    if (end == 0) {
        end = 1;
    }

    double cumulative = this->cumulativeMccList.at(end - 1);
    double avg = cumulative / (double) end;

    return avg;
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

void evolution::Metrics::lock() {
    std::sort(this->fitnessList.begin(), this->fitnessList.end());
    std::sort(this->accuracyList.begin(), this->accuracyList.end());
    std::sort(this->mccList.begin(), this->mccList.end());

    double cumulative = 0;
    for (unsigned long i = 0; i < this->fitnessList.size(); i++) {
        cumulative += this->fitnessList.at(i);
        this->cumulativeFitnessList.push_back(cumulative);
    }
    cumulative = 0;
    for (unsigned long i = 0; i < this->accuracyList.size(); i++) {
        cumulative += this->accuracyList.at(i);
        this->cumulativeAccuracyList.push_back(cumulative);
    }
    cumulative = 0;
    for (unsigned long i = 0; i < this->mccList.size(); i++) {
        cumulative += this->mccList.at(i);
        this->cumulativeMccList.push_back(cumulative);
    }
    this->locked = true;
}

double evolution::Metrics::getBestFitness() {
    if (!this->locked) {
        throw std::domain_error("This object is not yet locked!");
    }
    return this->fitnessList.at(this->fitnessList.size() - 1);
}

double evolution::Metrics::getBestAccuracy() {
    if (!this->locked) {
        throw std::domain_error("This object is not yet locked!");
    }
    return this->accuracyList.at(this->accuracyList.size() - 1);
}

double evolution::Metrics::getBestMcc() {
    if (!this->locked) {
        throw std::domain_error("This object is not yet locked!");
    }
    return this->mccList.at(this->mccList.size() - 1);
}
