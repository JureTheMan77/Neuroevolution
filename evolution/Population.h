//
// Created by jure on 1. 02. 21.
//

#ifndef NEUROEVOLUTION_POPULATION_H
#define NEUROEVOLUTION_POPULATION_H


#include <random>
#include "Agent.h"
#include "../enums/SelectionType.h"

namespace evolution {
    class Population {
    private:
        unsigned int populationSize;
        std::vector<std::shared_ptr<evolution::Agent>> population;
        std::vector<std::shared_ptr<evolution::Agent>> populationPlaceholder;
        std::random_device seeder;
        std::vector<std::string> inputLabels;
        std::vector<std::string> outputLabels;
        unsigned int numberOfInputs;
        unsigned int numberOfOutputs;
        std::vector<std::shared_ptr<data_structures::DataInstance>> trainingValues;
    public:

        explicit Population(const std::string &pathToDataSet);

        ~Population() = default;

        /**
         * Initialize a new population.
         * @param populationSize maximum number of agents
         * @param maxVertices maximum number of vertices an agent can have
         * @param maxEdges maximum number of edges an agent can have
         * @param keepDormantVextexesAndEdges keep vertices that don't have an out connection
         */
        void initialisePopulation(unsigned int populationSize, unsigned int maxDeepVertices, unsigned int maxEdges, unsigned int edgeTraverseLimit,
                                  bool keepDormantVerticesAndEdges);

        /**
         * Calculate the fitness of the popoulation.
         */
        void calculateFitness();

        /**
         * Sample the population.
         * @param type selection type
         * @param agentsToKeep number of agents to keep
         */
        void sample(enums::SelectionType type, unsigned int agentsToKeep);

        /**
         * Selects random parents from the population and
         */
        void crossover();

        std::string toString();

    private:
        std::shared_ptr<evolution::Agent>
        createAgent(unsigned int maxDeepVertices, unsigned int maxEdges, unsigned int edgeTraverseLimit, bool keepDormantVerticesAndEdges);

        /**
         * Adapted from @see https://en.wikipedia.org/wiki/Stochastic_universal_sampling
         * @param agentsToKeep number of agents to keep
         */
        std::vector<unsigned int> stochasticUniversalSampling(unsigned int agentsToKeep);

        std::shared_ptr<data_structures::Vertex> getDominantVertex(std::shared_ptr<data_structures::Vertex> v1, std::shared_ptr<data_structures::Vertex> v2);
    };
}


#endif //NEUROEVOLUTION_POPULATION_H
