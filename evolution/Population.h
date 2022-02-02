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
        unsigned int edgeTraverseLimit;
        unsigned int maxDeepVertices;
        unsigned int maxEdges;
        double maxMutationChance;
        std::vector<std::shared_ptr<data_structures::DataInstance>> trainingValues;
        std::vector<std::shared_ptr<data_structures::DataInstance>> testingValues;

        unsigned int const UINT_MAX = std::numeric_limits<unsigned int>::max();
    public:

        explicit Population(const std::string &pathToDataSet);

        ~Population() = default;

        /**
         * Initialize a new population.
         * @param populationSizeArg maximum number of agents
         * @param maxVertices maximum number of vertices an agent can have
         * @param maxEdgesArg maximum number of edges an agent can have
         * @param keepDormantVextexesAndEdges keep vertices that don't have an out connection
         */
        void
        initialisePopulation(unsigned int populationSizeArg, unsigned int maxDeepVerticesArg, unsigned int maxEdgesArg,
                             unsigned int edgeTraverseLimitArg, bool keepDormantVerticesAndEdges,
                             double maxMutationChanceArg);

        /**
         * Calculate the fitness of the population.
         */
        void calculateFitness(double vertexContribution, double edgeContribution);

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

        const std::vector<std::shared_ptr<evolution::Agent>> &getPopulation() const;

        std::shared_ptr<evolution::Agent> getFittestAgent();

        double getAverageFitness();

        void addAgent(std::shared_ptr<evolution::Agent> agent);

        std::vector<std::shared_ptr<data_structures::DataInstance>> getTestingValues();

        unsigned int getNumberOfOutputs();

        std::vector<std::string> getOutputLabels();

        std::string toString();

    private:
        std::shared_ptr<evolution::Agent> createAgent(bool keepDormantVerticesAndEdges);

        /**
         * Adapted from @see https://en.wikipedia.org/wiki/Stochastic_universal_sampling
         * @param agentsToKeep number of agents to keep
         */
        std::vector<unsigned int> stochasticUniversalSampling(unsigned int agentsToKeep);

        std::shared_ptr<data_structures::DeepVertex>
        getDominantVertex(const std::shared_ptr<data_structures::DeepVertex> &v1,
                          const std::shared_ptr<data_structures::DeepVertex> &v2);

        std::shared_ptr<data_structures::Edge> getDominantEdge(const std::shared_ptr<data_structures::Edge> &e1,
                                                               const std::shared_ptr<data_structures::Edge> &e2);

        unsigned int getDominantCrossoverable(const std::shared_ptr<data_structures::ICrossoverable> &c1,
                                              const std::shared_ptr<data_structures::ICrossoverable> &c2);

        std::vector<std::shared_ptr<data_structures::DeepVertex>>
        createDeepVertexChildren(const std::shared_ptr<data_structures::DeepVertex> &vertex);

        void crossoverDeepVertices(const std::shared_ptr<evolution::Agent> &childAgent,
                                   const std::shared_ptr<evolution::Agent> &parent1,
                                   const std::shared_ptr<evolution::Agent> &parent2,
                                   std::vector<unsigned int> &checkedIndices);


        void
        crossoverEdges(std::shared_ptr<evolution::Agent> &childAgent, const std::shared_ptr<evolution::Agent> &parent1,
                       const std::shared_ptr<evolution::Agent> &parent2, std::vector<unsigned int> &checkedIndices);

        void createAndAddEdgeChildren(const std::shared_ptr<data_structures::Edge> &edge,
                                      std::shared_ptr<evolution::Agent> &childAgent);

        void addNewRandomEdges(std::shared_ptr<evolution::Agent> const &childAgent);

        data_structures::ICrossoverable
        mutateCrossoverable(const std::shared_ptr<data_structures::ICrossoverable> &crossoverable);

        unsigned int childrenToProduce(const std::shared_ptr<data_structures::ICrossoverable> &crossoverable);

        void addRandomEdge(unsigned int index, std::shared_ptr<data_structures::Graph> const &graph);

        void calculateAgentFitness(double vertexContribution, double edgeContribution,
                                   const std::shared_ptr<evolution::Agent> &agent);

        std::shared_ptr<evolution::Agent> crossoverThreaded();
    };
}


#endif //NEUROEVOLUTION_POPULATION_H
