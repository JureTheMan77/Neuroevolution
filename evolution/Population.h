//
// Created by jure on 1. 02. 21.
//

#ifndef NEUROEVOLUTION_POPULATION_H
#define NEUROEVOLUTION_POPULATION_H


#include <random>
#include "Agent.h"
#include "../enums/SelectionType.h"
#include "../enums/FitnessMetric.h"
#include "Metrics.h"

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
        unsigned int edgeTraverseLimit{};
        unsigned int maxDeepVertices{};
        unsigned int maxEdges{};
        double maxMutationChance{};
        bool keepDormantVerticesAndEdges{};
        std::vector<std::shared_ptr<data_structures::DataInstance>> trainingValues;
        std::vector<std::shared_ptr<data_structures::DataInstance>> testingValues;

        unsigned int const UINT_MAX = std::numeric_limits<unsigned int>::max();
        std::uniform_int_distribution<unsigned int> mutationDistribution{0, 6};
        std::uniform_int_distribution<unsigned int> edgeTraverseLimitDistribution;
        std::uniform_int_distribution<unsigned int> inputVerticesDistribution;
        std::uniform_int_distribution<unsigned int> outputVerticesDistribution;

        std::shared_ptr<evolution::Agent> fittestAgent;
        double averageFittness{0};
        double averageAccuracy{0};
        double averageMatthewsCorrelationCoefficient{0};

    public:

        /**
         * Population constructor with path to dataset.
         * @param pathToDataSet path
         */
        explicit Population(const std::string &pathToDataSet);

        /**
         * Default destructor.
         */
        ~Population() = default;

        /**
         * Initialize a new population.
         * @param populationSizeArg maximum number of agents
         * @param maxVertices maximum number of vertices an agent can have
         * @param maxEdgesArg maximum number of edges an agent can have
         * @param keepDormantVextexesAndEdges keep vertices that don't have an out connection
         * @param mutationChanceArg percentage of children to be mutated
         */
        void
        initialisePopulation(unsigned int populationSizeArg, unsigned int maxDeepVerticesArg, unsigned int maxEdgesArg,
                             unsigned int edgeTraverseLimitArg, bool keepDormantVerticesAndEdgesArg,
                             double mutationChanceArg);

        /**
         * Calculate the fitness of the population.
         */
        evolution::Metrics calculateFitness(enums::FitnessMetric fitnessMetric, double vertexContribution, double edgeContribution);

        /**
         * Sample the population.
         * @param type selection type
         * @param agentsToKeep number of agents to keep
         * @param keepFittest should the fittest agent be kept
         */
        void sample(enums::SelectionType type, unsigned int agentsToKeep, bool keepFittest);

        /**
         * Crossover the population and mutate the children. Operations are threaded.
         */
        void crossoverAndMutate();

        /**
         * Get the population vector.
         * @return population vector
         */
        [[nodiscard]] std::vector<std::shared_ptr<evolution::Agent>> getPopulation() const;

        /**
         * Get the fittest agent of this population.
         * @return fittest agent
         */
        [[nodiscard]] std::shared_ptr<evolution::Agent> getFittestAgent() const;

        /**
         * Calculate average population fitness.
         * @return average fitness
         */
        [[nodiscard]] double getAverageFitness() const;

        /**
         * Get the testing set.
         * @return testing set
         */
        [[nodiscard]] std::vector<std::shared_ptr<data_structures::DataInstance>> getTestingValues() const;

        /**
         * Get the number of output vertices this population agents have.
         * @return number of outputs
         */
        [[nodiscard]] unsigned int getNumberOfOutputs() const;

        /**
         * Get this population's output labels.
         * @return output labels
         */
        [[nodiscard]] std::vector<std::string> getOutputLabels() const;

        [[nodiscard]] std::vector<std::string> getInputLabels() const;

        /**
         * Return a string representation of this population.
         * @param technical format for https://csacademy.com/app/graph_editor/
         * @return string representation
         */
        std::string toString(bool technical);

        /**
         * Create a new row for the exported csv file.
         * @param delimiter csv delimiter
         * @param iteration population iteration
         * @return csv line
         */
        std::string fitnessToCSVString(char delimiter, unsigned int iteration);

        /**
         * Minimize the agent. Removes stranded deep vertices, deep vertices without (or all recursive) input edges and
         * deep vertices without (or all recursive) output edges.
         * @param agentToMinimize agent
         * @return new minimized agent
         */
        std::shared_ptr<Agent> minimizeAgent(const std::shared_ptr<evolution::Agent> &agentToMinimize);

        /**
         * Calculate average population accuracy.
         * @return average accuracy
         */
        [[nodiscard]] double getAverageAccuracy() const;

        /**
         * Calculate average population mcc.
         * @return average mcc
         */
        [[nodiscard]] double getAverageMatthewsCorrelationCoefficient() const;

        void addAgent(const std::shared_ptr<evolution::Agent> &agent);

    private:
        /**
         * Create a new random agent with a random graph.
         * @return
         */
        std::shared_ptr<evolution::Agent> createAgent();

        /**
         * Adapted from @see https://en.wikipedia.org/wiki/Stochastic_universal_sampling
         * @param agentsToKeep number of agents to keep
         * @param keepFittest should the fittest agent be kept
         */
        std::vector<unsigned int> stochasticUniversalSampling(unsigned int agentsToKeep, bool keepFittest);

        /**
         * Creates new deep vertices for childAgent. More info in the implementation.
         * @param childAgent agent to modify
         * @param parent1 first parent
         * @param parent2 second parent
         * @param checkedIndices vector of checked indices
         */
        void crossoverDeepVertices(const std::shared_ptr<evolution::Agent> &childAgent,
                                   const std::shared_ptr<evolution::Agent> &parent1,
                                   const std::shared_ptr<evolution::Agent> &parent2,
                                   std::vector<unsigned int> &checkedIndices);

        /**
         * Create new edges for childAgent. More info in the implementation.
         * @param childAgent agent to modify
         * @param parent1 first parent
         * @param parent2 second parent
         * @param checkedIndices vector of checked indices
         */
        void crossoverEdges(const std::shared_ptr<evolution::Agent> &childAgent,
                            const std::shared_ptr<evolution::Agent> &parent1,
                            const std::shared_ptr<evolution::Agent> &parent2,
                            std::vector<unsigned int> &checkedIndices);

        /**
         * Add a new random edge to the graph.
         * @param index edge index
         * @param graph graph to modify
         */
        void addRandomEdge(unsigned int index, std::shared_ptr<data_structures::Graph> const &graph);

        /**
         * Calculate agent fitness. It can be based on accuracy or mcc, the number of vertices and edges also contribute to the score.
         * @param fitnessMetric either accuracy or mcc
         * @param vertexContribution vertex contribution, should be <=0
         * @param edgeContribution edge contribution, should be <=0
         * @param agent agent
         */
        void
        calculateAgentFitness(enums::FitnessMetric fitnessMetric, double vertexContribution, double edgeContribution,
                              const std::shared_ptr<evolution::Agent> &agent) const;

        /**
         * Threaded crossover operation. More info in the implementation.
         * @param populationDistribution random distribution of the population to prevent creating a new one
         * @return new child agent
         */
        std::shared_ptr<evolution::Agent>
        crossoverThreaded(std::uniform_int_distribution<unsigned int> populationDistribution);

        /**
         * Threaded mutation operation. More info in the implementation.
         * @param agentIndex child agent to mutate
         */
        void mutateThreaded(unsigned long agentIndex);

        /**
         * Sorting function - by fitness descending.
         * @param a1 first agent
         * @param a2 second agent
         * @return agent order
         */
        static bool
        sortByFitness(const std::shared_ptr<evolution::Agent> &a1, const std::shared_ptr<evolution::Agent> &a2);

        /**
         * Add a new random edge that is attached to deepVertex.
         * @param graph graph to modify
         * @param deepVertex deep vertex to attach to
         * @param isinputEdge does deepVertex act as input or output
         */
        void addRandomEdge(std::shared_ptr<data_structures::Graph> const &graph,
                           std::shared_ptr<data_structures::DeepVertex> const &deepVertex, bool isinputEdge);
    };
}


#endif //NEUROEVOLUTION_POPULATION_H
