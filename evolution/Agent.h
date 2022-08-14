//
// Created by jure on 11. 03. 21.
//

#ifndef NEUROEVOLUTION_AGENT_H
#define NEUROEVOLUTION_AGENT_H

#include "../data_structures/Graph.h"
#include "../data_structures/IDeepCloneable.h"

namespace evolution {
    class Agent : public data_structures::IDeepCloneable<Agent> {
    private:
        std::shared_ptr<data_structures::Graph> graph{};
        double fitness = 0;
        double accuracy = 0;
        double matthewsCorrelationCoefficient = 0;
        bool newAgent = true;

    public:
        /**
         * Construct a new agent with fitness, accuracy and mcc.
         * @param graph graph
         * @param fitness fitness
         * @param accuracy accuracy
         * @param matthewsCorrelationCoefficient matthewsCorrelationCoefficient
         */
        Agent(std::shared_ptr<data_structures::Graph> graph, double fitness, double accuracy,
              double matthewsCorrelationCoefficient) : graph(std::move(graph)), fitness(fitness), accuracy(accuracy),
                                                       matthewsCorrelationCoefficient(matthewsCorrelationCoefficient) {}

        /**
         * Construct a new agent with fitness, accuracy and mcc 0.
         * @param graph graph
         */
        explicit Agent(std::shared_ptr<data_structures::Graph> graph) : graph(std::move(graph)) {}

        /**
         * Create a new agent shared pointer.
         * @param graphArg graph
         * @return shared pointer
         */
        static std::shared_ptr<Agent> create(const std::shared_ptr<data_structures::Graph> &graphArg);

        /**
         * Create an agent with an empty graph.
         * @return a new agent
         */
        /**
         * Create an agent with an empty (unconnected) graph.
         * @param inputVertices number of input vertices
         * @param inputLabels input vertex labels
         * @param outputVertices number of output vertices
         * @param outputLabels output vertex labels
         * @return new agent shared pointer
         */
        static std::shared_ptr<Agent>
        create(unsigned int inputVertices, std::vector<std::string> &inputLabels, unsigned int outputVertices,
               std::vector<std::string> &outputLabels);

        /**
         * Deep clones this object.
         * @return a new instance
         */
        std::shared_ptr<Agent> deepClone() override;

        /**
         * Construct a new agent with fitness 0.
         */
        Agent() = default;

        /**
         * Default destructor.
         */
        virtual ~Agent() = default;

        /**
         * Get the agent's graph.
         * @return graph
         */
        [[nodiscard]] std::shared_ptr<data_structures::Graph> getGraph() const;

        /**
         * Get the agent's fitness.
         * @return fitness
         */
        [[nodiscard]] double getFitness() const;

        /**
         * Set the agent's newFitness.
         * @param newFitness newFitness
         */
        void setFitness(double newFitness);

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        [[nodiscard]] std::string toString(bool technical) const;

        /**
         * Is this agent new? Determines whether traversal happens.
         * @return if this agent is new
         */
        [[nodiscard]] bool isNewAgent() const;

        /**
         * Sets whether this agent is new.
         * @param newAgentArg is this agent new?
         */
        void setNewAgent(bool newAgentArg);

        /**
         * Get the accuracy of this agent.
         * @return accuracy
         */
        [[nodiscard]] double getAccuracy() const;

        /**
         * Sets the accuracyArg of this agent.
         * @param accuracyArg accuracyArg
         */
        void setAccuracy(double accuracyArg);

        /**
         * Gets the mcc of this agent.
         * @return mcc
         */
        [[nodiscard]] double getMatthewsCorrelationCoefficient() const;

        /**
         * Sets the mcc of this agent.
         * @param matthewsCorrelationCoefficientArg mcc
         */
        void setMatthewsCorrelationCoefficient(double matthewsCorrelationCoefficientArg);
    };
}


#endif //NEUROEVOLUTION_AGENT_H
