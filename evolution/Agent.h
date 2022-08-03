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
        bool newAgent = true;

    public:
        /**
         * Construct a new agent with fitness.
         * @param graph graph
         * @param fitness fitness
         */
        Agent(std::shared_ptr<data_structures::Graph> graph, double fitness) :
                graph(std::move(graph)), fitness(fitness) {}

        Agent(Agent const &agent) = default;

        /**
         * Construct a new agent with fitness 0.
         * @param graph graph
         */
        explicit Agent(std::shared_ptr<data_structures::Graph> graph) : graph(std::move(graph)) {}

        static std::shared_ptr<Agent>
        create(const std::shared_ptr<data_structures::Graph> &graphArg, double fitnessArg);

        static std::shared_ptr<Agent> create(const std::shared_ptr<data_structures::Graph> &graphArg);

        /**
         * Create an agent with an empty graph.
         * @return a new agent
         */
        static std::shared_ptr<Agent>
        create(unsigned int inputVertices, std::vector<std::string> &inputLabels, unsigned int outputVertices,
               std::vector<std::string> &outputLabels,double maxMutationChance);

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
         * Also deletes graph.
         */
        ~Agent() = default;

        /**
         * Get the agent's graph.
         * @return graph
         */
        std::shared_ptr<data_structures::Graph> getGraph();

        /**
         * Get the agent's fitness.
         * @return fitness
         */
        double getFitness();

        /**
         * Set the agent's newFitness.
         * @param newFitness newFitness
         */
        void setFitness(double newFitness);

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        std::string toString(bool technical);

        bool isNewAgent() const;

        void setNewAgent(bool newAgentArg);
    };
}


#endif //NEUROEVOLUTION_AGENT_H
