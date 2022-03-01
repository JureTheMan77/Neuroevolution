//
// Created by jure on 1. 02. 21.
//

#include <fstream>
#include <filesystem>
#include <algorithm>
#include <thread>
#include <future>
#include "Population.h"
#include "../util/util.h"
#include "../logging/logging.h"
#include "../data_structures/MulticlassConfusionMatrix.h"
#include "../enums/FitnessMetric.h"

evolution::Population::Population(const std::string &pathToDataSet) {
    if (!std::filesystem::exists(pathToDataSet)) {
        throw std::invalid_argument("File " + pathToDataSet + " does not exist.");
    }

    std::string delimiter = ",";
    std::vector<std::shared_ptr<data_structures::DataInstance>> fullDataSet;

    // read the file line by line
    unsigned int counter = 0;
    std::ifstream infile(pathToDataSet);
    std::ifstream file(pathToDataSet);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {

            if (counter == 0) {
                // counter = 0 -> input labels
                this->inputLabels = util::split(line, delimiter);
                this->numberOfInputs = this->inputLabels.size();
            } else if (counter == 1) {
                // counter = 1 -> output labels
                this->outputLabels = util::split(line, delimiter);
                this->numberOfOutputs = this->outputLabels.size();
            } else {
                auto dataInstance = data_structures::DataInstance::createDataInstance(
                        util::splitDouble(line, delimiter));
                fullDataSet.push_back(dataInstance);
            }

            counter++;
        }
        file.close();
    }

    // shuffle fullDataSet
    auto rng = std::default_random_engine{this->seeder()};
    std::shuffle(std::begin(fullDataSet), std::end(fullDataSet), rng);

    // split the full dataset into training and testing
    // use a 70/30 split
    auto splitMark = (unsigned int) round((double) fullDataSet.size() * 0.7);
    this->trainingValues.insert(trainingValues.begin(), fullDataSet.begin(), fullDataSet.begin() + splitMark);
    this->testingValues.insert(testingValues.begin(), fullDataSet.begin() + splitMark, fullDataSet.end());

    //this->trainingValues.insert(trainingValues.begin(), fullDataSet.begin(), fullDataSet.end());
    //this->testingValues.insert(testingValues.begin(), fullDataSet.begin(), fullDataSet.end());


    this->population = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationPlaceholder = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationSize = 0;
}

void evolution::Population::initialisePopulation(unsigned int populationSizeArg, unsigned int maxDeepVerticesArg,
                                                 unsigned int maxEdgesArg, unsigned int edgeTraverseLimitArg,
                                                 bool keepDormantVerticesAndEdgesArg, double maxMutationChanceArg) {
    this->population = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationSize = populationSizeArg;
    this->edgeTraverseLimit = edgeTraverseLimitArg;
    this->maxDeepVertices = maxDeepVerticesArg;
    this->maxEdges = maxEdgesArg;
    this->maxMutationChance = maxMutationChanceArg;
    this->keepDormantVerticesAndEdges = keepDormantVerticesAndEdgesArg;
    for (unsigned int i = 0; i < populationSizeArg; i++) {
        auto agent = this->createAgent();
        this->population.push_back(agent);
    }
}

std::shared_ptr<evolution::Agent>
evolution::Population::createAgent() {
    std::mt19937 rng(this->seeder());
    std::uniform_int_distribution<unsigned int> maxDeepVerticesDistribution(0, maxDeepVertices);
    std::uniform_int_distribution<unsigned int> edgesDistribution(1, maxEdges);
    unsigned int numberOfDeepVertices = maxDeepVerticesDistribution(rng);
    unsigned int numberOfEdges = edgesDistribution(rng);

    // generate vertices
    auto graph = data_structures::Graph::createGraph(this->numberOfInputs, this->inputLabels, numberOfDeepVertices,
                                                     this->numberOfOutputs, this->outputLabels,
                                                     this->maxMutationChance);

    // generate edges
    for (unsigned int i = 0; i < numberOfEdges; i++) {
        this->addRandomEdge(i, graph);
    }

    // normalize edge weights
    // find the largest weight
    double maxWeight = 0;
    for (const auto &edge: graph->getEdges()) {
        if (edge->getWeight() > maxWeight) {
            maxWeight = edge->getWeight();
        }
    }
    for (const auto &edge: graph->getEdges()) {
        double newWeight = edge->getWeight() / maxWeight;
        //newWeight = (int) (newWeight * 100000.0) / 100000.0;
        edge->setWeight(newWeight);
    }

    auto agent = evolution::Agent::create(graph);
    //if (!this->keepDormantVerticesAndEdges) {
    //    agent = this->minimizeAgent(agent);
    //}
    return agent;
}

void evolution::Population::addRandomEdge(unsigned int index, std::shared_ptr<data_structures::Graph> const &graph) {
    std::mt19937 rng(this->seeder());
    std::mt19937_64 rng64(this->seeder());
    std::uniform_int_distribution<unsigned int> inputVerticesDistribution(0, graph->getInputVertices().size() - 1);
    std::uniform_int_distribution<unsigned int> vertexTypeDistribution(0, 1);
    std::uniform_int_distribution<unsigned int> deepVerticesDistribution(0, graph->getDeepVertices().size() - 1);
    std::uniform_real_distribution<double> weightDistribution(0, 1);
    std::uniform_int_distribution<unsigned int> edgeTraverseLimitDistribution(1, this->edgeTraverseLimit);
    std::uniform_int_distribution<unsigned int> outputVerticesDistribution(0, graph->getOutputVertices().size() - 1);

    // pick two edges, they can be:
    // 1 input and 1 deep
    // 1 input and 1 output
    // 2 deep
    // 1 deep or 1 output
    unsigned int type1;
    unsigned int type2;

    if (!graph->getDeepVertices().empty()) {
        // 0 - input
        // 1 - deep
        type1 = vertexTypeDistribution(rng);
        if (type1 == 0) {
            type2 = vertexTypeDistribution(rng);
            // 0 - deep
            // 1 - output
            if (type2 == 0) {
                graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(rng), enums::VertexType::Deep,
                               deepVerticesDistribution(rng), index, weightDistribution(rng64),
                               edgeTraverseLimitDistribution(rng), this->maxMutationChance);
            } else {
                graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(rng), enums::VertexType::Output,
                               outputVerticesDistribution(rng), index, weightDistribution(rng64),
                               edgeTraverseLimitDistribution(rng), this->maxMutationChance);
            }
        } else {
            type2 = vertexTypeDistribution(rng);
            // 0 - deep
            // 1 - output
            if (type2 == 0) {
                graph->addEdge(enums::VertexType::Deep, deepVerticesDistribution(rng), enums::VertexType::Deep,
                               deepVerticesDistribution(rng), index, weightDistribution(rng64),
                               edgeTraverseLimitDistribution(rng), this->maxMutationChance);
            } else {
                graph->addEdge(enums::VertexType::Deep, deepVerticesDistribution(rng), enums::VertexType::Output,
                               outputVerticesDistribution(rng), index, weightDistribution(rng64),
                               edgeTraverseLimitDistribution(rng), this->maxMutationChance);
            }
        }
    } else {
        // no deep vertices
        // only legal configuration is input -> output
        graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(rng), enums::VertexType::Output,
                       outputVerticesDistribution(rng), index, weightDistribution(rng64),
                       edgeTraverseLimitDistribution(rng), this->maxMutationChance);
    }
}

std::string evolution::Population::toString(bool technical) {
    std::ostringstream result;
    result << "Number of agents: " << this->population.size() << std::endl;
    for (unsigned int i = 0; i < this->population.size(); i++) {
        result << "Agent " << i << ": " << std::endl;
        result << this->population.at(i)->toString(technical) << std::endl;
    }

    return result.str();
}

void evolution::Population::calculateFitness(enums::FitnessMetric fitnessMetric, double vertexContribution,
                                             double edgeContribution) {
    // run all data instances on all agents
    std::vector<std::future<void>> futures;
    for (const std::shared_ptr<evolution::Agent> &agent: this->population) {
        // calculate fitness only for new agents
        if (agent->isNewAgent()) {
            auto ftr = std::async(&evolution::Population::calculateAgentFitness, this, fitnessMetric,
                                  vertexContribution, edgeContribution, agent);
            futures.push_back(std::move(ftr));
        }
        //calculateAgentFitness(vertexContribution, edgeContribution, std::ref(agent));
    }
    // join the created threads
    for (auto &ftr: futures) {
        ftr.get();
    }

    // find the lowest fitness
    // if it's negative, then assign it's absolute value to the population

}

void evolution::Population::calculateAgentFitness(enums::FitnessMetric fitnessMetric, double vertexContribution,
                                                  double edgeContribution,
                                                  const std::shared_ptr<evolution::Agent> &agent) {
    unsigned int numCorrect = 0;
    unsigned int numIncorrect = 0;

//    for (const std::shared_ptr<data_structures::DataInstance> &di: trainingValues) {
//        agent->getGraph()->traverse(di);
//
//        // check if the prediction is correct
//        unsigned int predictedIndex = agent->getGraph()->getLargestOutputValueIndex();
//        if (predictedIndex == di->getCorrectIndex()) {
//            numCorrect += 1;
//        } else {
//            numIncorrect += 1;
//        }
//
//        // reset the agent
//        agent->getGraph()->reset();
//    }

    data_structures::MulticlassConfusionMatrix mcm(agent, this->trainingValues, this->numberOfOutputs);
    //logging::logs(mcm.toString(this->outputLabels));

    // size of the agent is a penalty, both edges and vertices contribute -0.1
    double sizeContribution = (double) agent->getGraph()->getDeepVertices().size() * vertexContribution +
                              (double) agent->getGraph()->getEdges().size() * edgeContribution;

    // calculate the fitness (accuracy)
// agent->setFitness((double) numCorrect / (double) this->trainingValues.size());
    // agent->setFitness((double) numCorrect + sizeContribution);
    double fitness = sizeContribution;

    switch (fitnessMetric) {
        case enums::FitnessMetric::Accuracy:
            fitness += mcm.getAccuracy();
            break;
        case enums::FitnessMetric::MatthewsCorrelationCoefficient:
            // mcc can be negative, up to -1
            fitness += mcm.getMatthewsCorrelationCoefficient() + 1;
            break;
        default:
            break;
    }

    if (fitness < 0) {
        fitness = 0;
    }
    agent->setFitness(fitness);
}

void evolution::Population::sample(enums::SelectionType type, unsigned int agentsToKeep) {
    std::vector<unsigned int> indexesToKeep{};
    if (this->populationSize < agentsToKeep) {
        throw std::invalid_argument("Number of agents to keep cannot be higher than the population size.");
    }

    if (type == enums::SelectionType::StochasticUniversalSampling) {
        indexesToKeep = stochasticUniversalSampling(agentsToKeep);
    } else {
        throw std::invalid_argument("This selection type is not supported.");
    }

    // copy all agents that were selected
    // the reason for copying is that if we select an agent twice and in this selection and don't select
    // only one reference the next time, both agents will be deleted
    for (unsigned int index: indexesToKeep) {
        // indexesToKeep will never be larger that the population
        auto agent = this->population.at(index)->deepClone();
        agent->setNewAgent(false);
        this->populationPlaceholder.push_back(agent);
    }

    this->population.clear();
    this->population.insert(this->population.begin(), this->populationPlaceholder.begin(),
                            this->populationPlaceholder.end());
    this->populationPlaceholder.clear();
}

std::vector<unsigned int> evolution::Population::stochasticUniversalSampling(unsigned int agentsToKeep) {
    // get highest fitness
    double highestFitness = -std::numeric_limits<double>::max();
    unsigned int fittestAgentIndex = 0;
    for (unsigned int i = 0; i < this->population.size(); i++) {
        if (highestFitness < this->population.at(i)->getFitness()) {
            highestFitness = this->population.at(i)->getFitness();
            fittestAgentIndex = i;
        }
    }

    // calculate the total fitness of the population
    double fitnessSum = 0;
    for (const std::shared_ptr<evolution::Agent> &agent: this->population) {
        fitnessSum += agent->getFitness();
    }

    // distance between selections
    double distance = fitnessSum / agentsToKeep;

    // determine the starting offset
    std::uniform_real_distribution<double> startDistribution(0, distance);
    std::mt19937_64 rng64(this->seeder());
    double start = startDistribution(rng64);

    // always sorted
    // std::sort(pointers.begin(), pointers.end());

    // roulette wheel selection
    // always keep the fittest agent
    std::vector<unsigned int> indexesToKeep{fittestAgentIndex};

    // loop until the vector of agents is filled
    double pointer = start;
    int counter = 0;
    double cumulativeFitness = this->population.at(0)->getFitness();
    while (indexesToKeep.size() < agentsToKeep) {
        if (cumulativeFitness >= pointer) {
            indexesToKeep.push_back(counter);
            pointer += distance;
        } else {
            counter += 1;
            cumulativeFitness += this->population.at(counter)->getFitness();
        }
    }

    // indexes will always be sorted
    // std::sort(indexesToKeep.begin(), indexesToKeep.end());

    return indexesToKeep;
}

void evolution::Population::crossover() {
    // sort the population according to fitness
    // std::sort(this->population.begin(), this->population.end(), sortByFitness);

    // create new agents until the population is full
    // save child agents in populationPlaceholder
    unsigned int agentsToCreate = this->populationSize - this->population.size();
    std::vector<std::future<std::shared_ptr<evolution::Agent>>> futures;

    for (int i = 0; i < agentsToCreate; i++) {
        auto ftr = std::async(&evolution::Population::crossoverThreaded, this);
        futures.push_back(std::move(ftr));
        //this->crossoverThreaded();
    }
    for (auto &ftr: futures) {
        populationPlaceholder.push_back(ftr.get());
    }

    // combine the populationPlaceholder with the population
    this->population.insert(this->population.end(),
                            std::make_move_iterator(this->populationPlaceholder.begin()),
                            std::make_move_iterator(this->populationPlaceholder.end()));

    this->populationPlaceholder.clear();
}

std::shared_ptr<evolution::Agent> evolution::Population::crossoverThreaded() {
    // initialize randomness
    std::mt19937 rng(this->seeder());
    // choose the maximum size here since we don't want to do crossover with child agents 0, (double)this->population.size() - 1)
    std::uniform_int_distribution<unsigned int> populationDistribution(0, this->population.size() - 1);
    // choose 2 random agents
    unsigned long firstAgentIndex = 0;
    unsigned long secondAgentIndex = 0;
    unsigned int firstRoll = populationDistribution(rng);
    unsigned int secondRoll = populationDistribution(rng);

    if (population.at(firstRoll)->getFitness() > population.at(secondRoll)->getFitness()) {
        firstAgentIndex = firstRoll;
    } else {
        firstAgentIndex = secondRoll;
    }
    firstRoll = populationDistribution(rng);
    secondRoll = populationDistribution(rng);
    if (population.at(firstRoll)->getFitness() > population.at(secondRoll)->getFitness()) {
        secondAgentIndex = firstRoll;
    } else {
        secondAgentIndex = secondRoll;
    }

    // create an empty child agent
    std::shared_ptr<evolution::Agent> childAgent = Agent::create(numberOfInputs, inputLabels,
                                                                 numberOfOutputs,
                                                                 outputLabels, this->maxMutationChance);

    std::shared_ptr<evolution::Agent> agent1 = population.at(firstAgentIndex);
    std::shared_ptr<evolution::Agent> agent2 = population.at(secondAgentIndex);
    //logging::logs(std::to_string(agent1->getGraph()->getEdges().at(0)->getMutationChance()));
    //logging::logs(std::to_string(agent2->getGraph()->getEdges().at(0)->getMutationChance()));

    // output and input vertices must be fixed, so only do crossover on deep vertices and their edges

    // 1. vertex crossover
    std::vector<unsigned int> checkedIndices;
    // - iterate over the deep vertices of the first agent
    crossoverDeepVertices(childAgent, agent1, agent2, checkedIndices);
    // - iterate over the deep vertices of the second agent
    crossoverDeepVertices(childAgent, agent2, agent1, checkedIndices);

    // 2. edge crossover on the populationPlaceholder vector
    checkedIndices.clear();
    // - iterate over the edges of the first agent
    crossoverEdges(childAgent, agent1, agent2, checkedIndices);
    // - iterate over the edges of the second agent
    crossoverEdges(childAgent, agent2, agent1, checkedIndices);
    // add random edges based on the number of vertices with index UINT_MAX
    if (childAgent->getGraph()->getEdges().size() < this->maxEdges) {
        addNewRandomEdges(childAgent);
    }

    // fix indices for vertices that are "twins" and randomly added edges
    childAgent->getGraph()->fixIndices();

    // populationPlaceholder.push_back(childAgent);

    // normalize edge weights
    // find the largest weight
    double maxWeight = 0;
    for (const auto &edge: childAgent->getGraph()->getEdges()) {
        if (edge->getWeight() > maxWeight) {
            maxWeight = edge->getWeight();
        }
    }
    for (const auto &edge: childAgent->getGraph()->getEdges()) {
        double newWeight = edge->getWeight() / maxWeight;
        //newWeight = (int) (newWeight * 100000.0) / 100000.0;
        edge->setWeight(newWeight);
    }

    if (!this->keepDormantVerticesAndEdges) {
        childAgent = this->minimizeAgent(childAgent);
    }
    return childAgent;
}

void
evolution::Population::crossoverDeepVertices(const std::shared_ptr<evolution::Agent> &childAgent,
                                             const std::shared_ptr<evolution::Agent> &parent1,
                                             const std::shared_ptr<evolution::Agent> &parent2,
                                             std::vector<unsigned int> &checkedIndices) {
    // loop through all vertices
    for (const auto &deepVertex1: parent1->getGraph()->getDeepVertices()) {

        if (childAgent->getGraph()->getDeepVertices().size() >= this->maxDeepVertices) {
            break;
        }

        // check if this vertex was already crossovered
        bool matchFound = false;
        for (auto index: checkedIndices) {
            if (deepVertex1->getIndex() == index) {
                matchFound = true;
                break;
            }
        }
        if (matchFound) {
            continue;
        }
        // get the vertex with the same index from parent2
        auto deepVertex2 = parent2->getGraph()->getDeepVertexByIndex(deepVertex1->getIndex());
        // if no vertex was found, then do the following:
        // - if the vertex is dominant, then always keep it
        // - if it's not dominant, then keep it with a chance of 1-deepVertex1->getChanceToGetDominated()
        std::vector<std::shared_ptr<data_structures::DeepVertex>> childVertices;
        if (deepVertex2 == nullptr) {
            childVertices = createDeepVertexChildren(deepVertex1);
        } else {
            // if vertex2 was found, then choose a dominant vertex
            auto dominant = this->getDominantVertex(deepVertex1, deepVertex2);
            childVertices = createDeepVertexChildren(dominant);
        }

        checkedIndices.push_back(deepVertex1->getIndex());

        childAgent->getGraph()->addDeepVertices(childVertices);
    }
}

std::shared_ptr<data_structures::DeepVertex>
evolution::Population::getDominantVertex(const std::shared_ptr<data_structures::DeepVertex> &v1,
                                         const std::shared_ptr<data_structures::DeepVertex> &v2) {
    unsigned int result = getDominantCrossoverable(v1, v2);
    return result == 1 ? v1 : v2;
}

std::shared_ptr<data_structures::Edge>
evolution::Population::getDominantEdge(const std::shared_ptr<data_structures::Edge> &e1,
                                       const std::shared_ptr<data_structures::Edge> &e2) {
    unsigned int result = getDominantCrossoverable(e1, e2);
    return result == 1 ? e1 : e2;
}

unsigned int evolution::Population::getDominantCrossoverable(const std::shared_ptr<data_structures::ICrossoverable> &c1,
                                                             const std::shared_ptr<data_structures::ICrossoverable> &c2) {
    double randomRoll = util::nextDouble();
    if (c1->isDominant() && c2->isDominant() || !c1->isDominant() && !c2->isDominant()) {
        // both vertices aren't dominant or recessive
        if (c1->getChanceToGetDominated() < c2->getChanceToGetDominated()) {
            // v1 has the least chance to get dominated
            return c1->getChanceToGetDominated() > randomRoll ? 2 : 1;
        } else {
            // v2 has the least chance to get dominated
            return c2->getChanceToGetDominated() > randomRoll ? 1 : 2;
        }
    } else if (c1->isDominant()) {
        // only v1 is dominant
        return c1->getChanceToGetDominated() > randomRoll ? 2 : 1;
    } else if (c2->isDominant()) {
        // only v2 is dominant, however there's a chance it's not
        return c2->getChanceToGetDominated() > randomRoll ? 1 : 2;
    } else {
        throw std::invalid_argument("Unable to determine the dominant crossoverable object.");
    }
}

std::vector<std::shared_ptr<data_structures::DeepVertex>>
evolution::Population::createDeepVertexChildren(const std::shared_ptr<data_structures::DeepVertex> &vertex) {
    // get the number of children to produce
    unsigned int childrenToProduce = this->childrenToProduce(vertex);

    std::vector<std::shared_ptr<data_structures::DeepVertex>> childDeepVertices;

    if (childrenToProduce == 0) {
        // this crossover operation produced no children
        return childDeepVertices;
    }

    // create children
    for (unsigned int i = 0; i < childrenToProduce; i++) {
        // check if the child has the max deep vertices allowed
        if (childDeepVertices.size() >= this->maxDeepVertices) {
            break;
        }

        auto childData = this->mutateCrossoverable(vertex);

        // create the child vertex
        // if this is not the first iteration, then assign UINT_MAX and fix it later
        childDeepVertices.emplace_back(
                data_structures::DeepVertex::createDeepVertex(i == 0 ? vertex->getIndex() : UINT_MAX,
                                                              childData.isDominant(),
                                                              childData.getChanceToGetDominated(),
                                                              childData.getMutationChance(),
                                                              childData.getMaxChildren()));
    }

    return childDeepVertices;
}

void evolution::Population::crossoverEdges(std::shared_ptr<evolution::Agent> &childAgent,
                                           const std::shared_ptr<evolution::Agent> &parent1,
                                           const std::shared_ptr<evolution::Agent> &parent2,
                                           std::vector<unsigned int> &checkedIndices) {
    // iterate over the edges of the first parent
    for (const auto &edge1: parent1->getGraph()->getEdges()) {
        if (childAgent->getGraph()->getEdges().size() >= this->maxEdges) {
            break;
        }

        // check if this edge1 was already crossovered
        bool matchFound = false;
        for (unsigned int index: checkedIndices) {
            if (edge1->getIndex() == index) {
                matchFound = true;
                break;
            }
        }
        if (matchFound) {
            continue;
        }
        // get the edge1 with the same input and output indices (and types) from the second parent
        auto edge2 = parent2->getGraph()->getEdgeByIndexAndType(edge1->getInput()->getIndex(),
                                                                edge1->getInput()->getType(),
                                                                edge1->getOutput()->getIndex(),
                                                                edge1->getOutput()->getType());
        std::vector<std::shared_ptr<data_structures::Edge>> childEdges;
        if (edge2 == nullptr) {
            // similar edge1 was not found in parent2
            // create children from just edge1
            this->createAndAddEdgeChildren(edge1, childAgent);
        } else {
            // a similar edge1 was found, then choose the dominant one
            auto dominant = this->getDominantEdge(edge1, edge2);
            this->createAndAddEdgeChildren(dominant, childAgent);
        }
        checkedIndices.push_back(edge1->getIndex());
    }
}

void evolution::Population::createAndAddEdgeChildren(const std::shared_ptr<data_structures::Edge> &edge,
                                                     std::shared_ptr<evolution::Agent> &childAgent) {
    std::uniform_real_distribution<double> distDouble(0, 1);
    // the maximum weight is the current weight * 2
    double maxWeight = edge->getWeight() < 1 ? 1 : edge->getWeight() * 2;
    std::uniform_real_distribution<double> distWeight(0, maxWeight);
    std::mt19937_64 rng64(seeder());

    // get the number of children to produce
    unsigned int childrenToProduce = this->childrenToProduce(edge);
    if (childrenToProduce == 0) {
        // this crossover operation produced no children
        return;
    }

    // create children
    auto childData = this->mutateCrossoverable(edge);
    // mutate weight
    double mutationRoll = distDouble(rng64);

    // calculate the weight of the first child
    double childWeight = edge->getMutationChance() > mutationRoll ? distWeight(rng64) : edge->getWeight();

    // mutate traversal limit
    mutationRoll = distDouble(rng64);
    unsigned int childTraversalLimit =
            edge->getMutationChance() > mutationRoll ? util::nextUnsignedInt(1, this->edgeTraverseLimit)
                                                     : edge->getTraverseLimit();

    // add the new edge
    auto childEdge = childAgent->getGraph()->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(),
                                                     edge->getOutput()->getType(), edge->getOutput()->getIndex(),
                                                     edge->getIndex(), childWeight, childTraversalLimit, childData);
    // end the procedure if no edge was created
    if (childEdge == nullptr) {
        return;
    }
    // combine weights if multiple children must be created
    for (int i = 1; i < childrenToProduce; i++) {
        mutationRoll = distDouble(rng64);
        childWeight = edge->getMutationChance() > mutationRoll ? distWeight(rng64) : edge->getWeight();
        childEdge->combineWeight(childWeight);
    }
}

data_structures::ICrossoverable
evolution::Population::mutateCrossoverable(const std::shared_ptr<data_structures::ICrossoverable> &crossoverable) {
    std::uniform_real_distribution<double> distDouble(0, 1);
    std::uniform_real_distribution<double> distNewMutationChance(0, this->maxMutationChance);
    std::mt19937_64 rng64(seeder());

    double mutationRoll = distDouble(rng64);
    // take the values of the dominant vertex with a possible mutation
    // dominant
    bool childDominant = crossoverable->getMutationChance() > mutationRoll ? !crossoverable->isDominant()
                                                                           : crossoverable->isDominant();
    // chanceToGetDominated
    mutationRoll = distDouble(rng64);
    double childChanceToGetDominated = crossoverable->getMutationChance() > mutationRoll ? distDouble(rng64)
                                                                                         : crossoverable->getChanceToGetDominated();
    // mutationChance
    mutationRoll = distDouble(rng64);
    double childMutationChance =
            crossoverable->getMutationChance() > mutationRoll ? distNewMutationChance(rng64)
                                                              : crossoverable->getMutationChance();
    // maxChildren
    mutationRoll = distDouble(rng64);
    unsigned int childMaxChildren = crossoverable->getMutationChance() > mutationRoll ? util::nextUnsignedInt(1, 2)
                                                                                      : crossoverable->getMaxChildren();

    return data_structures::ICrossoverable(childMutationChance, childDominant, childChanceToGetDominated,
                                           childMaxChildren);
}

unsigned int
evolution::Population::childrenToProduce(const std::shared_ptr<data_structures::ICrossoverable> &crossoverable) {
    unsigned int childrenToProduce = crossoverable->getMaxChildren();
    // choose whether this vertex produces one less child (either 0 or 1)
    double mutationRoll = util::nextDouble();
    childrenToProduce = crossoverable->getMutationChance() > mutationRoll ? childrenToProduce - 1 : childrenToProduce;
    return childrenToProduce;
}

const std::vector<std::shared_ptr<evolution::Agent>> &evolution::Population::getPopulation() const {
    return this->population;
}

std::shared_ptr<evolution::Agent> evolution::Population::getFittestAgent() {
    std::shared_ptr<evolution::Agent> fittestAgent = this->population.at(0);
    for (auto const &agent: this->population) {
        if (fittestAgent->getFitness() < agent->getFitness()) {
            fittestAgent = agent;
        }
    }
    return fittestAgent;
}

void evolution::Population::addAgent(std::shared_ptr<evolution::Agent> agent) {
    this->population.push_back(agent);
}

double evolution::Population::getAverageFitness() {
    double average = 0;
    for (const auto &agent: this->population) {
        average += agent->getFitness();
    }
    return average / (double) this->population.size();
}

void evolution::Population::addNewRandomEdges(std::shared_ptr<evolution::Agent> const &childAgent) {
    // add random edges equal to the amount of vertices with index value of UINT_MAX
    for (auto const &vertex: childAgent->getGraph()->getDeepVertices()) {
        // check of the child already has the maximum amount of edges allowed
        if (childAgent->getGraph()->getEdges().size() >= this->maxEdges) {
            break;
        }
        if (vertex->getIndex() == UINT_MAX) {
            this->addRandomEdge(UINT_MAX, childAgent->getGraph());
        }
    }

}

std::vector<std::shared_ptr<data_structures::DataInstance>> evolution::Population::getTestingValues() {
    return this->testingValues;
}

unsigned int evolution::Population::getNumberOfOutputs() {
    return this->numberOfOutputs;
}

std::vector<std::string> evolution::Population::getOutputLabels() {
    return this->outputLabels;
}

std::string evolution::Population::fitnessToCSVString(char delimiter, unsigned int iteration) {
    std::string averageFitnessStr = std::to_string(this->getAverageFitness());
    std::replace(averageFitnessStr.begin(), averageFitnessStr.end(), '.', ',');

    std::string fittestAgentFitnessStr = std::to_string(this->getFittestAgent()->getFitness());
    std::replace(fittestAgentFitnessStr.begin(), fittestAgentFitnessStr.end(), '.', ',');

    return std::to_string(iteration) + delimiter + averageFitnessStr + delimiter + fittestAgentFitnessStr;
}

bool evolution::Population::sortByFitness(const std::shared_ptr<evolution::Agent> &a1,
                                          const std::shared_ptr<evolution::Agent> &a2) {
    return a1->getFitness() > a2->getFitness();
}

std::shared_ptr<evolution::Agent> evolution::Population::minimizeFittestAgent() {
    // remove all deep vertices that don't have output edges
    // also remove their input edges
    auto fittestAgent = this->getFittestAgent();
    return minimizeAgent(fittestAgent);

}

std::shared_ptr<evolution::Agent>
evolution::Population::minimizeAgent(const std::shared_ptr<evolution::Agent> &agentToMinimize) {
    // Step 1. if the deep vertex has no input or output edges, then flag it for deletion
    bool keepRunning = true;
    while (keepRunning) {
        keepRunning = false;
        for (const auto &deepVertex: agentToMinimize->getGraph()->getDeepVertices()) {
            if (deepVertex->isFlaggedForDeletion()) {
                continue;
            }
            if (deepVertex->getOutputEdges().empty() && deepVertex->getInputEdges().empty()) {
                // remove all isolated edges
                deepVertex->setFlaggedForDeletion(true);
            } else if (deepVertex->getOutputEdges().empty() || deepVertex->allOutputEdgesFlaggedForDeletion()) {
                deepVertex->setFlaggedForDeletion(true);
                keepRunning = true;
                if (deepVertex->allInputEdgesFlaggedForDeletion()) {
                    continue;
                }
                for (const auto &edge: deepVertex->getInputEdges()) {
                    edge->setFlaggedForDeletion(true);
                }
            } else if (deepVertex->getInputEdges().empty() || deepVertex->allInputEdgesFlaggedForDeletion()) {
                deepVertex->setFlaggedForDeletion(true);
                keepRunning = true;
                if (deepVertex->allOutputEdgesFlaggedForDeletion()) {
                    continue;
                }
                for (const auto &edge: deepVertex->getOutputEdges()) {
                    edge->setFlaggedForDeletion(true);
                }
            }
        }
    }

    // construct the agent
    std::shared_ptr<Agent> minimizedAgent = Agent::create(numberOfInputs, inputLabels,
                                                          numberOfOutputs,
                                                          outputLabels, maxMutationChance);

    // add the deep vertices
    for (const auto &deepVertex: agentToMinimize->getGraph()->getDeepVertices()) {
        if (!deepVertex->isFlaggedForDeletion()) {
            auto vertexClone = deepVertex->deepClone();
            minimizedAgent->getGraph()->addDeepVertex(vertexClone);
        }
    }
    // add the edges
    for (const auto &edge: agentToMinimize->getGraph()->getEdges()) {
        if (!edge->isFlaggedForDeletion()) {
            minimizedAgent->getGraph()->addEdge(edge);
        }
    }

    // normalize weights from 0 to 1
    // get the largest weight
    double maxWeight = 0;
    for (const auto &edge: minimizedAgent->getGraph()->getEdges()) {
        if (edge->getWeight() > maxWeight) {
            maxWeight = edge->getWeight();
        }
    }

    // re-index deep vertices
    unsigned int newIndex = 0;
    for (const auto &deepVertex: minimizedAgent->getGraph()->getDeepVertices()) {
        deepVertex->setIndex(newIndex);
        newIndex += 1;
    }

    // copy over the fitness
    minimizedAgent->setFitness(agentToMinimize->getFitness());


    return minimizedAgent;
}

std::vector<std::string> evolution::Population::getInputLabels() {
    return this->inputLabels;
}
