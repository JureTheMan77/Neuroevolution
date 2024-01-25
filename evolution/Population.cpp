//
// Created by jure on 1. 02. 21.
//

#include <fstream>
#include <filesystem>
#include <algorithm>
#include <thread>
#include <future>
#include <sstream>
#include "Population.h"
#include "../util/util.h"
#include "../data_structures/MulticlassConfusionMatrix.h"
#include "Metrics.h"

evolution::Population::Population(const std::string &pathToDataSet) {
    if (!std::filesystem::exists(pathToDataSet)) {
        throw std::invalid_argument("File " + pathToDataSet + " does not exist.");
    }

    std::string delimiter = ",";
    std::vector<std::shared_ptr<data_structures::DataInstance>> fullDataSet;

    // read the file line by line
    unsigned int counter = 0;
    std::ifstream infile(pathToDataSet);
    if (std::ifstream file(pathToDataSet); file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {

            if (counter == 0) {
                // counter = 0 -> input labels
                this->inputLabels = util::split(line, delimiter);
                this->numberOfInputs = (unsigned int) this->inputLabels.size();
            } else if (counter == 1) {
                // counter = 1 -> output labels
                this->outputLabels = util::split(line, delimiter);
                this->numberOfOutputs = (unsigned int) this->outputLabels.size();
            } else {
                auto dataInstance = data_structures::DataInstance::createDataInstance(
                        util::splitDouble(line, delimiter));
                fullDataSet.push_back(dataInstance);
            }

            counter++;
        }
        file.close();
    }

    // split fullDataSet into subsets according to class
    std::vector<std::vector<std::shared_ptr<data_structures::DataInstance>>> splitDataSet(this->numberOfOutputs);
    for (std::shared_ptr<data_structures::DataInstance> &instance: fullDataSet) {
        splitDataSet.at(instance->getCorrectIndex()).push_back(instance);
    }

    // shuffle splitDataSet
    auto rng = std::default_random_engine{this->seeder()};
    std::shuffle(std::begin(fullDataSet), std::end(fullDataSet), rng);
    for (std::vector<std::shared_ptr<data_structures::DataInstance>> &subset: splitDataSet) {
        std::shuffle(std::begin(subset), std::end(subset), rng);
        // split the full dataset into training and testing
        // use a 70/30 split
        auto splitMark = (unsigned int) round((double) subset.size() * 0.7);
        this->trainingValues.insert(trainingValues.end(), subset.begin(), subset.begin() + splitMark);
        this->testingValues.insert(testingValues.end(), subset.begin() + splitMark, subset.end());
    }

    //this->trainingValues.insert(trainingValues.begin(), fullDataSet.begin(), fullDataSet.end());
    //this->testingValues.insert(testingValues.begin(), fullDataSet.begin(), fullDataSet.end());

    this->population = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationPlaceholder = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationSize = 0;
    this->inputVerticesDistribution = std::uniform_int_distribution<unsigned int>(0, numberOfInputs - 1);
    this->outputVerticesDistribution = std::uniform_int_distribution<unsigned int>(0, numberOfOutputs - 1);
}

void evolution::Population::initialisePopulation(unsigned int populationSizeArg, unsigned int maxDeepVerticesArg,
                                                 unsigned int maxEdgesArg, unsigned int edgeTraverseLimitArg,
                                                 bool keepDormantVerticesAndEdgesArg, double mutationChanceArg) {
    this->population = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationSize = populationSizeArg;
    this->edgeTraverseLimit = edgeTraverseLimitArg;
    this->maxDeepVertices = maxDeepVerticesArg;
    this->maxEdges = maxEdgesArg;
    this->maxMutationChance = mutationChanceArg;
    this->keepDormantVerticesAndEdges = keepDormantVerticesAndEdgesArg;
    this->edgeTraverseLimitDistribution = std::uniform_int_distribution<unsigned int>(1, edgeTraverseLimitArg);
    while (this->population.size() < populationSizeArg) {
        auto agent = this->createAgent();
        if(agent != nullptr) {
            this->population.push_back(agent);
        }
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
                                                     this->numberOfOutputs, this->outputLabels);

    // generate edges
    unsigned int edgeIndex = 0;
    while (graph->getEdges().size() < numberOfEdges && graph->getEdges().size() < graph->getNumEdgesPossible()) {
        this->addRandomEdge(graph->getEdges().size(), graph);
    }

    for(const auto &vertex : graph->getOutputVertices()){
        if(vertex->getInputEdges().empty()){
            return nullptr;
        }
    }

    // normalize edge weights
    graph->normalizeEdgeWeights();

    auto agent = evolution::Agent::create(graph);
    return agent;
}

void evolution::Population::addRandomEdge(std::shared_ptr<data_structures::Graph> const &graph,
                                          std::shared_ptr<data_structures::DeepVertex> const &deepVertex,
                                          bool isinputEdge) {
    // if "type" is true, then an input edge will be created
    // otherwise, an output edge will be created

    // choose an input vertex type: either input (true) or deep (false)
    bool choice = util::nextBool();

    // randomize weight
    double weight = util::nextWeight();

    // randomize traverse limit
    unsigned int traverseLimit = this->edgeTraverseLimitDistribution(util::rng);

    if (isinputEdge) {
        enums::VertexType inputVertexType = choice ? enums::VertexType::Input : enums::VertexType::Deep;
        // choose a vertex to serve as input
        std::shared_ptr<data_structures::Vertex> inputVertex;
        if (inputVertexType == enums::VertexType::Input) {
            unsigned int position = this->inputVerticesDistribution(util::rng);
            inputVertex = graph->getInputVertices().at(position);
        } else {
            unsigned long position = util::nextUnsignedLong(0, graph->getDeepVertices().size() - 1);
            inputVertex = graph->getDeepVertices().at(position);
        }
        // create and add edge
        graph->addEdge(inputVertexType, inputVertex->getIndex(), enums::VertexType::Deep, deepVertex->getIndex(),
                       std::numeric_limits<unsigned int>::max(), weight, traverseLimit);
    } else {
        enums::VertexType outputVertexType = choice ? enums::VertexType::Deep : enums::VertexType::Output;
        // choose a vertex to serve as output
        std::shared_ptr<data_structures::Vertex> outputVertex;
        if (outputVertexType == enums::VertexType::Deep) {
            unsigned long position = util::nextUnsignedLong(0, graph->getDeepVertices().size() - 1);
            outputVertex = graph->getDeepVertices().at(position);
        } else {
            unsigned long position = this->outputVerticesDistribution(util::rng);
            outputVertex = graph->getOutputVertices().at(position);
        }
        // create and add edge
        graph->addEdge(enums::VertexType::Deep, deepVertex->getIndex(), outputVertexType,
                       outputVertex->getIndex(), std::numeric_limits<unsigned int>::max(), weight, traverseLimit);
    }
}

void evolution::Population::addRandomEdge(unsigned int index, std::shared_ptr<data_structures::Graph> const &graph) {
    std::uniform_int_distribution<unsigned int> vertexTypeDistribution(0, 1);
    std::uniform_int_distribution<unsigned int> deepVerticesDistribution(0,
                                                                         (unsigned int) graph->getDeepVertices().size() -
                                                                         1);

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
        type1 = vertexTypeDistribution(util::rng);
        if (type1 == 0) {
            type2 = vertexTypeDistribution(util::rng);
            // 0 - deep
            // 1 - output
            if (type2 == 0) {
                unsigned int deepVertexIndex = graph->getDeepVertices().at(deepVerticesDistribution(util::rng))->getIndex();
                graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(util::rng), enums::VertexType::Deep,
                               deepVertexIndex, index, util::nextWeight(),
                               edgeTraverseLimitDistribution(util::rng));
            } else {
                graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(util::rng),
                               enums::VertexType::Output,
                               outputVerticesDistribution(util::rng), index, util::nextWeight(),
                               edgeTraverseLimitDistribution(util::rng));
            }
        } else {
            type2 = vertexTypeDistribution(util::rng);
            // 0 - deep
            // 1 - output
            if (type2 == 0) {
                unsigned int deepVertexIndexIn = graph->getDeepVertices().at(deepVerticesDistribution(util::rng))->getIndex();
                unsigned int deepVertexIndexOut = graph->getDeepVertices().at(deepVerticesDistribution(util::rng))->getIndex();
                graph->addEdge(enums::VertexType::Deep, deepVertexIndexIn, enums::VertexType::Deep,
                               deepVertexIndexOut, index, util::nextWeight(),
                               edgeTraverseLimitDistribution(util::rng));
            } else {
                unsigned int deepVertexIndex = graph->getDeepVertices().at(deepVerticesDistribution(util::rng))->getIndex();
                graph->addEdge(enums::VertexType::Deep, deepVertexIndex, enums::VertexType::Output,
                               outputVerticesDistribution(util::rng), index, util::nextWeight(),
                               edgeTraverseLimitDistribution(util::rng));
            }
        }
    } else {
        // no deep vertices
        // only legal configuration is input -> output
        graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(util::rng), enums::VertexType::Output,
                       outputVerticesDistribution(util::rng), index, util::nextWeight(),
                       edgeTraverseLimitDistribution(util::rng));
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

evolution::Metrics
evolution::Population::calculateFitness(enums::FitnessMetric fitnessMetric, double complexityContribution) {
    // run all data instances on all agents
    std::vector<std::future<void>> futures;
    for (const std::shared_ptr<evolution::Agent> &agent: this->population) {
        // calculate fitness only for new agents
        if (agent->isNewAgent()) {
            auto ftr = std::async(&evolution::Population::calculateAgentFitness, this, fitnessMetric,
                                  complexityContribution, agent);
            futures.push_back(std::move(ftr));
        }
    }
    // join the created threads
    for (auto &ftr: futures) {
        ftr.get();
    }
    // cache the fittest agent
    this->fittestAgent = this->population.at(0);
    for (auto const &agent: this->population) {
        if (this->fittestAgent->getFitness() < agent->getFitness()) {
            this->fittestAgent = agent;
        }
    }

    // calculate average fitness
    double average = 0;
    evolution::Metrics m{};
    for (const auto &agent: this->population) {
        average += agent->getFitness();
        m.addFitness(agent->getFitness());
    }
    this->averageFittness = average / (double) this->population.size();

    // calculate average accuracy
    average = 0;
    for (const auto &agent: this->population) {
        average += agent->getAccuracy();
        m.addAccuracy(agent->getAccuracy());
    }
    this->averageAccuracy = average / (double) this->population.size();

    // calculate average mcc
    average = 0;
    for (const auto &agent: this->population) {
        average += agent->getMatthewsCorrelationCoefficient();
        m.addMcc(agent->getMatthewsCorrelationCoefficient());
    }
    this->averageMatthewsCorrelationCoefficient = average / (double) this->population.size();
    m.lock();
    return m;
}

void evolution::Population::calculateAgentFitness(enums::FitnessMetric fitnessMetric, double complexityContribution,
                                                  const std::shared_ptr<evolution::Agent> &agent) const {
    data_structures::MulticlassConfusionMatrix mcm(agent, this->trainingValues, this->numberOfOutputs);

    // size of the agent can be a penalty, both edges and vertices contribute
    //double sizeContribution = (double) agent->getGraph()->getDeepVertices().size() * vertexContribution +
    //                          (double) agent->getGraph()->getEdges().size() * edgeContribution;

    double sizeContribution = (double) agent->getGraph()->getNumOfPropagations() * complexityContribution;

    // calculate the fitness
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

    // all output vertices have at least one input edge
    bool noInput = false;
    for(const auto &vertex : agent->getGraph()->getOutputVertices()){
        if(vertex->getInputEdges().empty()){
            noInput = true;
            break;
        }
    }

    //for(const auto &vertex : agent->getGraph()->getInputVertices()){
    //    if(vertex->getOutputEdges().empty()){
    //        noInput = true;
    //        break;
    //    }
    //}

    if(noInput) {
        fitness = 0;
    }

    if (fitness < 0.001) {
        fitness = 0.001;
    }
    agent->setFitness(fitness);
    agent->setAccuracy(mcm.getAccuracy());
    agent->setMatthewsCorrelationCoefficient(mcm.getMatthewsCorrelationCoefficient());
}

void evolution::Population::sample(enums::SelectionType type, unsigned int agentsToKeep, bool keepFittest,
                                   enums::FitnessMetric fitnessMetric) {
    std::vector<unsigned int> indexesToKeep{};
    if (this->populationSize < agentsToKeep) {
        throw std::invalid_argument("Number of agents to keep cannot be higher than the population size.");
    }

    if (type == enums::SelectionType::StochasticUniversalSampling) {
        indexesToKeep = stochasticUniversalSampling(agentsToKeep, keepFittest, fitnessMetric);
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

std::vector<unsigned int>
evolution::Population::stochasticUniversalSampling(unsigned int agentsToKeep, bool keepFittest,
                                                   enums::FitnessMetric fitnessMetric) {
    // get highest fitness
    unsigned int fittestAgentIndex = 0;
    double highestFitness = this->population.at(fittestAgentIndex)->getFitness();
    for (unsigned int i = 1; i < this->population.size(); i++) {
        if (highestFitness < this->population.at(i)->getFitness()) {
            highestFitness = this->population.at(i)->getFitness();
            fittestAgentIndex = i;
        } else if (highestFitness == this->population.at(i)->getFitness()) {
            switch (fitnessMetric) {
                case enums::FitnessMetric::Accuracy:
                    if (this->population.at(fittestAgentIndex)->getAccuracy() < this->population.at(i)->getAccuracy()) {
                        fittestAgentIndex = i;
                    }
                    break;
                case enums::FitnessMetric::MatthewsCorrelationCoefficient:
                    if (this->population.at(fittestAgentIndex)->getMatthewsCorrelationCoefficient() <
                        this->population.at(i)->getMatthewsCorrelationCoefficient()) {
                        fittestAgentIndex = i;
                    }
                    break;
            }
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

    // roulette wheel selection
    std::vector<unsigned int> indexesToKeep;
    if (keepFittest) {
        indexesToKeep.push_back(fittestAgentIndex);
    }

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

    return indexesToKeep;
}

void evolution::Population::crossoverAndMutate() {
    // create new agents until the population is full
    // save child agents in populationPlaceholder
    unsigned int agentsToCreate = this->populationSize - (unsigned int) this->population.size();
    std::vector<std::future<std::shared_ptr<evolution::Agent>>> crossoverFutures;
    // create the population distribution here since it's the same for all threads
    std::uniform_int_distribution<unsigned int> populationDistribution(0, (unsigned int) this->population.size() - 1);

    for (int i = 0; i < agentsToCreate; i++) {
        auto ftr = std::async(&evolution::Population::crossoverThreaded, this, populationDistribution);
        crossoverFutures.push_back(std::move(ftr));
    }
    for (auto &ftr: crossoverFutures) {
        this->populationPlaceholder.push_back(ftr.get());
    }

    // mutate children
    // mutate the amount of children proportional to the global mutation chance
    std::vector<std::future<void>> mutationFutures;
    auto childrenToMutate = (unsigned int) std::round(agentsToCreate * this->maxMutationChance);
    std::uniform_int_distribution<unsigned int> childPopulationDistribution(0, agentsToCreate - 1);
    // choose unique agent indexes to mutate
    std::vector<unsigned long> indexesToMutate;
    while (indexesToMutate.size() < childrenToMutate) {
        unsigned long index = childPopulationDistribution(util::rng);
        bool sameIndexFound = false;
        for (const auto &generatedIndex: indexesToMutate) {
            if (generatedIndex == index) {
                sameIndexFound = true;
                break;
            }
        }
        if (!sameIndexFound) {
            indexesToMutate.push_back(index);
        }
    }

    for (const auto &index: indexesToMutate) {
        auto ftr = std::async(&evolution::Population::mutateThreaded, this, index);
        mutationFutures.push_back(std::move(ftr));
    }
    for (auto &ftr: mutationFutures) {
        ftr.get();
    }

    // normalize child agent weights
    std::vector<std::future<void>> normalizationFutures;
    for (const auto &agent: this->populationPlaceholder) {
        auto ftr = std::async(&data_structures::Graph::normalizeEdgeWeights, agent->getGraph());
        normalizationFutures.push_back(std::move(ftr));
    }
    for (auto &ftr: normalizationFutures) {
        ftr.get();
    }

    // combine the populationPlaceholder with the population
    this->population.insert(this->population.end(),
                            std::make_move_iterator(this->populationPlaceholder.begin()),
                            std::make_move_iterator(this->populationPlaceholder.end()));

    this->populationPlaceholder.clear();
}

std::shared_ptr<evolution::Agent>
evolution::Population::crossoverThreaded(std::uniform_int_distribution<unsigned int> populationDistribution) {
    // choose 2 random agents
    unsigned long firstAgentIndex = populationDistribution(util::rng);
    unsigned long secondAgentIndex;
    // keep rolling the second agent until the index is different from the first agent
    do {
        secondAgentIndex = populationDistribution(util::rng);
    } while (firstAgentIndex == secondAgentIndex);

    // create an empty child agent
    std::shared_ptr<evolution::Agent> childAgent = Agent::create(numberOfInputs, inputLabels,
                                                                 numberOfOutputs,
                                                                 outputLabels);

    std::shared_ptr<evolution::Agent> agent1 = population.at(firstAgentIndex);
    std::shared_ptr<evolution::Agent> agent2 = population.at(secondAgentIndex);

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

    // normalize edge weights
    // childAgent->getGraph()->normalizeEdgeWeights();

    if (!this->keepDormantVerticesAndEdges) {
        childAgent = this->minimizeAgent(childAgent);
    }
    return childAgent;
}

void evolution::Population::mutateThreaded(unsigned long agentIndex) {
    // choose a random agent
    auto childAgent = this->populationPlaceholder.at(agentIndex);

    // sanitize agent


    // choose a property to mutate
    // 0 - add a deep vertex with a random input and output edge (random dominance)
    // 1 - remove a deep vertex and all connecting edges
    // 2 - add a random edge
    // 3 - remove a random edge
    // 4 - change the weight of an edge
    // 5 - change traverse limit of an edge
    // 6 - minimize agent
    unsigned int choice = this->mutationDistribution(util::rng);
    switch (choice) {
        case 0: {
            if (childAgent->getGraph()->getDeepVertices().size() == this->maxDeepVertices) {
                return;
            }
            // create vertex
            auto newDeepVertex = data_structures::DeepVertex::createDeepVertex(std::numeric_limits<unsigned int>::max(), util::nextBool());
            // add vertex
            childAgent->getGraph()->addDeepVertex(newDeepVertex);
            // edges can have childAgent as input or output
            // add input edge
            this->addRandomEdge(childAgent->getGraph(), newDeepVertex, true);
            // add output edge
            this->addRandomEdge(childAgent->getGraph(), newDeepVertex, false);
            // fix indices for vertices
            childAgent->getGraph()->fixIndices();
            break;
        }
        case 1: {
            if (childAgent->getGraph()->getDeepVertices().empty()) {
                return;
            }
            // CONSTRUCT A NEW AGENT
            // choose a random deep vertex
            std::uniform_int_distribution<unsigned long> deepVertexDistribution(0,
                                                                                childAgent->getGraph()->getDeepVertices().size() -
                                                                                1);
            unsigned long positionToRemove = deepVertexDistribution(util::rng);

            // flag vertex and connected edges for deletion
            auto deepVertexToRemove = childAgent->getGraph()->getDeepVertices().at(positionToRemove);
            deepVertexToRemove->setFlaggedForDeletion(true);
            for (const auto &edge: deepVertexToRemove->getInputEdges()) {
                edge->setFlaggedForDeletion(true);
            }
            for (const auto &edge: deepVertexToRemove->getOutputEdges()) {
                edge->setFlaggedForDeletion(true);
            }
            // construct the agent
            std::shared_ptr<Agent> mutatedAgent = Agent::create(numberOfInputs, inputLabels,
                                                                numberOfOutputs,
                                                                outputLabels);

            // add the deep vertices
            for (const auto &deepVertex: childAgent->getGraph()->getDeepVertices()) {
                if (!deepVertex->isFlaggedForDeletion()) {
                    auto vertexClone = deepVertex->deepClone();
                    mutatedAgent->getGraph()->addDeepVertex(vertexClone);
                }
            }
            // add the edges
            for (const auto &edge: childAgent->getGraph()->getEdges()) {
                if (/*edge != nullptr &&*/ !edge->isFlaggedForDeletion()) {
                    mutatedAgent->getGraph()->addEdge(edge);
                }
            }

            // replace the child agent
            this->populationPlaceholder.at(agentIndex) = mutatedAgent;
            break;
        }
        case 2: {
            if (childAgent->getGraph()->getEdges().size() >= maxEdges ||
                childAgent->getGraph()->getEdges().size() >= childAgent->getGraph()->getNumEdgesPossible()) {
                return;
            }
            unsigned long oldEdgeSize = childAgent->getGraph()->getEdges().size();
            while (childAgent->getGraph()->getEdges().size() <= oldEdgeSize) {
                this->addRandomEdge(std::numeric_limits<unsigned int>::max(), childAgent->getGraph());
            }
            // fix indices for randomly added edges
            childAgent->getGraph()->fixIndices();
            break;
        }
        case 3: {
            if (childAgent->getGraph()->getEdges().empty()) {
                return;
            }
            // choose a random edge
            unsigned long position = util::nextUnsignedLong(0, childAgent->getGraph()->getEdges().size() - 1);
            childAgent->getGraph()->removeEdge(position);
            break;
        }
        case 4: {
            if (childAgent->getGraph()->getEdges().empty()) {
                return;
            }
            // choose a random edge
            unsigned long maxPos = childAgent->getGraph()->getEdges().size() - 1;
            unsigned long position = util::nextUnsignedLong(0, maxPos);
            // choose a random weight
            double weight = util::nextWeight();
            // get edge
            auto edge = childAgent->getGraph()->getEdges().at(position);
            edge->setWeight(weight);
            break;
        }
        case 5: {
            if (childAgent->getGraph()->getEdges().empty()) {
                return;
            }
            // choose a random edge
            unsigned long maxPos = childAgent->getGraph()->getEdges().size() - 1;
            unsigned long position = util::nextUnsignedLong(0, maxPos);
            // choose a random weight
            unsigned int traverseLimit = this->edgeTraverseLimitDistribution(util::rng);
            // get edge
            auto edge = childAgent->getGraph()->getEdges().at(position);
            edge->setTraverseLimit(traverseLimit);
            break;
        }
        case 6: {
            auto minimizedAgent = this->minimizeAgent(childAgent);
            // replace the child agent
            this->populationPlaceholder.at(agentIndex) = minimizedAgent;
            break;
        }
        default:
            throw std::invalid_argument("Invalid mutation choice: " + std::to_string(choice));
    }
}

void evolution::Population::crossoverDeepVertices(const std::shared_ptr<evolution::Agent> &childAgent,
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
        // if no second vertex was found, clone vertex1 with a 50% chance
        std::shared_ptr<data_structures::DeepVertex> childVertex;
        if (deepVertex2 == nullptr) {
            bool choice = util::nextBool();
            if (!choice) {
                continue;
            }
            childVertex = deepVertex1->deepClone();
        } else {
            // choose a random vertex
            bool choice = util::nextBool();
            childVertex = choice ? deepVertex1->deepClone() : deepVertex2->deepClone();
            // randomize the child dominant trait to prevent convergence
            childVertex->setDominant(util::nextBool());
        }

        checkedIndices.push_back(deepVertex1->getIndex());
        childAgent->getGraph()->addDeepVertex(childVertex);
    }
}

void evolution::Population::crossoverEdges(const std::shared_ptr<evolution::Agent> &childAgent,
                                           const std::shared_ptr<evolution::Agent> &parent1,
                                           const std::shared_ptr<evolution::Agent> &parent2,
                                           std::vector<unsigned int> &checkedIndices) {
    // iterate over the edges of the first parent
    for (const auto &edge1: parent1->getGraph()->getEdges()) {
        if (childAgent->getGraph()->getEdges().size() >= this->maxEdges) {
            break;
        }

        // check if this edge1 was already crossovered
        //bool matchFound = false;
        //for (unsigned int index: checkedIndices) {
        //    if (edge1->getIndex() == index) {
        //        matchFound = true;
        //        break;
        //    }
        //}
        //if (matchFound) {
        //    continue;
        //}
        // get the edge1 with the same input and output indices (and types) from the second parent
        auto edge2 = parent2->getGraph()->getEdgeByIndexAndType(edge1->getInput()->getIndex(),
                                                                edge1->getInput()->getType(),
                                                                edge1->getOutput()->getIndex(),
                                                                edge1->getOutput()->getType());
        if (edge2 == nullptr) {
            // similar edge1 was not found in parent2
            // 50% chance to keep edge1
            bool choice = util::nextBool();
            if (!choice) {
                continue;
            }
            childAgent->getGraph()->addEdge(edge1->getInput()->getType(), edge1->getInput()->getIndex(),
                                            edge1->getOutput()->getType(), edge1->getOutput()->getIndex(),
                                            edge1->getIndex(), edge1->getWeight(),
                                            edge1->getTraverseLimit(),
                                            data_structures::ICrossoverable(edge1->isDominant()));
        } else {
            // choose a random edge
            bool choice = util::nextBool();
            auto dominantEdge = choice ? edge1 : edge2;
            childAgent->getGraph()->addEdge(dominantEdge->getInput()->getType(), dominantEdge->getInput()->getIndex(),
                                            dominantEdge->getOutput()->getType(), dominantEdge->getOutput()->getIndex(),
                                            dominantEdge->getIndex(), dominantEdge->getWeight(),
                                            dominantEdge->getTraverseLimit(),
                                            data_structures::ICrossoverable(dominantEdge->isDominant()));
        }
        checkedIndices.push_back(edge1->getIndex());
    }
}

std::vector<std::shared_ptr<evolution::Agent>> evolution::Population::getPopulation() const {
    return this->population;
}

std::shared_ptr<evolution::Agent> evolution::Population::getFittestAgent() const {
    return this->fittestAgent;
}

double evolution::Population::getAverageFitness() const {
    return this->averageFittness;
}

std::vector<std::shared_ptr<data_structures::DataInstance>> evolution::Population::getTestingValues() const {
    return this->testingValues;
}

unsigned int evolution::Population::getNumberOfOutputs() const {
    return this->numberOfOutputs;
}

std::vector<std::string> evolution::Population::getOutputLabels() const {
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
            } else {
                // check if all non-flagged input edges are recursive
                bool allInputsRecursive = true;
                for (const auto &edge: deepVertex->getInputEdges()) {
                    if (!edge->isFlaggedForDeletion() &&
                        (edge->getInput()->getIndex() != edge->getOutput()->getIndex() ||
                         edge->getInput()->getType() != edge->getOutput()->getType())) {
                        allInputsRecursive = false;
                        break;
                    }
                }
                // check if all outputs are recursive
                bool allOutputsRecursive = true;
                // if all inputs are already recurse, then this check is redundant
                if (!allInputsRecursive) {
                    for (const auto &edge: deepVertex->getOutputEdges()) {
                        if (!edge->isFlaggedForDeletion() &&
                            (edge->getInput()->getIndex() != edge->getOutput()->getIndex() ||
                             edge->getInput()->getType() != edge->getOutput()->getType())) {
                            allOutputsRecursive = false;
                            break;
                        }
                    }
                }
                // if all inputs or outputs are recursive, then flag vertex and all its edges for deletion
                if (allInputsRecursive || allOutputsRecursive) {
                    deepVertex->setFlaggedForDeletion(true);
                    for (const auto &edge: deepVertex->getInputEdges()) {
                        edge->setFlaggedForDeletion(true);
                    }
                    for (const auto &edge: deepVertex->getOutputEdges()) {
                        edge->setFlaggedForDeletion(true);
                    }
                    keepRunning = true;
                    continue;
                }
            }
        }
    }

    // construct the agent
    std::shared_ptr<Agent> minimizedAgent = Agent::create(numberOfInputs, inputLabels, numberOfOutputs, outputLabels);

    // add the deep vertices
    for (const auto &deepVertex: agentToMinimize->getGraph()->getDeepVertices()) {
        if (!deepVertex->isFlaggedForDeletion()) {
            auto vertexClone = deepVertex->deepClone();
            minimizedAgent->getGraph()->addDeepVertex(vertexClone);
        }
    }
    // add the edges
    for (const auto &edge: agentToMinimize->getGraph()->getEdges()) {
        if (/*edge != nullptr &&*/ !edge->isFlaggedForDeletion()) {
            minimizedAgent->getGraph()->addEdge(edge);
        }
    }

    // re-index deep vertices
    unsigned int newIndex = 0;
    for (const auto &deepVertex: minimizedAgent->getGraph()->getDeepVertices()) {
        deepVertex->setIndex(newIndex);
        newIndex += 1;
    }

    // copy over the fitness, accuracy and mcc
    minimizedAgent->setFitness(agentToMinimize->getFitness());
    minimizedAgent->setAccuracy(agentToMinimize->getAccuracy());
    minimizedAgent->setMatthewsCorrelationCoefficient(agentToMinimize->getMatthewsCorrelationCoefficient());

    return minimizedAgent;
}

double evolution::Population::getAverageAccuracy() const {
    return this->averageAccuracy;
}

double evolution::Population::getAverageMatthewsCorrelationCoefficient() const {
    return this->averageMatthewsCorrelationCoefficient;
}

std::vector<std::string> evolution::Population::getInputLabels() const {
    return this->inputLabels;
}

void evolution::Population::addAgent(const std::shared_ptr<evolution::Agent> &agent) {
    this->population.push_back(agent);
}

void evolution::Population::uncacheFittestAgent() {
    this->fittestAgent = nullptr;
}
