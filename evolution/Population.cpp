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
    this->inputVerticesDistribution = std::uniform_int_distribution<unsigned long>(0, numberOfInputs - 1);
    this->outputVerticesDistribution = std::uniform_int_distribution<unsigned long>(0, numberOfOutputs - 1);
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
    this->edgeTraverseLimitDistribution = std::uniform_int_distribution<unsigned int>(1, edgeTraverseLimitArg);
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
    graph->normalizeEdgeWeights();

    auto agent = evolution::Agent::create(graph);
    //if (!this->keepDormantVerticesAndEdges) {
    //    agent = this->minimizeAgent(agent);
    //}
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
            unsigned long position = this->inputVerticesDistribution(util::rng);
            inputVertex = graph->getInputVertices().at(position);
        } else {
            unsigned long position = util::nextUnsignedLong(0, graph->getDeepVertices().size() - 1);
            inputVertex = graph->getDeepVertices().at(position);
        }
        // create and add edge
        graph->addEdge(inputVertexType, inputVertex->getIndex(), enums::VertexType::Deep, deepVertex->getIndex(),
                       UINT_MAX, weight, traverseLimit, this->maxMutationChance);
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
                       outputVertex->getIndex(), UINT_MAX, weight, traverseLimit,
                       this->maxMutationChance);
    }
}

void evolution::Population::addRandomEdge(unsigned int index, std::shared_ptr<data_structures::Graph> const &graph) {
    std::mt19937 rng(this->seeder());
    std::uniform_int_distribution<unsigned int> vertexTypeDistribution(0, 1);
    std::uniform_int_distribution<unsigned int> deepVerticesDistribution(0, graph->getDeepVertices().size() - 1);

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
                               deepVerticesDistribution(rng), index, util::nextWeight(),
                               edgeTraverseLimitDistribution(rng), this->maxMutationChance);
            } else {
                graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(rng), enums::VertexType::Output,
                               outputVerticesDistribution(rng), index, util::nextWeight(),
                               edgeTraverseLimitDistribution(rng), this->maxMutationChance);
            }
        } else {
            type2 = vertexTypeDistribution(rng);
            // 0 - deep
            // 1 - output
            if (type2 == 0) {
                graph->addEdge(enums::VertexType::Deep, deepVerticesDistribution(rng), enums::VertexType::Deep,
                               deepVerticesDistribution(rng), index, util::nextWeight(),
                               edgeTraverseLimitDistribution(rng), this->maxMutationChance);
            } else {
                graph->addEdge(enums::VertexType::Deep, deepVerticesDistribution(rng), enums::VertexType::Output,
                               outputVerticesDistribution(rng), index, util::nextWeight(),
                               edgeTraverseLimitDistribution(rng), this->maxMutationChance);
            }
        }
    } else {
        // no deep vertices
        // only legal configuration is input -> output
        graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(rng), enums::VertexType::Output,
                       outputVerticesDistribution(rng), index, util::nextWeight(),
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
        //calculateAgentFitness(fitnessMetric, vertexContribution, edgeContribution, std::ref(agent));
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

    // size of the agent is a penalty, both edges and vertices contribute
    //unsigned int maxTraverseSum = 0;
    //for (const auto &edge: agent->getGraph()->getEdges()) {
    //    maxTraverseSum += edge->getTraverseLimit();
    //}
    //double sizeContribution = (double) agent->getGraph()->getDeepVertices().size() * vertexContribution +
    //                          (double) maxTraverseSum * edgeContribution;
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
    agent->setAccuracy(mcm.getAccuracy());
    agent->setMatthewsCorrelationCoefficient(mcm.getMatthewsCorrelationCoefficient());
}

void evolution::Population::sample(enums::SelectionType type, unsigned int agentsToKeep, bool keepFittest) {
    std::vector<unsigned int> indexesToKeep{};
    if (this->populationSize < agentsToKeep) {
        throw std::invalid_argument("Number of agents to keep cannot be higher than the population size.");
    }

    if (type == enums::SelectionType::StochasticUniversalSampling) {
        indexesToKeep = stochasticUniversalSampling(agentsToKeep, keepFittest);
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
evolution::Population::stochasticUniversalSampling(unsigned int agentsToKeep, bool keepFittest) {
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
    std::vector<unsigned int> indexesToKeep;
    if (keepFittest) {
        indexesToKeep.push_back(fittestAgentIndex);
    }
    //std::vector<unsigned int> indexesToKeep;

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

void evolution::Population::crossoverAndMutate() {
    // sort the population according to fitness
    // std::sort(this->population.begin(), this->population.end(), sortByFitness);

    // create new agents until the population is full
    // save child agents in populationPlaceholder
    unsigned int agentsToCreate = this->populationSize - this->population.size();
    std::vector<std::future<std::shared_ptr<evolution::Agent>>> crossoverFutures;
    // create the population distribution here since it's the same for all threads
    std::uniform_int_distribution<unsigned int> populationDistribution(0, this->population.size() - 1);

    for (int i = 0; i < agentsToCreate; i++) {
        auto ftr = std::async(&evolution::Population::crossoverThreaded, this, populationDistribution);
        crossoverFutures.push_back(std::move(ftr));
        //this->crossoverThreaded();
    }
    for (auto &ftr: crossoverFutures) {
        this->populationPlaceholder.push_back(ftr.get());
    }

    // mutate children
    // mutate the amount of children proportional to the global mutation chance
    std::vector<std::future<void>> mutationFutures;
    unsigned int childrenToMutate = (unsigned int) std::round(agentsToCreate * this->maxMutationChance);
    std::uniform_int_distribution<unsigned int> childPopulationDistribution(0, agentsToCreate - 1);
    // choose unique agent indexes to mutate
    std::vector<unsigned long> indexesToMutate;
    while (indexesToMutate.size() < childrenToMutate) {
        unsigned long index = childPopulationDistribution(util::rng);
        for (const auto &generatedIndex: indexesToMutate) {
            if (generatedIndex == index) {
                continue;
            }
        }
        indexesToMutate.push_back(index);
    }

    for (const auto &index: indexesToMutate) {
        //this->mutateThreaded(index);
        auto ftr = std::async(&evolution::Population::mutateThreaded, this, index);
        mutationFutures.push_back(std::move(ftr));
    }
    for (auto &ftr: mutationFutures) {
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
    //if (childAgent->getGraph()->getEdges().size() < this->maxEdges) {
    //    addNewRandomEdges(childAgent);
    //}

    // fix indices for vertices that are "twins" and randomly added edges
    //childAgent->getGraph()->fixIndices();

    // populationPlaceholder.push_back(childAgent);

    // normalize edge weights
    childAgent->getGraph()->normalizeEdgeWeights();

    if (!this->keepDormantVerticesAndEdges) {
        childAgent = this->minimizeAgent(childAgent);
    }
    return childAgent;
}

void evolution::Population::mutateThreaded(unsigned long agentIndex) {
    // choose a random agent
    auto childAgent = this->populationPlaceholder.at(agentIndex);

    // choose a property to mutate
    // 0 - add a deep vertex with a random input and output edge (random dominance)
    // 1 - remove a deep vertex and all connecting edges
    // 2 - add a random edge
    // 3 - remove a random edge
    // 4 - change the weight of an edge
    // 5 - change traverse limit of an edge
    // 6 - minimize agent
    unsigned int choice = this->mutationDistribution(util::rng);
    //unsigned int choice = 6;

    //unsigned int choice = 1;
    switch (choice) {
        case 0: {
            if (childAgent->getGraph()->getDeepVertices().size() == this->maxDeepVertices) {
                return;
            }
            // create vertex
            auto newDeepVertex = data_structures::DeepVertex::createDeepVertex(UINT_MAX, util::nextBool());
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
                                                                outputLabels, maxMutationChance);

            // add the deep vertices
            for (const auto &deepVertex: childAgent->getGraph()->getDeepVertices()) {
                if (!deepVertex->isFlaggedForDeletion()) {
                    auto vertexClone = deepVertex->deepClone();
                    mutatedAgent->getGraph()->addDeepVertex(vertexClone);
                }
            }
            // add the edges
            for (const auto &edge: childAgent->getGraph()->getEdges()) {
                if (!edge->isFlaggedForDeletion()) {
                    mutatedAgent->getGraph()->addEdge(edge);
                }
            }

            // replace the child agent
            this->populationPlaceholder.at(agentIndex) = mutatedAgent;

            //childAgent->getGraph()->removeDeepVertex(positionToRemove);
            break;
        }
        case 2: {
            if (childAgent->getGraph()->getEdges().size() == maxEdges) {
                return;
            }
            this->addRandomEdge(UINT_MAX, childAgent->getGraph());
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
            unsigned long position = util::nextUnsignedLong(0, childAgent->getGraph()->getEdges().size() - 1);
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
            unsigned long position = util::nextUnsignedLong(0, childAgent->getGraph()->getEdges().size() - 1);
            // choose a random weight
            unsigned int traverseLimit = this->edgeTraverseLimitDistribution(util::rng);
            // get edge
            auto edge = childAgent->getGraph()->getEdges().at(position);
            edge->setTraverseLimit(traverseLimit);
            break;
        }
        case 6: {
            auto minimizedAgent = this->minimizeAgent(childAgent);
            childAgent = minimizedAgent;
            break;
        }
        default:
            throw std::invalid_argument("Invalid mutation choice: " + std::to_string(choice));
    }
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
        // if no second vertex was found, clone vertex1
        std::shared_ptr<data_structures::DeepVertex> childVertex;
        if (deepVertex2 == nullptr) {
            childVertex = deepVertex1->deepClone();
        } else {
            // if vertex2 was found, then choose a dominant vertex
            // if both vertices are dominant or recessive, then always pick the first one
            childVertex = deepVertex1->isDominant() ? deepVertex1->deepClone() : deepVertex2->deepClone();
            // randomize the child dominant trait to prevent convergence
            childVertex->setDominant(util::nextBool());
        }

        checkedIndices.push_back(deepVertex1->getIndex());
        childAgent->getGraph()->addDeepVertex(childVertex);
    }
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
            // keep edge1
            childAgent->getGraph()->addEdge(edge1->getInput()->getType(), edge1->getInput()->getIndex(),
                                            edge1->getOutput()->getType(), edge1->getOutput()->getIndex(),
                                            edge1->getIndex(), edge1->getWeight(),
                                            edge1->getTraverseLimit(),
                                            data_structures::ICrossoverable(edge1->isDominant()));
        } else {
            // if edge2 was found, then choose a dominant vertex
            // if both edges are dominant or recessive, then always pick the first one
            auto dominantEdge = edge1->isDominant() ? edge1 : edge2;
            childAgent->getGraph()->addEdge(dominantEdge->getInput()->getType(), dominantEdge->getInput()->getIndex(),
                                            dominantEdge->getOutput()->getType(), dominantEdge->getOutput()->getIndex(),
                                            dominantEdge->getIndex(), dominantEdge->getWeight(),
                                            dominantEdge->getTraverseLimit(),
                                            data_structures::ICrossoverable(dominantEdge->isDominant()));
            //this->createAndAddEdgeChildren(dominant, childAgent);
        }
        checkedIndices.push_back(edge1->getIndex());
    }
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
                for (const auto &edge: deepVertex->getOutputEdges()) {
                    if (!edge->isFlaggedForDeletion() &&
                        (edge->getInput()->getIndex() != edge->getOutput()->getIndex() ||
                         edge->getInput()->getType() != edge->getOutput()->getType())) {
                        allOutputsRecursive = false;
                        break;
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
    //double maxWeight = 0;
    //for (const auto &edge: minimizedAgent->getGraph()->getEdges())
    //    if (edge->getWeight() > maxWeight) {
    //        maxWeight = edge->getWeight();
    //    }


    // re-index deep vertices
    //unsigned int newIndex = 0;
    //for (const auto &deepVertex: minimizedAgent->getGraph()->getDeepVertices()) {
    //    deepVertex->setIndex(newIndex);
    //    newIndex += 1;
    //}

    // copy over the fitness, accuracy and mcc
    minimizedAgent->setFitness(agentToMinimize->getFitness());
    minimizedAgent->setAccuracy(agentToMinimize->getAccuracy());
    minimizedAgent->setMatthewsCorrelationCoefficient(agentToMinimize->getMatthewsCorrelationCoefficient());


    return minimizedAgent;
}

std::vector<std::string> evolution::Population::getInputLabels() {
    return this->inputLabels;
}

double evolution::Population::getAverageAccuracy() {
    double average = 0;
    for (const auto &agent: this->population) {
        average += agent->getAccuracy();
    }
    return average / (double) this->population.size();
}

double evolution::Population::getAverageMatthewsCorrelationCoefficient() {
    double average = 0;
    for (const auto &agent: this->population) {
        average += agent->getMatthewsCorrelationCoefficient();
    }
    return average / (double) this->population.size();
}
