//
// Created by jure on 1. 02. 21.
//

#include <fstream>
#include <filesystem>
#include <algorithm>
#include "Population.h"
#include "../util/util.h"

evolution::Population::Population(const std::string &pathToDataSet) {
    if (!std::filesystem::exists(pathToDataSet)) {
        throw std::invalid_argument("File " + pathToDataSet + " does not exist.");
    }

    std::string delimiter = ",";
    this->trainingValues = std::vector<std::shared_ptr<data_structures::DataInstance>>();

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
                auto dataInstance = data_structures::DataInstance::createDataInstance(util::splitDouble(line, delimiter));
                this->trainingValues.push_back(dataInstance);
            }

            counter++;
        }
        file.close();
    }

    this->population = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationPlaceholder = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationSize = 0;
}

//evolution::Population::~Population() {
//    for (auto graph : *this->population) {
//        delete graph;
//    }
//
//    // should always be empty
//    delete this->populationPlaceholder;
//
//    for (auto value : *this->trainingValues) {
//        delete value;
//    }
//    delete this->inputLabels;
//    delete this->outputLabels;
//    delete this->population;
//    delete this->trainingValues;
//}

void evolution::Population::initialisePopulation(unsigned int populationSize, unsigned int maxDeepVertices,
                                                 unsigned int maxEdges,
                                                 unsigned int edgeTraverseLimit, bool keepDormantVerticesAndEdges) {
    this->population = std::vector<std::shared_ptr<evolution::Agent>>();
    this->populationSize = populationSize;
    for (unsigned int i = 0; i < populationSize; i++) {
        this->population.push_back(
                this->createAgent(maxDeepVertices, maxEdges, edgeTraverseLimit, keepDormantVerticesAndEdges));
    }
}

std::shared_ptr<evolution::Agent>
evolution::Population::createAgent(unsigned int maxDeepVertices, unsigned int maxEdges, unsigned int edgeTraverseLimit,
                                   bool keepDormantVerticesAndEdges) {
    std::mt19937 rng(this->seeder());
    std::mt19937_64 rng64(this->seeder());
    std::uniform_int_distribution<unsigned int> inputVerticesDistribution(0, this->numberOfInputs - 1);
    std::uniform_int_distribution<unsigned int> maxDeepVerticesDistribution(0, maxDeepVertices - 1);
    std::uniform_int_distribution<unsigned int> outputVerticesDistribution(0, this->numberOfOutputs - 1);
    std::uniform_int_distribution<unsigned int> edgeTraverseLimitDistribution(1, edgeTraverseLimit);
    std::uniform_real_distribution<double> weightDistribution(0, 1);
    std::uniform_int_distribution<unsigned int> edgesDistribution(1, maxEdges);
    std::uniform_int_distribution<unsigned int> vertexTypeDistribution(0, 1);
    unsigned int numberOfDeepVertices = maxDeepVerticesDistribution(rng);
    std::uniform_int_distribution<unsigned int> deepVerticesDistribution(0, numberOfDeepVertices - 1);
    unsigned int numberOfEdges = edgesDistribution(rng);

    // generate vertices
    auto graph = data_structures::Graph::createGraph(this->numberOfInputs, this->inputLabels, numberOfDeepVertices,
                                             this->numberOfOutputs, this->outputLabels);

    // generate edges
    unsigned int type1;
    unsigned int type2;
    for (unsigned int i = 0; i < numberOfEdges; i++) {
        // pick two edges, they can be:
        // 1 input and 1 deep
        // 1 input and 1 output
        // 2 deep
        // 1 deep or 1 output

        if (numberOfDeepVertices > 0) {
            // 0 - input
            // 1 - deep
            type1 = vertexTypeDistribution(rng);
            if (type1 == 0) {
                type2 = vertexTypeDistribution(rng);
                // 0 - deep
                // 1 - output
                if (type2 == 0) {
                    graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(rng), enums::VertexType::Deep,
                                   deepVerticesDistribution(rng), weightDistribution(rng64),
                                   edgeTraverseLimitDistribution(rng));
                } else {
                    graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(rng), enums::VertexType::Output,
                                   outputVerticesDistribution(rng), weightDistribution(rng64),
                                   edgeTraverseLimitDistribution(rng));
                }
            } else {
                type2 = vertexTypeDistribution(rng);
                // 0 - deep
                // 1 - output
                if (type2 == 0) {
                    graph->addEdge(enums::VertexType::Deep, deepVerticesDistribution(rng), enums::VertexType::Deep,
                                   deepVerticesDistribution(rng), weightDistribution(rng64),
                                   edgeTraverseLimitDistribution(rng));
                } else {
                    graph->addEdge(enums::VertexType::Deep, deepVerticesDistribution(rng), enums::VertexType::Output,
                                   outputVerticesDistribution(rng), weightDistribution(rng64),
                                   edgeTraverseLimitDistribution(rng));
                }
            }
        } else {
            // no deep vertices
            // only legal configuration is input -> output
            graph->addEdge(enums::VertexType::Input, inputVerticesDistribution(rng), enums::VertexType::Output,
                           outputVerticesDistribution(rng), weightDistribution(rng64),
                           edgeTraverseLimitDistribution(rng));
        }
    }

    return evolution::Agent::create(graph);
}

std::string evolution::Population::toString() {
    std::ostringstream result;
    result << "Number of agents: " << this->population.size() << std::endl;
    for (unsigned int i = 0; i < this->population.size(); i++) {
        result << "Agent " << i << ": " << std::endl;
        result << this->population.at(i)->toString() << std::endl;
    }

    return result.str();
}

void evolution::Population::calculateFitness() {
    // run all data instances on all agents
    unsigned int numCorrect;
    unsigned int numIncorrect;

    for (const std::shared_ptr<evolution::Agent>& agent : this->population) {
        numCorrect = 0;
        numIncorrect = 0;

        for (const std::shared_ptr<data_structures::DataInstance>& di : this->trainingValues) {
            agent->getGraph()->traverse(di);

            // check if the prediction is correct
            if (agent->getGraph()->getLargestOutputValueIndex() == di->getCorrectIndex()) {
                numCorrect += 1;
            } else {
                numIncorrect += 1;
            }

            // reset the agent
            agent->getGraph()->reset();
        }
        // calculate the fitness (accuracy)
        agent->setFitness((double) numCorrect / (double) this->trainingValues.size());
    }
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
    for (int i = 0; i < this->populationSize; i++) {
        if (std::find(indexesToKeep.begin(), indexesToKeep.end(), i) != indexesToKeep.end()) {
            this->populationPlaceholder.push_back(this->population.at(i));
        }
    }

    for(unsigned int index : indexesToKeep){
        // indexesToKeep will never be larger that the population
        this->populationPlaceholder.push_back(evolution::Agent::create(this->population.at(index)));
    }

    this->population->swap(*this->populationPlaceholder);
    this->populationPlaceholder->clear();
}

std::vector<unsigned int> evolution::Population::stochasticUniversalSampling(unsigned int agentsToKeep) {
    // calculate the total fitness of the population
    double fitnessSum = 0;
    for (evolution::Agent *agent : *this->population) {
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
    std::vector<unsigned int> indexesToKeep{};

    // loop until the vector of agents is filled
    double pointer = start;
    int counter = 0;
    double cumulativeFitness = this->population->at(0)->getFitness();
    while (indexesToKeep.size() < agentsToKeep) {
        if (cumulativeFitness >= pointer) {
            indexesToKeep.push_back(counter);
            pointer += distance;
        } else {
            counter += 1;
            cumulativeFitness += this->population->at(counter)->getFitness();
        }
    }

    // indexes will always be sorted
    // std::sort(indexesToKeep.begin(), indexesToKeep.end());

    return indexesToKeep;
}

void evolution::Population::crossover() {
    // initialize randomness
    std::mt19937 rng(this->seeder());
    // choose the maximum size here since we don't want to do crossover with child agents
    std::uniform_int_distribution<unsigned int> populationDistribution(0, this->population->size());

    // create new agents until the population is full
    while (this->population->size() <= this->populationSize) {
        // choose 2 random agents
        // the same agent can be rolled twice
        evolution::Agent *agent1 = this->population->at(populationDistribution(rng));
        evolution::Agent *agent2 = this->population->at(populationDistribution(rng));
        auto *childAgent = new evolution::Agent();
        // childAgent->getGraph()->addInputVertices()

        // output vertices must be fixed, so only do crossover on input and deep vertices
        // only select the dominant vertices from both agents

        // input vertices
        for (int i = 0; i < agent1->getGraph()->getInputVertices()->size(); i++) {
            data_structures::InputVertex *v1 = agent1->getGraph()->getInputVertices()->at(i);
            data_structures::InputVertex *v2 = agent2->getGraph()->getInputVertices()->at(i);
            // get the dominant vertex
            auto* dominant = (data_structures::DeepVertex *)this->getDominantVertex(v1, v2);
            // add the dominant vertex to the child agent
            //childAgent->getGraph()->addInputVertices();
        }

        // deep vertices
        unsigned int leastDeepVertices = std::min(agent1->getGraph()->getDeepVertices()->size(),
                                                  agent2->getGraph()->getDeepVertices()->size());
        for (int i = 0; i < leastDeepVertices; i++) {
            data_structures::DeepVertex *v1 = agent1->getGraph()->getDeepVertices()->at(i);
            data_structures::DeepVertex *v2 = agent2->getGraph()->getDeepVertices()->at(i);
            // get the dominant vertex
            auto* dominant = (data_structures::DeepVertex *)this->getDominantVertex(v1, v2);
            // add the dominant vertex to the child agent
            //childAgent->getGraph().add;
        }

    }
}

data_structures::Vertex * evolution::Population::getDominantVertex(data_structures::Vertex *v1, data_structures::Vertex *v2) {
    if (v1->isDominant() && v2->isDominant() || !v1->isDominant() && !v2->isDominant()) {
        // if both vertices are dominant or recessive, then choose the one with the least chance to get dominated
        return v1->getChanceToGetDominated() > v2->getChanceToGetDominated() ? v2 : v1;
    } else if (v1->isDominant()) {
        // only v1 is dominant, however there's a chance it's not
        return v1->getChanceToGetDominated() > util::nextDouble() ? v2 : v1;
    } else {
        // only v2 is dominant, however there's a chance it's not
        return v2->getChanceToGetDominated() > util::nextDouble() ? v1 : v2;
    }
}
