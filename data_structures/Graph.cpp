//
// Created by jure on 7/13/20.
//

#include "Graph.h"
#include <vector>
#include <string>
#include <sstream>
#include "../enums/EnumUtil.h"
#include "../logging/logging.h"
#include "../util/util.h"
#include "../nlohmann/json.hpp"


data_structures::Graph::Graph(unsigned int inputVertices, std::vector<std::string> &inputLabels,
                              unsigned int deepVertices, unsigned int outputVertices,
                              std::vector<std::string> &outputLabels, double maxMutationChance) {
    this->addInputVertices(inputVertices, inputLabels);
    this->addOutputVertices(outputVertices, outputLabels);
    this->addDeepVertices(deepVertices, maxMutationChance);
}

void data_structures::Graph::addInputVertices(unsigned int numberOfInputVertices, std::vector<std::string> &labels) {
    if (numberOfInputVertices != labels.size()) {
        throw std::invalid_argument("The number of input vertices and labels does not match.");
    }
    for (unsigned int i = 0; i < numberOfInputVertices; i++) {
        auto inputVertex = data_structures::InputVertex::createInputVertex(i, labels.at(i));
        this->inputVertices.push_back(inputVertex);
    }
}

void data_structures::Graph::addInputVertex(const std::shared_ptr<data_structures::InputVertex> &inputVertex) {
    this->inputVertices.push_back(inputVertex);
}

void data_structures::Graph::addOutputVertices(unsigned int numberOfOutputVertices, std::vector<std::string> &labels) {
    for (unsigned int i = 0; i < numberOfOutputVertices; i++) {
        this->outputVertices.push_back(data_structures::OutputVertex::createOutputVertex(i, labels.at(i)));
    }
}

void data_structures::Graph::addOutputVertex(const std::shared_ptr<data_structures::OutputVertex> &outputVertex) {
    this->outputVertices.push_back(outputVertex);
}

void data_structures::Graph::addDeepVertices(unsigned int numberOfDeepVertices, double maxMutationChance) {
    for (unsigned int i = 0; i < numberOfDeepVertices; i++) {
        // initialise chances as values from 0 to 1
        auto deepVertex = data_structures::DeepVertex::createDeepVertex(i, util::nextBool(),
                                                                        util::nextDouble(maxMutationChance),
                                                                        util::nextUnsignedInt(1, 2));
        this->deepVertices.push_back(deepVertex);
        this->largestDeepVertexIndex = i;
    }
}

void data_structures::Graph::addDeepVertices(
        const std::vector<std::shared_ptr<data_structures::DeepVertex>> &deepVerticesArg) {
    for (const auto &vertex: deepVerticesArg) {
        this->addDeepVertex(vertex);
    }
}

void data_structures::Graph::addDeepVertex(const std::shared_ptr<data_structures::DeepVertex> &deepVertex) {
    // assume that there are not vertices with the same index
    if (this->largestDeepVertexIndex < deepVertex->getIndex() && deepVertex->getIndex() != UINT_MAX) {
        this->largestDeepVertexIndex = deepVertex->getIndex();
    }
    this->deepVertices.push_back(deepVertex);
}

void data_structures::Graph::addEdge(enums::VertexType inputVertexType, unsigned int inputVertexIndex,
                                     enums::VertexType outputVertexType, unsigned int outputVertexIndex,
                                     unsigned int index, double weight, unsigned int traversalLimit,
                                     double maxMutationChance) {
    std::shared_ptr<data_structures::Vertex> input;
    std::shared_ptr<data_structures::Vertex> output;
    std::shared_ptr<data_structures::Edge> commonEdge;
    // if the result of the  function is false, then don't add the edge
    if (!edgeBeforeAdd(inputVertexType, inputVertexIndex, outputVertexType, outputVertexIndex, traversalLimit, input,
                       output,
                       commonEdge)) {
        return;
    }

    // if a common edge has not been found, create a new one
    // otherwise, combine the weights
    if (commonEdge == nullptr) {
        // assume that there are no edges with the same index
        if (this->largestEdgeIndex < index && index != UINT_MAX) {
            this->largestEdgeIndex = index;
        }
        std::shared_ptr<data_structures::Edge> newEdge = data_structures::Edge::createEdge(input, output,
                                                                                           index,
                                                                                           weight,
                                                                                           traversalLimit,
                                                                                           util::nextDouble(
                                                                                                   maxMutationChance),
                                                                                           util::nextBool(),
                                                                                           util::nextUnsignedInt(1, 2));
        this->edges.push_back(newEdge);
        addEdgeAfterAdd(inputVertexType, outputVertexType, input, output, newEdge);
    } else {
        commonEdge->combineWeight(weight);
    }
}

void data_structures::Graph::addEdgeAfterAdd(const enums::VertexType &inputVertexType,
                                             const enums::VertexType &outputVertexType,
                                             const std::shared_ptr<data_structures::Vertex> &input,
                                             const std::shared_ptr<data_structures::Vertex> &output,
                                             const std::shared_ptr<data_structures::Edge> &newEdge) const {
    switch (inputVertexType) {
        case enums::VertexType::Input:
        case enums::VertexType::Deep:
            input->addOutputEdge(newEdge);
            break;
        default:
            throw std::invalid_argument("impossible");
    }

    switch (outputVertexType) {
        case enums::VertexType::Deep:
        case enums::VertexType::Output:
            output->addInputEdge(newEdge);
            break;
        default:
            throw std::invalid_argument("impossible");
    }
}

bool data_structures::Graph::edgeBeforeAdd(enums::VertexType &inputVertexType, unsigned int inputVertexIndex,
                                           enums::VertexType &outputVertexType, unsigned int outputVertexIndex,
                                           unsigned int traversalLimit, std::shared_ptr<data_structures::Vertex> &input,
                                           std::shared_ptr<data_structures::Vertex> &output,
                                           std::shared_ptr<data_structures::Edge> &commonEdge) {
    commonEdge = nullptr;
    if (traversalLimit < 1) {
        throw std::invalid_argument(
                "Traversal limit has to be at least 1. The supplied value is " + std::to_string(traversalLimit) + ".");
    }
    switch (inputVertexType) {
        case enums::VertexType::Input:
            input = inputVertices.at(inputVertexIndex);
            break;
        case enums::VertexType::Deep:
            input = this->getDeepVertexByIndex(inputVertexIndex);
            break;
        default:
            throw std::invalid_argument("Vertex of type " + enums::EnumUtil::VertexTypeToString(inputVertexType) +
                                        " cannot be used as an input vertex.");
    }

    // if input vertex doesn't exist, return false
    if (input == nullptr) {
        return false;
    }

    switch (outputVertexType) {
        case enums::VertexType::Output:
            output = outputVertices.at(outputVertexIndex);
            break;
        case enums::VertexType::Deep:
            output = this->getDeepVertexByIndex(outputVertexIndex);
            break;
        default:
            throw std::invalid_argument("Vertex of type " + enums::EnumUtil::VertexTypeToString(outputVertexType) +
                                        " cannot be used as an output vertex.");
    }

    // if output vertex doesn't exist, return false
    if (output == nullptr) {
        return false;
    }

    // check if input and output are already connected
    for (const std::shared_ptr<Edge> &e: input->getOutputEdges()) {
        if (e->getOutput() == output) {
            commonEdge = e;
            break;
        }
    }

    // both input and output vertices exist, return true
    return true;

}

std::shared_ptr<data_structures::Edge>
data_structures::Graph::addEdge(enums::VertexType inputVertexType, unsigned int inputVertexIndex,
                                enums::VertexType outputVertexType, unsigned int outputVertexIndex, unsigned int index,
                                double weight, unsigned int traversalLimit,
                                const data_structures::ICrossoverable &crossoverableData) {
    std::shared_ptr<data_structures::Vertex> input;
    std::shared_ptr<data_structures::Vertex> output;
    std::shared_ptr<data_structures::Edge> commonEdge;
    // if the result of the  function is false, then don't add the edge
    if (!edgeBeforeAdd(inputVertexType, inputVertexIndex, outputVertexType, outputVertexIndex, traversalLimit, input,
                       output, commonEdge)) {
        return nullptr;
    }

    // if a common edge has not been found, create a new one
    // otherwise, combine the weights
    if (commonEdge == nullptr) {
        // assume that there are no edges with the same index
        if (this->largestEdgeIndex < index && index != UINT_MAX) {
            this->largestEdgeIndex = index;
        }
        std::shared_ptr<data_structures::Edge> newEdge = data_structures::Edge::createEdge(input, output, index, weight,
                                                                                           traversalLimit,
                                                                                           crossoverableData.getMutationChance(),
                                                                                           crossoverableData.isDominant(),
                                                                                           crossoverableData.getMaxChildren());
        this->edges.push_back(newEdge);
        addEdgeAfterAdd(inputVertexType, outputVertexType, input, output, newEdge);
        return newEdge;
    } else {
        commonEdge->combineWeight(weight);
        return commonEdge;
    }
}


void data_structures::Graph::traverse(const std::shared_ptr<data_structures::DataInstance> &dataInstance) {
    if (dataInstance->getValues().size() != this->inputVertices.size()) {
        throw std::invalid_argument("Number of values in dataInstance does not match the number of input vertices.");
    }

    // assign values to input vertices
    for (unsigned int i = 0; i < this->inputVertices.size(); i++) {
        this->inputVertices.at(i)->combineValue(dataInstance->getValues().at(i));
    }

    // add the input vertices to the set of pending vertices
    // enqueue them only if the have any output edges
    for (const auto &v: this->inputVertices) {
        if (!v->getOutputEdges().empty()) {
            this->pendingVertices.enqueue(v);
        }
    }

    // propagate the values through the graph
    // repeat until the set of pending vertices is empty
    std::shared_ptr<data_structures::Vertex> vertex;
    std::vector<std::shared_ptr<data_structures::Edge>> outputEdges;
    while (!this->pendingVertices.empty()) {
        vertex = this->pendingVertices.dequeue();
        vertex->setVisited(true);

        // skip deep vertices without input edges
        if (vertex->getType() == enums::VertexType::Deep && vertex->getInputEdges().empty()) {
            continue;
        }

        // vertex->getOutputEdges() cannot be empty at this stage
        for (const std::shared_ptr<data_structures::Edge> &edge: vertex->getOutputEdges()) {
            if (!edge->isAtTaversalLimit()) {
                // propagate the values to adjacent vertices and add the output vertex to the pending vertices set
                // if it has not been visited yet
                edge->propagateValue();
                // only add the output vertex if:
                // - all input edges of the output vertex have been traversed (the value stored in the output vertex is final)
                // - output vertex has not yet been visited
                // - output vertex has at lest one output edge
                // - cyclic connections are handled further on
                if (edge->getOutput()->getType() == enums::VertexType::Deep &&
                    edge->getOutput()->allInputEdgesTraversed() && !edge->getOutput()->isVisited() &&
                    !edge->getOutput()->getOutputEdges().empty()) {
                    this->pendingVertices.enqueue(edge->getOutput());
                }
            }
        }

        // in case of recursive connections, check if any edge has not been traversed. If such an edge is found, then
        // check if its output vertex has any input edges that have been traversed and enqueue it (should always be a deep vertex)
        // otherwise, end
        if (!this->pendingVertices.empty()) {
            continue;
        }
        std::shared_ptr<data_structures::Vertex> vertexToEnqueue = nullptr;
        unsigned int leastNonTraversedEdges = UINT32_MAX;
//        for (const std::shared_ptr<data_structures::Edge> &edge: this->edges) {
//            if (edge->isTraversed() || edge->getOutput()->getType() != enums::VertexType::Deep) {
//                continue;
//            }
//            for (const std::shared_ptr<data_structures::Edge> &outputVertexEdge: edge->getOutput()->getInputEdges()) {
//                if (outputVertexEdge->isTraversed()) {
//                    vertexToEnqueue = edge->getOutput();
//                    break;
//                }
//            }
//            if (vertexToEnqueue != nullptr) {
//                break;
//            }
//        }

        // in case of recursive connections, check if any deep vertex with out/input vertices has not been visited.
        // If such a vertex is found, then check if it has any input edges that have been traversed and save it
        // enqueue the vertex with the least amount of non-traversed edges
        // otherwise, end
        for (const std::shared_ptr<data_structures::DeepVertex> &deepVertex: this->deepVertices) {
            if (deepVertex->isVisited() || deepVertex->getOutputEdges().empty() ||
                deepVertex->getInputEdges().empty()) {
                continue;
            }
            // count non-traversed edges
            unsigned int numNonTraversed = 0;
            for (const std::shared_ptr<data_structures::Edge> &inputEdge: deepVertex->getInputEdges()) {
                if (!inputEdge->isTraversed()) {
                    numNonTraversed += 1;
                }
            }
            // all edges have been traversed, continue
            if (numNonTraversed == 0) {
                continue;
            }
            // handle the computed value
            if (leastNonTraversedEdges > numNonTraversed) {
                leastNonTraversedEdges = numNonTraversed;
                vertexToEnqueue = deepVertex;
            }
        }

        if (vertexToEnqueue != nullptr) {
            this->pendingVertices.enqueue(vertexToEnqueue);
        }
    }
}

unsigned int data_structures::Graph::getLargestOutputValueIndex() {
    unsigned int largestIndex = 0;
    double largestValue = -std::numeric_limits<double>::max();
    for (const std::shared_ptr<data_structures::OutputVertex> &v: this->outputVertices) {
        if (v->getValue() > largestValue) {
            largestValue = v->getValue();
            largestIndex = v->getIndex();
        }
    }

    //if (largestValue == 0) {
    //    return -1;
    //}

    return largestIndex;
}

void data_structures::Graph::reset() {
    for (const std::shared_ptr<data_structures::InputVertex> &v: this->inputVertices) {
        v->reset();
    }
    for (const std::shared_ptr<data_structures::DeepVertex> &v: this->deepVertices) {
        v->reset();
    }
    for (const std::shared_ptr<data_structures::OutputVertex> &v: this->outputVertices) {
        v->reset();
    }

    for (std::shared_ptr<data_structures::Edge> e: this->edges) {
        e->reset();
    }

    this->pendingVertices.clear();
}

std::vector<std::shared_ptr<data_structures::DeepVertex>> data_structures::Graph::getDeepVertices() {
    return this->deepVertices;
}

std::vector<std::shared_ptr<data_structures::InputVertex >> data_structures::Graph::getInputVertices() {
    return this->inputVertices;
}

std::vector<std::shared_ptr<data_structures::OutputVertex>> data_structures::Graph::getOutputVertices() {
    return this->outputVertices;
}

std::string data_structures::Graph::toString(bool technical) {
    std::ostringstream result;
    result << "Input vertices: " << std::to_string(this->inputVertices.size()) << std::endl;
    result << "Output vertices: " << std::to_string(this->outputVertices.size()) << std::endl;
    result << "Deep vertices: " << std::to_string(this->deepVertices.size()) << std::endl;
    result << "Edges: " << std::to_string(this->edges.size());

    for (const std::shared_ptr<data_structures::Edge> &edge: this->edges) {
        if (technical) {
            result << std::endl << edge->toString(technical);
        } else {
            result << std::endl << "\t" << edge->toString(technical);
        }
    }

    return result.str();
}

std::shared_ptr<data_structures::Graph>
data_structures::Graph::createGraph(unsigned int inputVertices, std::vector<std::string> &inputLabels,
                                    unsigned int deepVertices, unsigned int outputVertices,
                                    std::vector<std::string> &outputLabels, double maxMutationChance) {
    return std::make_shared<data_structures::Graph>(inputVertices, inputLabels, deepVertices, outputVertices,
                                                    outputLabels, maxMutationChance);
}

std::shared_ptr<data_structures::Graph> data_structures::Graph::deepClone() {
    /*
     * unsigned int inputVertices, std::vector<std::string> &inputLabels,
                              unsigned int deepVertices,
                              unsigned int outputVertices, std::vector<std::string> &outputLabels
     */
    // make a new empty instance
    auto newGraph = std::make_shared<data_structures::Graph>();

    // clone the input vertices
    for (const std::shared_ptr<data_structures::InputVertex> &inputVertex: this->inputVertices) {
        newGraph->addInputVertex(inputVertex->deepClone());
    }

    // clone output vertices
    for (const std::shared_ptr<data_structures::OutputVertex> &outputVertex: this->outputVertices) {
        auto deepVertex = outputVertex->deepClone();
        newGraph->addOutputVertex(deepVertex);
    }

    // clone deep vertices
    for (const std::shared_ptr<data_structures::DeepVertex> &deepVertex: this->deepVertices) {
        newGraph->addDeepVertex(deepVertex->deepClone());
    }

    // recreate edges of input vertices
    for (const std::shared_ptr<data_structures::InputVertex> &inputVertex: this->inputVertices) {
        // input vertices have no input edges
        for (const std::shared_ptr<data_structures::Edge> &edge: inputVertex->getOutputEdges()) {
            data_structures::ICrossoverable crossoverable(edge->getMutationChance(), edge->isDominant(),
                                                          edge->getMaxChildren());
            newGraph->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(),
                              edge->getOutput()->getType(), edge->getOutput()->getIndex(),
                              edge->getIndex(), edge->getWeight(), edge->getTraverseLimit(), crossoverable);
        }
    }

    // recreate edges of deep vertices
    for (const std::shared_ptr<data_structures::DeepVertex> &deepVertex: this->deepVertices) {
        for (const std::shared_ptr<data_structures::Edge> &edge: deepVertex->getOutputEdges()) {
            data_structures::ICrossoverable crossoverable(edge->getMutationChance(), edge->isDominant(),
                                                          edge->getMaxChildren());
            newGraph->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(),
                              edge->getOutput()->getType(), edge->getOutput()->getIndex(),
                              edge->getIndex(), edge->getWeight(), edge->getTraverseLimit(), crossoverable);
        }
    }

    // output vertices have no output edges

    // pendingVertices stays empty
    return newGraph;
}

unsigned int data_structures::Graph::getLargestDeepVertexIndex() {
    if (this->largestDeepVertexIndex == -1) {
        for (const auto &vertex: this->deepVertices) {
            if (this->largestDeepVertexIndex < vertex->getIndex()) {
                this->largestDeepVertexIndex = vertex->getIndex();
            }
        }
    }
    return this->largestDeepVertexIndex;
}

void data_structures::Graph::fixIndices() {
    for (const auto &deepVertex: this->deepVertices) {
        if (deepVertex->getIndex() == UINT_MAX) {
            this->largestDeepVertexIndex += 1;
            deepVertex->setIndex(this->largestDeepVertexIndex);
        }
    }
    for (const auto &edge: this->edges) {
        if (edge->getIndex() == UINT_MAX) {
            this->largestEdgeIndex += 1;
            edge->setIndex(this->largestEdgeIndex);
        }
    }
}

std::vector<std::shared_ptr<data_structures::Edge>> data_structures::Graph::getEdges() const {
    return this->edges;
}

std::shared_ptr<data_structures::Edge>
data_structures::Graph::getEdgeByIndexAndType(unsigned int inputVertexIndex, enums::VertexType inputVertexType,
                                              unsigned int outputVertexIndex, enums::VertexType outputVertexType) {
    for (const auto &edge: this->edges) {
        if (edge->getInput()->getIndex() == inputVertexIndex &&
            edge->getInput()->getType() == inputVertexType &&
            edge->getOutput()->getIndex() == outputVertexIndex &&
            edge->getOutput()->getType() == outputVertexType) {
            return edge;
        }
    }
    return nullptr;
}

std::shared_ptr<data_structures::DeepVertex> data_structures::Graph::getDeepVertexByIndex(unsigned int index) {
    for (auto vertex: this->deepVertices) {
        if (vertex->getIndex() == index) {
            return vertex;
        }
    }
    return nullptr;
}

void data_structures::Graph::addEdge(std::shared_ptr<data_structures::Edge> edge) {
    ICrossoverable crossoverable(edge->getMutationChance(), edge->isDominant(),
                                 edge->getMaxChildren());

    this->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(), edge->getOutput()->getType(),
                  edge->getOutput()->getIndex(), edge->getIndex(), edge->getWeight(), edge->getTraverseLimit(),
                  crossoverable);
}

void data_structures::Graph::normalizeEdgeWeights() {
    // find the largest weight
    double maxWeight = 0;
    for (const auto &edge: this->getEdges()) {
        if (fabs(edge->getWeight()) > maxWeight) {
            maxWeight = fabs(edge->getWeight());
        }
    }
    // normalize
    for (const auto &edge: this->getEdges()) {
        double newWeight = edge->getWeight() / maxWeight;
        edge->setWeight(newWeight);
    }
}

void data_structures::Graph::removeEdge(unsigned long position) {
    // get edge
    auto edgeToRemove = this->edges.at(position);

    // remove edge from vertex vectors
    auto vertex = edgeToRemove->getInput();
    // create a new output vertex vector and copy all edges except the one to delete
    //std::vector<std::shared_ptr<data_structures::Edge>> newOutputEdgesVector;
    for (unsigned long i = 0; i < vertex->getOutputEdges().size(); i++) {
        if (vertex->getOutputEdges().at(i)->getIndex() == edgeToRemove->getIndex()) {
            vertex->eraseOutputEdge(i);
            break;
        }
    }


    vertex = edgeToRemove->getOutput();
    // create a new input vertex vector and copy all edges except the one to delete
    //std::vector<std::shared_ptr<data_structures::Edge>> newInputEdgesVector;
    for (unsigned long i = 0; i < vertex->getInputEdges().size(); i++) {
        if (vertex->getInputEdges().at(i)->getIndex() == edgeToRemove->getIndex()) {
            //newInputEdgesVector.push_back(vertex->getInputEdges().at(i));
            vertex->eraseInputEdge(i);
        }
    }
    //vertex->replaceInputEdges(newInputEdgesVector);


    // remove edge from graph vertex
    this->edges.erase(this->edges.begin() + position);

}

std::string data_structures::Graph::toForceGraphJson() {
    auto nodes = nlohmann::json::array();
    for (const auto &vertex: this->inputVertices) {
        nlohmann::json j;
        j["id"] =  vertex->getLabel();
        j["name"] = vertex->getLabel();
        j["color"] = "#000";
        j["pos"] = vertex->getIndex();
        j["group"] = "input";
        nodes.push_back(j);
    }
    for (const auto &vertex: this->deepVertices) {
        nlohmann::json j;
        auto label = std::to_string(vertex->getIndex());
        j["id"] = label;
        j["name"] = label;
        j["color"] = "#0000FF";
        //j["pos"] = vertex->getIndex();
        j["group"] = "deep";
        nodes.push_back(j);
    }
    for (const auto &vertex: this->outputVertices) {
        nlohmann::json j;
        j["id"] = vertex->getLabel();
        j["name"] = vertex->getLabel();
        j["color"] = "#000";
        j["pos"] = vertex->getIndex();
        j["group"] = "output";
        nodes.push_back(j);
    }

    auto links = nlohmann::json::array();
    for (const auto &edge: this->edges) {
        std::string sourceId;
        std::string targetId;
        if (edge->getInput()->getType() == enums::VertexType::Input) {
            auto inputVertex = std::dynamic_pointer_cast<data_structures::InputVertex>(edge->getInput());
            sourceId = inputVertex->getLabel();
        } else {
            sourceId = std::to_string(edge->getInput()->getIndex());
        }

        if (edge->getOutput()->getType() == enums::VertexType::Output) {
            auto outputVertex = std::dynamic_pointer_cast<data_structures::OutputVertex>(edge->getOutput());
            targetId = outputVertex->getLabel();
        } else {
            targetId = std::to_string(edge->getOutput()->getIndex());
        }
        nlohmann::json j;
        j["source"] = sourceId;
        j["target"] = targetId;
        j["value"] = edge->getWeight();
        links.push_back(j);
    }

    nlohmann::json graph = nlohmann::json::object();
    graph["nodes"] = nodes;
    graph["links"] = links;

    return graph.dump();
}
