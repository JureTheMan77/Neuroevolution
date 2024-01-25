//
// Created by jure on 7/13/20.
//

#include "Graph.h"
#include <vector>
#include <string>
#include <sstream>
#include "../enums/EnumUtil.h"
#include "../util/util.h"
#include "../nlohmann/json.hpp"


data_structures::Graph::Graph(unsigned int inputVertices, std::vector<std::string> &inputLabels,
                              unsigned int deepVertices, unsigned int outputVertices,
                              std::vector<std::string> &outputLabels) {
    this->addInputVertices(inputVertices, inputLabels);
    this->addOutputVertices(outputVertices, outputLabels);
    this->addDeepVertices(deepVertices);
    this->calculateNumEdgesPossible();
}

void data_structures::Graph::addInputVertices(unsigned int numberOfInputVertices, std::vector<std::string> &labels) {
    if (numberOfInputVertices != labels.size()) {
        throw std::invalid_argument("The number of input vertices and labels does not match.");
    }
    for (unsigned int i = 0; i < numberOfInputVertices; i++) {
        auto inputVertex = data_structures::InputVertex::createInputVertex(i, labels.at(i));
        this->inputVertices.push_back(inputVertex);
    }
    this->calculateNumEdgesPossible();
}

void data_structures::Graph::addInputVertex(const std::shared_ptr<data_structures::InputVertex> &inputVertex) {
    this->inputVertices.push_back(inputVertex);
    this->calculateNumEdgesPossible();
}

void data_structures::Graph::addOutputVertices(unsigned int numberOfOutputVertices, std::vector<std::string> &labels) {
    for (unsigned int i = 0; i < numberOfOutputVertices; i++) {
        this->outputVertices.push_back(data_structures::OutputVertex::createOutputVertex(i, labels.at(i)));
    }
    this->calculateNumEdgesPossible();
}

void data_structures::Graph::addOutputVertex(const std::shared_ptr<data_structures::OutputVertex> &outputVertex) {
    this->outputVertices.push_back(outputVertex);
    this->calculateNumEdgesPossible();
}

void data_structures::Graph::addDeepVertices(unsigned int numberOfDeepVertices) {
    for (unsigned int i = 0; i < numberOfDeepVertices; i++) {
        // initialise chances as values from 0 to 1
        auto deepVertex = data_structures::DeepVertex::createDeepVertex(i, util::nextBool());
        this->deepVertices.push_back(deepVertex);
        this->largestDeepVertexIndex = i;
    }
    this->calculateNumEdgesPossible();
}

void data_structures::Graph::addDeepVertex(const std::shared_ptr<data_structures::DeepVertex> &deepVertex) {
    // assume that there are not vertices with the same index
    if (this->largestDeepVertexIndex < deepVertex->getIndex() && deepVertex->getIndex() != UINT_MAX) {
        this->largestDeepVertexIndex = deepVertex->getIndex();
    }
    this->deepVertices.push_back(deepVertex);
    this->calculateNumEdgesPossible();
}

void data_structures::Graph::addEdge(enums::VertexType inputVertexType, unsigned int inputVertexIndex,
                                     enums::VertexType outputVertexType, unsigned int outputVertexIndex,
                                     unsigned int index, double weight, unsigned int traversalLimit) {
    std::shared_ptr<data_structures::Vertex> input;
    std::shared_ptr<data_structures::Vertex> output;
    std::shared_ptr<data_structures::Edge> commonEdge;
    // if the result of the  function is false, then don't add the edge
    if (!edgeBeforeAdd(inputVertexType, inputVertexIndex, outputVertexType, outputVertexIndex, traversalLimit, input,
                       output, commonEdge)) {
        return;
    }

    // if a common edge has not been found, create a new one
    // otherwise, return
    if (commonEdge == nullptr) {
        // assume that there are no edges with the same index
        if (this->largestEdgeIndex < index && index != UINT_MAX) {
            this->largestEdgeIndex = index;
        }
        std::shared_ptr<data_structures::Edge> newEdge = data_structures::Edge::createEdge(input, output, index, weight,
                                                                                           traversalLimit,
                                                                                           util::nextBool());
        if(newEdge == nullptr){
            std::cout<<"Graph.cpp:86";
            std::invalid_argument("Edge is null!");
        }
        this->edges.push_back(newEdge);
        addEdgeAfterAdd(inputVertexType, outputVertexType, input, output, newEdge);
    } else {
        //commonEdge->combineWeight(weight);
    }
}

void data_structures::Graph::addEdgeAfterAdd(const enums::VertexType &inputVertexType,
                                             const enums::VertexType &outputVertexType,
                                             const std::shared_ptr<data_structures::Vertex> &input,
                                             const std::shared_ptr<data_structures::Vertex> &output,
                                             const std::shared_ptr<data_structures::Edge> &newEdge) {
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

bool data_structures::Graph::edgeBeforeAdd(enums::VertexType inputVertexType, unsigned int inputVertexIndex,
                                           enums::VertexType outputVertexType, unsigned int outputVertexIndex,
                                           unsigned int traversalLimit, std::shared_ptr<data_structures::Vertex> &input,
                                           std::shared_ptr<data_structures::Vertex> &output,
                                           std::shared_ptr<data_structures::Edge> &commonEdge) {
    commonEdge = nullptr;
    if (traversalLimit < 1) {
        throw std::invalid_argument(
                "Traversal limit has to be at least 1. The supplied value is " + std::to_string(traversalLimit) + ".");
    }

    // already at max edges?
    if(this->edges.size() == this->numEdgesPossible){
        return false;
    }

    switch (inputVertexType) {
        case enums::VertexType::Input:
            input = inputVertices.at(inputVertexIndex);
            break;
        case enums::VertexType::Deep:
            input = this->getDeepVertexByIndex(inputVertexIndex);
            //if (input == nullptr) {
            //    throw std::invalid_argument("[INPUT] Deep vertex with index " + std::to_string(inputVertexIndex) + " doesn't exist.");
            //}
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
            //if (output == nullptr) {
            //    throw std::invalid_argument("[OUTPUT] Deep vertex with index " + std::to_string(outputVertexIndex) + " doesn't exist.");
            //}
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
    // if the result of the function is false, then don't add the edge
    if (!edgeBeforeAdd(inputVertexType, inputVertexIndex, outputVertexType, outputVertexIndex, traversalLimit, input,
                       output, commonEdge)) {
        return nullptr;
    }

    // if a common edge has not been found, create a new one
    // otherwise, return it
    if (commonEdge == nullptr) {
        // assume that there are no edges with the same index
        if (this->largestEdgeIndex < index && index != UINT_MAX) {
            this->largestEdgeIndex = index;
        }
        std::shared_ptr<data_structures::Edge> newEdge = data_structures::Edge::createEdge(input, output, index, weight,
                                                                                           traversalLimit,
                                                                                           crossoverableData.isDominant());
        if(newEdge == nullptr){
            std::cout<<"Graph.cpp:202";
            std::invalid_argument("Edge is null!");
        }
        this->edges.push_back(newEdge);
        addEdgeAfterAdd(inputVertexType, outputVertexType, input, output, newEdge);
        return newEdge;
    } else {
        //commonEdge->combineWeight(weight);
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

    // queue of pending vertices
    data_structures::UniqueVertexQueue pendingVertices;
    // set of unique backup vertices in case the above one is empty
    std::unordered_set<std::shared_ptr<data_structures::Vertex>> backupVertices;

    // add the input vertices to the queue of pending vertices
    // enqueue them only if they can be traversed out from
    for (const auto &v: this->inputVertices) {
        if (v->canTraverseOut()) {
            pendingVertices.enqueue(v);
        }
    }

    // propagate the values through the graph
    // repeat until the set of pending vertices is empty
    //std::vector<std::shared_ptr<data_structures::Edge>> outputEdges;
    this->numOfPropagations = 0;
    while (!pendingVertices.empty()) {
        std::shared_ptr<data_structures::Vertex> vertex = pendingVertices.dequeue();

        // vertex->getOutputEdges() cannot be empty at this stage
        for (const std::shared_ptr<data_structures::Edge> &edge: vertex->getOutputEdges()) {
            if (!edge->isAtTraverseLimit()) {
                // propagate the value to adjacent vertex
                edge->propagateValue();
                this->numOfPropagations += 1;

                // only add the output vertex to the pending vertices set queue if:
                // - its type is "deep"
                // - has output edges that can still be traversed
                // - its input edges have been fully traversed; if only this condition is false, then add the
                // vertex to the set of unique backup vertices (possible cycle)
                if (edge->getOutput()->getType() == enums::VertexType::Deep &&
                    edge->getOutput()->canTraverseOut()) {
                    if (edge->getOutput()->inputEdgeTraversalsRemaining() == 0) {
                        pendingVertices.enqueue(edge->getOutput());
                        // just in case the backup set already contains the same vertex
                        backupVertices.erase(edge->getOutput());
                    } else {
                        backupVertices.insert(edge->getOutput());
                    }
                }
            }
        }

        // continue if there are still vertices pending
        if (!pendingVertices.empty()) {
            continue;
        }

        // if the set of backup vertices is empty here, break the loop as there is no more processing to be done
        if (backupVertices.empty()) {
            break;
        }

        // if the queue of pending vertices is empty, select a pending vertex with the least remaining input traversals
        // and add it to the pending queue
        // this can happen if the graph has cycles or the input edges have too high of a traversal limit
        // if multiple vertices have the same remaining traversals, choose the one with the lowest index
        std::shared_ptr<data_structures::Vertex> vertexToEnqueue = nullptr;
        unsigned int leastTraversalsRemaining = UINT32_MAX;

        for (const std::shared_ptr<data_structures::Vertex> &backupVertex: backupVertices) {
            // the above loop has already checked that these vertices can be traversed out of
            // no additional validation is necessary
            unsigned int traversalsRemaining = backupVertex->inputEdgeTraversalsRemaining();
            if (traversalsRemaining < leastTraversalsRemaining) {
                vertexToEnqueue = backupVertex;
                leastTraversalsRemaining = traversalsRemaining;
            } else if (traversalsRemaining == leastTraversalsRemaining &&
                       backupVertex->getIndex() < vertexToEnqueue->getIndex()) {
                vertexToEnqueue = backupVertex;
            }
        }

        if (vertexToEnqueue != nullptr) {
            pendingVertices.enqueue(vertexToEnqueue);
            backupVertices.erase(vertexToEnqueue);
        }
    }
}

unsigned int data_structures::Graph::getLargestOutputValueIndex() const {
    unsigned int largestIndex = 0;
    double largestValue = -std::numeric_limits<double>::max();
    for (const std::shared_ptr<data_structures::OutputVertex> &v: this->outputVertices) {
        if (v->getActivationValue() > largestValue) {
            largestValue = v->getActivationValue();
            largestIndex = v->getIndex();
        }
    }

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

    for (const auto &e: this->edges) {
        //if(e != nullptr) {
            e->reset();
        //}
    }

    this->numOfPropagations = 0;
}

std::vector<std::shared_ptr<data_structures::DeepVertex>> data_structures::Graph::getDeepVertices() const {
    return this->deepVertices;
}

std::vector<std::shared_ptr<data_structures::InputVertex >> data_structures::Graph::getInputVertices() const {
    return this->inputVertices;
}

std::vector<std::shared_ptr<data_structures::OutputVertex>> data_structures::Graph::getOutputVertices() const {
    return this->outputVertices;
}

std::string data_structures::Graph::toString(bool technical) const {
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
                                    std::vector<std::string> &outputLabels) {
    return std::make_shared<data_structures::Graph>(inputVertices, inputLabels, deepVertices, outputVertices,
                                                    outputLabels);
}

std::shared_ptr<data_structures::Graph> data_structures::Graph::deepClone() {
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
            data_structures::ICrossoverable crossoverable(edge->isDominant());
            newGraph->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(), edge->getOutput()->getType(),
                              edge->getOutput()->getIndex(), edge->getIndex(), edge->getWeight(),
                              edge->getTraverseLimit(), crossoverable);
        }
        newGraph->calculateNumEdgesPossible();
    }

    // recreate edges of deep vertices
    for (const std::shared_ptr<data_structures::DeepVertex> &deepVertex: this->deepVertices) {
        for (const std::shared_ptr<data_structures::Edge> &edge: deepVertex->getOutputEdges()) {
            data_structures::ICrossoverable crossoverable(edge->isDominant());
            newGraph->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(), edge->getOutput()->getType(),
                              edge->getOutput()->getIndex(), edge->getIndex(), edge->getWeight(),
                              edge->getTraverseLimit(), crossoverable);
        }
    }

    // output vertices have no output edges

    // pendingVertices stays empty
    return newGraph;
}

void data_structures::Graph::fixIndices() {
    for (const auto &deepVertex: this->deepVertices) {
        if (deepVertex->getIndex() == UINT_MAX) {
            this->largestDeepVertexIndex += 1;
            deepVertex->setIndex(this->largestDeepVertexIndex);
        }
    }
    for (const auto &edge: this->edges) {
        if (/*edge != nullptr &&*/ edge->getIndex() == UINT_MAX) {
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
                                              unsigned int outputVertexIndex,
                                              enums::VertexType outputVertexType) const {
    for (const auto &edge: this->edges) {
        if (edge->getInput()->getIndex() == inputVertexIndex && edge->getInput()->getType() == inputVertexType &&
            edge->getOutput()->getIndex() == outputVertexIndex && edge->getOutput()->getType() == outputVertexType) {
            return edge;
        }
    }
    return nullptr;
}

std::shared_ptr<data_structures::DeepVertex> data_structures::Graph::getDeepVertexByIndex(unsigned int index) const {
    for (auto vertex: this->deepVertices) {
        if (vertex->getIndex() == index) {
            return vertex;
        }
    }
    return nullptr;
}

void data_structures::Graph::addEdge(const std::shared_ptr<data_structures::Edge> &edge) {
    ICrossoverable crossoverable(edge->isDominant());

    this->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(), edge->getOutput()->getType(),
                  edge->getOutput()->getIndex(), edge->getIndex(), edge->getWeight(), edge->getTraverseLimit(),
                  crossoverable);
}

void data_structures::Graph::normalizeEdgeWeights() const {
    // find the largest weight
    double maxWeight = 0;
    for (const auto &edge: this->getEdges()) {
        //if(edge != nullptr) {
            double absoluteWeight = std::fabs(edge->getWeight());
            if (absoluteWeight > maxWeight) {
                maxWeight = absoluteWeight;
            }
        //}
    }
    // normalize
    for (const auto &edge: this->getEdges()) {
        //if (edge != nullptr) {
            double newWeight = edge->getWeight() / maxWeight;
            edge->setWeight(newWeight);
        //}
    }
}

void data_structures::Graph::removeEdge(unsigned long position) {
    // get edge
    auto edgeToRemove = this->edges.at(position);

    // remove edge from vertex vectors
    auto vertex = edgeToRemove->getInput();

    for (unsigned long i = 0; i < vertex->getOutputEdges().size(); i++) {
        if (vertex->getOutputEdges().at(i)->getIndex() == edgeToRemove->getIndex()) {
            vertex->eraseOutputEdge(i);
            break;
        }
    }

    vertex = edgeToRemove->getOutput();
    for (unsigned long i = 0; i < vertex->getInputEdges().size(); i++) {
        if (vertex->getInputEdges().at(i)->getIndex() == edgeToRemove->getIndex()) {
            vertex->eraseInputEdge(i);
        }
    }

    // remove edge from graph vertex
    this->edges.erase(this->edges.begin() + (long) position);

}

std::string data_structures::Graph::toForceGraphJson() const {
    auto nodes = nlohmann::json::array();
    for (const auto &vertex: this->inputVertices) {
        nlohmann::json j;
        j["id"] = vertex->getLabel();
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
        j["value"] = std::to_string(edge->getWeight()) + "|" + std::to_string(edge->getTraverseLimit());
        links.push_back(j);
    }

    nlohmann::json graph = nlohmann::json::object();
    graph["nodes"] = nodes;
    graph["links"] = links;

    return graph.dump();
}

void data_structures::Graph::calculateNumEdgesPossible() {
    unsigned long a = this->inputVertices.size();
    unsigned long b = this->deepVertices.size();
    unsigned long c = this->outputVertices.size();
    // a*c+a*b+b*c+b*b
    // - inputs directly to outputs: a*c
    // - inputs to deep: a*b
    // - deep to outputs: b*c
    // - deep to all other deep, including themselves: b*b
    // shorter form: (a+b)(b+c)

    this->numEdgesPossible = (a + b) * (b + c);
}

unsigned long data_structures::Graph::getNumEdgesPossible() const {
    return this->numEdgesPossible;
}

unsigned int data_structures::Graph::getNumOfPropagations() const {
    return this->numOfPropagations;
}
