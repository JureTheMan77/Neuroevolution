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


data_structures::Graph::Graph(unsigned int inputVertices, std::vector<std::string> &inputLabels,
                              unsigned int deepVertices,
                              unsigned int outputVertices, std::vector<std::string> &outputLabels) {
    this->addInputVertices(inputVertices, inputLabels);
    this->addOutputVertices(outputVertices, outputLabels);
    this->addDeepVertices(deepVertices);
}

void data_structures::Graph::addInputVertices(unsigned int numberOfInputVertices, std::vector<std::string> &labels) {
    if (numberOfInputVertices != labels.size()) {
        throw std::invalid_argument("The number of input vertices and labels does not match.");
    }
    for (unsigned int i = 0; i < numberOfInputVertices; i++) {
        this->inputVertices.push_back(
                data_structures::InputVertex::createInputVertex(i, util::nextBool(), util::nextDouble(),
                                                                util::nextInt(0, 2), labels.at(i)));
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

void data_structures::Graph::addDeepVertices(unsigned int numberOfDeepVertices) {
    for (unsigned int i = 0; i < numberOfDeepVertices; i++) {
        // initialise chances as values from 0 to 1
        this->deepVertices.push_back(
                data_structures::DeepVertex::createDeepVertex(i, util::nextBool(), util::nextDouble(),
                                                              util::nextInt(0, 2)));
    }
}

void data_structures::Graph::addDeepVertex(const std::shared_ptr<data_structures::DeepVertex> &deepVertex) {
    this->deepVertices.push_back(deepVertex);
}

void data_structures::Graph::addEdge(enums::VertexType inputVertexType, unsigned int inputVertexIndex,
                                     enums::VertexType outputVertexType, unsigned int outputVertexIndex, double weight,
                                     unsigned int traversalLimit) {
    if (traversalLimit < 1) {
        throw std::invalid_argument(
                "Traversal limit has to be at least 1. The supplied value is " + std::to_string(traversalLimit) + ".");
    }

    std::shared_ptr<data_structures::Vertex> input;
    std::shared_ptr<data_structures::Vertex> output;

    switch (inputVertexType) {
        case enums::VertexType::Input:
            input = this->inputVertices.at(inputVertexIndex);
            break;
        case enums::VertexType::Deep:
            input = this->deepVertices.at(inputVertexIndex);
            break;
        default:
            throw std::invalid_argument("Vertex of type " + enums::EnumUtil::VertexTypeToString(inputVertexType) +
                                        " cannot be used as an input vertex.");
    }

    switch (outputVertexType) {
        case enums::VertexType::Output:
            output = this->outputVertices.at(outputVertexIndex);
            break;
        case enums::VertexType::Deep:
            output = this->deepVertices.at(outputVertexIndex);
            break;
        default:
            throw std::invalid_argument("Vertex of type " + enums::EnumUtil::VertexTypeToString(outputVertexType) +
                                        " cannot be used as an output vertex.");
    }

    // check if input and output are already connected
    std::shared_ptr<data_structures::Edge> commonEdge = nullptr;
    for (const std::shared_ptr<data_structures::Edge> &e: input->getOutputEdges()) {
        if (e->getOutput() == output) {
            commonEdge = e;
            break;
        }
    }

    // if a common edge has not been found, create a new one
    // otherwise, combine the weights
    if (commonEdge == nullptr) {
        std::shared_ptr<data_structures::Edge> newEdge = data_structures::Edge::createEdge(input, output, weight,
                                                                                           traversalLimit);
        this->edges.push_back(newEdge);

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
    } else {
        commonEdge->combineWeight(weight);
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
    for (const std::shared_ptr<data_structures::InputVertex> &v: this->inputVertices) {
        this->pendingVertices.enqueue(v);
    }

    // propagate the values through the graph
    // repeat until the set of pending vertices is empty
    std::shared_ptr<data_structures::Vertex> vertex;
    std::vector<std::shared_ptr<data_structures::Edge>> outputEdges;
    while (!this->pendingVertices.empty()) {
        vertex = this->pendingVertices.dequeue();

        // if vertex is an OutputVertex, then getOutputEdges() will return an empty array
        if (vertex->getOutputEdges().empty()) {
            continue;
        }
        for (const std::shared_ptr<data_structures::Edge> &edge: vertex->getOutputEdges()) {
            if (!edge->isAtTaversalLimit()) {
                // propagate the values to adjacent vertices and add the output vertex to the pending vertices set
                edge->propagateValue();
                if (edge->getOutput()->allInputEdgesTraversed()) {
                    this->pendingVertices.enqueue(edge->getOutput());
                }
            }
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

std::string data_structures::Graph::toString() {
    std::ostringstream result;
    result << "Input vertices: " << std::to_string(this->inputVertices.size()) << std::endl;
    result << "Output vertices: " << std::to_string(this->outputVertices.size()) << std::endl;
    result << "Deep vertices: " << std::to_string(this->deepVertices.size()) << std::endl;
    result << "Edges: " << std::to_string(this->edges.size());

    for (const std::shared_ptr<data_structures::Edge> &edge: this->edges) {
        result << std::endl << "\t" << edge->toString();
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
    // recreate edges
    for (const std::shared_ptr<data_structures::InputVertex> &inputVertex: this->inputVertices) {
        // input vertices have no input edges
        for (const std::shared_ptr<data_structures::Edge>& edge: inputVertex->getOutputEdges()) {
            newGraph->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(),
                              edge->getOutput()->getType(), edge->getOutput()->getIndex(),
                              edge->getWeight(), edge->getTraverseLimit());
        }
    }


    // clone output vertices
    for (const std::shared_ptr<data_structures::OutputVertex> &outputVertex: this->outputVertices) {
        newGraph->addOutputVertex(outputVertex->deepClone());
    }
    // recreate edges
    for (const std::shared_ptr<data_structures::OutputVertex> &outputVertex: this->outputVertices) {
        // output vertices have no output edges
        for (const std::shared_ptr<data_structures::Edge>& edge: outputVertex->getInputEdges()) {
            newGraph->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(),
                              edge->getOutput()->getType(), edge->getOutput()->getIndex(),
                              edge->getWeight(), edge->getTraverseLimit());
        }
    }


    // clone deep vertices
    for (const std::shared_ptr<data_structures::DeepVertex> &deepVertex: this->deepVertices) {
        newGraph->addDeepVertex(deepVertex->deepClone());
    }
    // recreate edges
    for (const std::shared_ptr<data_structures::DeepVertex> &deepVertex: this->deepVertices) {
        for (const std::shared_ptr<data_structures::Edge>& edge: deepVertex->getInputEdges()) {
            newGraph->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(),
                              edge->getOutput()->getType(), edge->getOutput()->getIndex(),
                              edge->getWeight(), edge->getTraverseLimit());
        }
        for (const std::shared_ptr<data_structures::Edge>& edge: deepVertex->getOutputEdges()) {
            newGraph->addEdge(edge->getInput()->getType(), edge->getInput()->getIndex(),
                              edge->getOutput()->getType(), edge->getOutput()->getIndex(),
                              edge->getWeight(), edge->getTraverseLimit());
        }
    }

    // pendingVertices stays empty
    return newGraph;
}
