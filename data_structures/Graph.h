//
// Created by jure on 7/13/20.
//

#ifndef NEUROEVOLUTION_GRAPH_H
#define NEUROEVOLUTION_GRAPH_H

#include <vector>
#include <memory>
#include "InputVertex.h"
#include "OutputVertex.h"
#include "DeepVertex.h"
#include "UniqueVertexQueue.h"
#include "DataInstance.h"
#include "../enums/VertexType.h"

namespace data_structures {
    class Graph {
    private:
        std::vector<std::shared_ptr<data_structures::InputVertex>> inputVertices{};
        std::vector<std::shared_ptr<data_structures::OutputVertex>> outputVertices{};
        std::vector<std::shared_ptr<data_structures::DeepVertex>> deepVertices{};
        data_structures::UniqueVertexQueue pendingVertices{};

        std::vector<std::shared_ptr<data_structures::Edge>> edges{};

    public:
        /**
         * Initializes inputVertices, outputVertices, deepVertices and edges vectors.
         */
        Graph() = default;

        /**
         * Creates a new Graph.
         * @param inputVertices number of input vertices
         * @param inputLabels labels of input vertices
         * @param deepVertices number of deep vertices
         * @param outputVertices number of output vertices
         * @param outputLabels labels of output vertices
         */
        Graph(unsigned int inputVertices, std::vector<std::string> &inputLabels, unsigned int deepVertices,
              unsigned int outputVertices, std::vector<std::string> &outputLabels);

        static std::shared_ptr<Graph> createGraph(unsigned int inputVertices, std::vector<std::string> &inputLabels,
                                                  unsigned int deepVertices, unsigned int outputVertices,
                                                  std::vector<std::string> &outputLabels);

        /**
         * Calls destructors on all edges and vertices.
         */
        ~Graph() = default;

        /**
         * Adds new input vertices.
         * @param numberOfInputVertices number of input vertices to add
         * @param labels labels of vertices
         */
        void addInputVertices(unsigned int numberOfInputVertices, std::vector<std::string> &labels);

        void addInputVertex(const std::shared_ptr<data_structures::InputVertex> &inputVertex);

        std::vector<std::shared_ptr<data_structures::InputVertex>> getInputVertices();

        /**
         * Adds new output vertices.
         * @param numberOfOutputVertices number of output vertices to add
         * @param labels labels of vertices
         */
        void addOutputVertices(unsigned int numberOfOutputVertices, std::vector<std::string> &labels);

        std::vector<std::shared_ptr<data_structures::OutputVertex>> getOutputVertices();

        /**
         * Adds new deep vertices.
         * @param numberOfDeepVertices number of deep vertices to add
         */
        void addDeepVertices(unsigned int numberOfDeepVertices);

        /**
         * Adds an edge between the input and output vertex.
         * @param inputVertexType type of the input vertex
         * @param inputVertexIndex index of the input vertex to connect from
         * @param outputVertexType type of the output vertex
         * @param outputVertexIndex index of the input vertex to connect to
         * @param weight weight of the edge
         * @param traversalLimit how many times this edge can be traversed, minimum value is 1
         */
        void
        addEdge(enums::VertexType inputVertexType, unsigned int inputVertexIndex, enums::VertexType outputVertexType,
                unsigned int outputVertexIndex, double weight, unsigned int traversalLimit);

        /**
         * Propagate the values of the @param dataInstance through the graph.
         * @param dataInstance values to propagate
         */
        void traverse(const std::shared_ptr<data_structures::DataInstance> &dataInstance);

        /**
         * Get the index of the output vertex with the largest value.
         * @return index of the vertex
         */
        unsigned int getLargestOutputValueIndex();

        /**
         * Resets the values of all vertices and edges.
         */
        void reset();

        std::vector<std::shared_ptr<data_structures::DeepVertex>> getDeepVertices();

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        std::string toString();
    };
}


#endif //NEUROEVOLUTION_GRAPH_H
