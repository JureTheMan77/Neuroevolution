//
// Created by jure on 7/13/20.
//

#ifndef NEUROEVOLUTION_GRAPH_H
#define NEUROEVOLUTION_GRAPH_H

#include <vector>
#include <memory>
#include <limits>
#include "InputVertex.h"
#include "OutputVertex.h"
#include "DeepVertex.h"
#include "UniqueVertexQueue.h"
#include "DataInstance.h"
#include "../enums/VertexType.h"
#include "IDeepCloneable.h"

namespace data_structures {
    class Graph : public data_structures::IDeepCloneable<Graph> {
    private:
        std::vector<std::shared_ptr<data_structures::InputVertex>> inputVertices{};
        std::vector<std::shared_ptr<data_structures::OutputVertex>> outputVertices{};
        std::vector<std::shared_ptr<data_structures::DeepVertex>> deepVertices{};

        std::vector<std::shared_ptr<data_structures::Edge>> edges{};

        unsigned int largestDeepVertexIndex{0};
        unsigned int largestEdgeIndex{0};
        // unsigned int const UINT_MAX = std::numeric_limits<unsigned int>::max();

        /**
         * First step of adding an edge.
         */
        bool edgeBeforeAdd(enums::VertexType inputVertexType, unsigned int inputVertexIndex,
                           enums::VertexType outputVertexType, unsigned int outputVertexIndex,
                           unsigned int traversalLimit, std::shared_ptr<data_structures::Vertex> &input,
                           std::shared_ptr<data_structures::Vertex> &output,
                           std::shared_ptr<data_structures::Edge> &commonEdge);

        /**
         * Last step of adding an edge.
         */
        static void addEdgeAfterAdd(const enums::VertexType &inputVertexType, const enums::VertexType &outputVertexType,
                                    const std::shared_ptr<data_structures::Vertex> &input,
                                    const std::shared_ptr<data_structures::Vertex> &output,
                                    const std::shared_ptr<data_structures::Edge> &newEdge);

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

        /**
         * Create a new graph shared pointer.
         * @param inputVertices number of input vertices
         * @param inputLabels labels of input vertices
         * @param deepVertices number of deep vertices
         * @param outputVertices number of output vertices
         * @param outputLabels labels of output vertices
         * @return
         */
        static std::shared_ptr<Graph> createGraph(unsigned int inputVertices, std::vector<std::string> &inputLabels,
                                                  unsigned int deepVertices, unsigned int outputVertices,
                                                  std::vector<std::string> &outputLabels);

        /**
         * Create a copy of this graph.
         * @return graph clone
         */
        std::shared_ptr<Graph> deepClone() override;

        /**
         * Calls destructors on all edges and vertices.
         */
        virtual ~Graph() = default;

        /**
         * Adds new input vertices.
         * @param numberOfInputVertices number of input vertices to add
         * @param labels labels of vertices
         */
        void addInputVertices(unsigned int numberOfInputVertices, std::vector<std::string> &labels);

        /**
         * Add a new input vertex.
         * @param inputVertex vertex to add
         */
        void addInputVertex(const std::shared_ptr<data_structures::InputVertex> &inputVertex);

        /**
         * Get the vector of input vertices.
         * @return input vertex vector
         */
        std::vector<std::shared_ptr<data_structures::InputVertex>> getInputVertices() const;

        /**
         * Adds new output vertices.
         * @param numberOfOutputVertices number of output vertices to add
         * @param labels labels of vertices
         */
        void addOutputVertices(unsigned int numberOfOutputVertices, std::vector<std::string> &labels);

        /**
         * Get the vector of output vertices.
         * @param outputVertex vertex to add
         */
        void addOutputVertex(const std::shared_ptr<data_structures::OutputVertex> &outputVertex);

        /**
         * Get the vector of output vertices.
         * @return output vertex vector
         */
        std::vector<std::shared_ptr<data_structures::OutputVertex>> getOutputVertices() const;

        /**
         * Adds new deep vertices.
         * @param numberOfDeepVertices number of deep vertices to add
         */
        void addDeepVertices(unsigned int numberOfDeepVertices);

        /**
         * Add a new deep vertex.
         * @param deepVertex vertex to add
         */
        void addDeepVertex(const std::shared_ptr<data_structures::DeepVertex> &deepVertex);

        /**
         * Adds an edge between the input and output vertex.
         * @param inputVertexType type of the input vertex
         * @param inputVertexIndex index of the input vertex to connect from
         * @param outputVertexType type of the output vertex
         * @param outputVertexIndex index of the input vertex to connect to
         * @param index edge index
         * @param weight weight of the edge
         * @param traversalLimit how many times this edge can be traversed, minimum value is 1
         */
        void
        addEdge(enums::VertexType inputVertexType, unsigned int inputVertexIndex, enums::VertexType outputVertexType,
                unsigned int outputVertexIndex, unsigned int index, double weight, unsigned int traversalLimit);

        /**
         * Adds an edge between the input and output vertex with crossoverable data.
         * @param inputVertexType type of the input vertex
         * @param inputVertexIndex index of the input vertex to connect from
         * @param outputVertexType type of the output vertex
         * @param outputVertexIndex index of the input vertex to connect to
         * @param index edge index
         * @param weight weight of the edge
         * @param traversalLimit how many times this edge can be traversed, minimum value is 1
         * @param crossoverableData data
         * @return
         */
        std::shared_ptr<data_structures::Edge>
        addEdge(enums::VertexType inputVertexType, unsigned int inputVertexIndex, enums::VertexType outputVertexType,
                unsigned int outputVertexIndex, unsigned int index, double weight, unsigned int traversalLimit,
                const data_structures::ICrossoverable &crossoverableData);

        /**
         * Add a new edge.
         * @param edge edge to add
         */
        void addEdge(const std::shared_ptr<data_structures::Edge> &edge);

        /**B
         * Propagate the values of the @param dataInstance through the graph.
         * @param dataInstance values to propagate
         */
        void traverse(const std::shared_ptr<data_structures::DataInstance> &dataInstance);

        /**
         * Get the index of the output vertex with the largest value.
         * @return index of the vertex
         */
        unsigned int getLargestOutputValueIndex() const;

        /**
         * Resets the values of all vertices and edges.
         */
        void reset();

        /**
         * Get the vector of deep vertices.
         * @return deep vertex vector
         */
        std::vector<std::shared_ptr<data_structures::DeepVertex>> getDeepVertices() const;

        /**
         * Get a vertex with the correct index.
         * @param index index to search for
         * @return the vertex
         */
        std::shared_ptr<data_structures::DeepVertex> getDeepVertexByIndex(unsigned int index) const;

        /**
         * Get a vector od all graph's edges.
         * @return edge vector
         */
        std::vector<std::shared_ptr<data_structures::Edge>> getEdges() const;

        /**
         * Find an edge with the supplied index and input and output vertex types.
         * @param inputVertexIndex index of the input vertex
         * @param inputVertexType type of the input vertex
         * @param outputVertexIndex index of the output vertex
         * @param outputVertexType type of the output vertex
         * @return edge
         */
        std::shared_ptr<data_structures::Edge>
        getEdgeByIndexAndType(unsigned int inputVertexIndex, enums::VertexType inputVertexType,
                              unsigned int outputVertexIndex, enums::VertexType outputVertexType) const;

        /**
         * Assigns proper index values to vertices and edges with index UINT_MAX.
         */
        void fixIndices();

        /**
         * Get a string with useful information about this object.
         * @param technical format for https://csacademy.com/app/graph_editor/
         * @return this object's information
         */
        std::string toString(bool technical) const;

        /**
         * Normalize edge values from -1 to 1.
         */
        void normalizeEdgeWeights() const;

        /**
         * Remove an edge from this graph.
         * @param position position in edge vector
         */
        void removeEdge(unsigned long position);

        /**
         * Create a string for force graph rendering.
         * @return
         */
        std::string toForceGraphJson() const;
    };
}


#endif //NEUROEVOLUTION_GRAPH_H
