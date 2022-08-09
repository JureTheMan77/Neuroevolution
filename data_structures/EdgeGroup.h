//
// Created by jure on 11. 08. 21.
//

#ifndef NEUROEVOLUTION_EDGEGROUP_H
#define NEUROEVOLUTION_EDGEGROUP_H

#include <vector>
#include "Edge.h"

namespace data_structures {
    class EdgeGroup {
    protected:
        std::vector<std::shared_ptr<data_structures::Edge>> edges{};
    public:
        /**
         * Default simple constructor.
         */
        explicit EdgeGroup() = default;

        /**
         * Default simple destructor.
         */
        ~EdgeGroup() = default;

        /**
         * Attaches an edge.
         * @param edge the edge to attach
         */
        void addEdge(const std::shared_ptr<data_structures::Edge> &edge);

        /**
         * Edges getter.
         * @return the vector of input edges
         */
        std::vector<std::shared_ptr<data_structures::Edge>> getEdges();

        void replaceEdges(std::vector<std::shared_ptr<data_structures::Edge>> const &newEdges);

        void clearEdges();

        void eraseEdge(unsigned long position);
    };
}


#endif //NEUROEVOLUTION_EDGEGROUP_H
