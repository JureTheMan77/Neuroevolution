//
// Created by jure on 18. 08. 21.
//

#ifndef NEUROEVOLUTION_OUTPUTEDGES_H
#define NEUROEVOLUTION_OUTPUTEDGES_H

#include "EdgeGroup.h"

namespace data_structures {
    class OutputEdges : private data_structures::EdgeGroup {
    public:
        /**
         * Use base class constructor
         */
        using EdgeGroup::EdgeGroup;

        /**
         * Default constructor.
         */
        virtual ~OutputEdges() = default;

        /**
         * Attaches an output edge.
         * @param edge the edge to attach
         */
        virtual void addOutputEdge(const std::shared_ptr<data_structures::Edge>& edge);

        /**
         * Output edges getter.
         * @return the vector of input edges
         */
        virtual std::vector<std::shared_ptr<data_structures::Edge>> getOutputEdges();
    };
}


#endif //NEUROEVOLUTION_OUTPUTEDGES_H
