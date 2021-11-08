//
// Created by jure on 11. 08. 21.
//

#ifndef NEUROEVOLUTION_INPUTEDGES_H
#define NEUROEVOLUTION_INPUTEDGES_H


#include "EdgeGroup.h"

namespace data_structures {
    class InputEdges : private data_structures::EdgeGroup {
    public:
        /**
         * Use base class constructor
         */
        using EdgeGroup::EdgeGroup;

        /**
         * Default simple destructor.
         */
        virtual ~InputEdges() = default;

        /**
         * Attaches an input edge.
         * @param edge the edge to attach
         */
        virtual void addInputEdge(const std::shared_ptr<data_structures::Edge>& edge);

        /**
         * Have all input edges been traversed at least once?
         * @return true if all input edges have been traversed at least once
         */
        virtual bool allInputEdgesTraversed();

        /**
         * Input edges getter.
         * @return the vector of input edges
         */
        virtual std::vector<std::shared_ptr<data_structures::Edge>> getInputEdges();
    };
}


#endif //NEUROEVOLUTION_INPUTEDGES_H
