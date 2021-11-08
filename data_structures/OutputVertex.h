//
// Created by jure on 4. 02. 21.
//

#ifndef NEUROEVOLUTION_OUTPUTVERTEX_H
#define NEUROEVOLUTION_OUTPUTVERTEX_H

#include "Vertex.h"

namespace data_structures {
    class OutputVertex : public data_structures::Vertex {
    private:
        std::string label;

    public:
        /**
         * Creates a new output vertex object.
         * @param index index of the vertex
         * @param label label of the vertex
         */
        OutputVertex(unsigned int index, std::string label) :
                Vertex(index, false, 0, 0), label(std::move(label)) {}

        static std::shared_ptr<OutputVertex> createOutputVertex(unsigned int index, std::string label);

        /**
         * Destroys this object.
         */
        ~OutputVertex() override = default;

        /**
         * Output vertices cannot have output edges so this method does nothing.
         * @param edge edge to not add
         */
        void addOutputEdge(const std::shared_ptr<data_structures::Edge> &edge) override;

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        std::string toString() override;
    };
}


#endif //NEUROEVOLUTION_OUTPUTVERTEX_H
