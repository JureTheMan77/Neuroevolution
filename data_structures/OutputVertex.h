//
// Created by jure on 4. 02. 21.
//

#ifndef NEUROEVOLUTION_OUTPUTVERTEX_H
#define NEUROEVOLUTION_OUTPUTVERTEX_H

#include "Vertex.h"

namespace data_structures {
    class OutputVertex : public data_structures::Vertex, public data_structures::IDeepCloneable<OutputVertex> {
    private:
        std::string label;

    public:
        /**
         * Creates a new output vertex object.
         * @param index index of the vertex
         * @param label label of the vertex
         */
        OutputVertex(unsigned int index, std::string label) : Vertex(index), label(std::move(label)) {}

        /**
         * Create a new output vertex shared pointer.
         * @param index index of the vertex
         * @param label label of the vertex
         * @return a new shared pointer
         */
        static std::shared_ptr<OutputVertex> createOutputVertex(unsigned int index, std::string &label);

        /**
         * Clone this object. Value is not preserved.
         * @return
         */
        std::shared_ptr<OutputVertex> deepClone() override;

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
         * @param technical format for https://csacademy.com/app/graph_editor/
         * @return this object's information
         */
        std::string toString(bool technical) override;

        /**
         * Get type of this vertex. For easier casting.
         * @return always enums::VertexType::OutputVertex
         */
        enums::VertexType getType() override;

        /**
         * Get this vertex's label.
         * @return label
         */
        [[nodiscard]] const std::string &getLabel() const;
    };
}


#endif //NEUROEVOLUTION_OUTPUTVERTEX_H
