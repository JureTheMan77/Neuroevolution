//
// Created by jure on 3. 02. 21.
//

#ifndef NEUROEVOLUTION_INPUTVERTEX_H
#define NEUROEVOLUTION_INPUTVERTEX_H


#include "Vertex.h"

namespace data_structures {
    class InputVertex : public data_structures::Vertex, public data_structures::IDeepCloneable<InputVertex> {
    private:
        std::string label{};

    public:
        /**
         * Creates a new input vertex object.
         * @param index index of the vertex
         * @param label label of the vertex
         */
        InputVertex(unsigned int index, std::string labelArg) : Vertex(index), label(std::move(labelArg)) {}

        /**
         * Create a new input vertex shared pointer.
         * @param index index of the vertex
         * @param labelArg label of the vertex
         * @return shared pointer
         */
        static std::shared_ptr<InputVertex> createInputVertex(unsigned int index, std::string &labelArg);

        /**
         * Clone this object. Vertex value and edges are not kept.
         * @return
         */
        std::shared_ptr<InputVertex> deepClone() override;

        /**
         * Destroys this object.
         */
        ~InputVertex() override = default;

        /**
         * Input vertex cannot have input edges so this method does nothing.
         * @param edge edge to not add
         */
        void addInputEdge(const std::shared_ptr<data_structures::Edge> &edge) override;

        /**
         * Always returns true.
         * @return true
         */
        bool allInputEdgesTraversedOnce() override;

        /**
         * Get a string with useful information about this object.
         * @param technical format for https://csacademy.com/app/graph_editor/
         * @return this object's information
         */
        std::string toString(bool technical) override;

        /**
         * Get type of this vertex. For easier casting.
         * @return always enums::VertexType::InputVertex
         */
        enums::VertexType getType() override;

        /**
         * Get the label of this vertex.
         * @return label
         */
        [[nodiscard]] const std::string &getLabel() const;
    };
}


#endif //NEUROEVOLUTION_INPUTVERTEX_H
