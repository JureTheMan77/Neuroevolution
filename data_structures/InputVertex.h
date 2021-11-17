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
        InputVertex(unsigned int index, bool dominant, double chanceToGetDominated, unsigned int maxChildren,
                    std::string labelArg) :
                Vertex(index, dominant, chanceToGetDominated, maxChildren), label(std::move(labelArg)) {}

        static std::shared_ptr<InputVertex>
        createInputVertex(unsigned int index, bool dominant, double chanceToGetDominated, unsigned int maxChildren,
                          std::string &labelArg);

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
         * @return
         */
        bool allInputEdgesTraversed() override;

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        std::string toString() override;

        enums::VertexType getType() override;
    };
}


#endif //NEUROEVOLUTION_INPUTVERTEX_H
