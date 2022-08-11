//
// Created by jure on 4. 02. 21.
//

#ifndef NEUROEVOLUTION_DEEPVERTEX_H
#define NEUROEVOLUTION_DEEPVERTEX_H

#include "Vertex.h"
#include "ICrossoverable.h"

namespace data_structures {
    class DeepVertex
            : public data_structures::Vertex,
              public data_structures::IDeepCloneable<DeepVertex>,
              public data_structures::ICrossoverable {
    private:
        bool flaggedForDeletion = false;
    public:
        /**
         * Creates a new deep vertex object without edges, calls the parent constructor.
         * @param index vertex index
         * @param dominant is this vertex dominant
         */
        explicit DeepVertex(unsigned int index, bool dominant) :
                Vertex(index),
                ICrossoverable(dominant) {};

        /**
         * Creates a new deep vertex object.
         * @param index vertex index
         * @param dominant is this vertex dominant
         * @return smart pointer
         */
        static std::shared_ptr<data_structures::DeepVertex>
        createDeepVertex(unsigned int index, bool dominant);

        /**
         * Clones this object. Does not copy the vertex value.
         * @return a new shared pointer with the same field values
         */
        std::shared_ptr<DeepVertex> deepClone() override;

        /**
         * Nothing to see here.
         */
        ~DeepVertex() override = default;

        /**
         * Get a string with useful information about this object.
         * @param technical format for https://csacademy.com/app/graph_editor/
         * @return this object's information
         */
        std::string toString(bool technical) override;

        /**
         * Get vertex type.
         * @return enums::VertexType::DeepVertex
         */
        enums::VertexType getType() override;

        /**
         * Needed for agent minimization step.
         * @return should this vertex be kept
         */
        [[nodiscard]] bool isFlaggedForDeletion() const;

        /**
         * Needed for agent minimization step.
         * @param flaggedForDeletionArg should this vertex be kept
         */
        void setFlaggedForDeletion(bool flaggedForDeletionArg);

        /**
         * Needed for agent minimization step.
         * @return should any input edges be kept
         */
        bool allInputEdgesFlaggedForDeletion();

        /**
         * Needed for agent minimization step.
         * @return should any output edges be kept
         */
        bool allOutputEdgesFlaggedForDeletion();

        /**
         * Use the leaky RELU activation function.
         * @param argValue value to add
         */
        void combineValue(double argValue) override;


    };
}

#endif //NEUROEVOLUTION_DEEPVERTEX_H
