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
         * @param chanceToGetDominated chance to get dominated
         * @param maxChildren maximum amount of children this vertex can produce
         */
        explicit DeepVertex(unsigned int index, bool dominant, double chanceToGetDominated, double mutationChance,
                            unsigned int maxChildren) :
                Vertex(index),
                ICrossoverable(mutationChance, dominant, chanceToGetDominated, maxChildren) {};

        /**
         * Creates a new deep vertex object.
         * @param index vertex index
         * @param dominant is this vertex dominant
         * @param chanceToGetDominated chance to get dominated
         * @param maxChildren maximum amount of children this vertex can produce
         * @return smart pointer
         */
        static std::shared_ptr<data_structures::DeepVertex>
        createDeepVertex(unsigned int index, bool dominant, double chanceToGetDominated, double mutationChanceArg,
                         unsigned int maxChildren);

        std::shared_ptr<DeepVertex> deepClone() override;

        /**
         * Deletes the input and output vertex vectors.
         */
        ~DeepVertex() override = default;

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        std::string toString(bool technical) override;

        enums::VertexType getType() override;

        bool isFlaggedForDeletion() const;

        void setFlaggedForDeletion(bool flaggedForDeletion);

        bool allInputEdgesFlaggedForDeletion();

        bool allOutputEdgesFlaggedForDeletion();
    };
}


#endif //NEUROEVOLUTION_DEEPVERTEX_H
