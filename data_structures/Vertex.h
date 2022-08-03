//
// Created by jure on 7/13/20.
//

#ifndef NEUROEVOLUTION_VERTEX_H
#define NEUROEVOLUTION_VERTEX_H

#include <string>
#include "Vertex.fwd.h"
#include "Edge.fwd.h"
#include "InputEdges.h"
#include "OutputEdges.h"
#include "IDeepCloneable.h"
#include "IVertexType.h"

namespace data_structures {
    class Vertex
            : public data_structures::InputEdges,
              public data_structures::OutputEdges,
              public data_structures::IVertexType {
    private:
        unsigned int index;
        double value{0};

        bool visited{};
    public:
        /**
         * Creates a new vertex object.
         * @param index
         * @return smart pointer
         */
        // static std::unique_ptr<Vertex> createVertex(unsigned int index);

        /**
         * Create a nev vertex object.
         * @param index vertex index
         * @param dominant is this vertex dominant
         * @param chanceToGetDominated chance to get dominated
         * @param maxChildren maximum amount of children this vertex can produce
         */
        explicit Vertex(unsigned int index) : index(index){};

        /**
         * Copy constructor.
         * @param v object to copy
         */
        Vertex(Vertex const &v) = default;

        /**
         * Destroys this object.
         */
        ~Vertex() override = default;

        /**
         * Combines @param argValue with the stored value. Default behaviour is addition.
         * @param argValue value to combine with the stored value
         */
        virtual void combineValue(double argValue);

        /**
         * Returns the stored value
         * @return stored value
         */
        double getValue();

        void setValue(double value);

        /**
         * Sets the stored value to 0.
         */
        virtual void reset();

        /**
         * Returns the index of this vertex.
         * @return index
         */
        virtual unsigned int getIndex();

        void setIndex(unsigned int index);

        bool isVisited() const;

        void setVisited(bool visited);

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        virtual std::string toString(bool technical);

        enums::VertexType getType() override = 0;
    };
}

#endif //NEUROEVOLUTION_VERTEX_H
