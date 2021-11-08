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

namespace data_structures {
    class Vertex : public data_structures::InputEdges, public data_structures::OutputEdges {
    private:
        unsigned int index;
        double value{};

        bool dominant{};
        double chanceToGetDominated{};
        unsigned int maxChildren{};

        /**
         * Creates a new vertex object.
         * @param index index of the vertex
         */
        explicit Vertex(unsigned int index) : data_structures::InputEdges(), data_structures::OutputEdges(), index(index) {}

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
        Vertex(unsigned int index, bool dominant, double chanceToGetDominated, unsigned int maxChildren) :
                index(index), dominant(dominant), chanceToGetDominated(chanceToGetDominated),
                maxChildren(maxChildren) {};

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

        /**
         * Sets the stored value to 0.
         */
        virtual void reset();

        /**
         * Returns the index of this vertex.
         * @return index
         */
        virtual unsigned int getIndex();

        /**
         * Returns true if this vertex is dominant.
         * @return dominance
         */
        bool isDominant();

        /**
         * Returns the chance of this vertex to be dominated
         * @return domination chance
         */
        double getChanceToGetDominated();

        /**
         * Returns the maximum number of children this vertex can produce.
         * @return max children
         */
        unsigned int getMaxChildren();

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        virtual std::string toString();
    };
}

#endif //NEUROEVOLUTION_VERTEX_H
