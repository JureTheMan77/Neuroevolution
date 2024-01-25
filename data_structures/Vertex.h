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
#include "../enums/VertexType.h"

namespace data_structures {
    class Vertex
            : public data_structures::InputEdges,
              public data_structures::OutputEdges{
    private:
        unsigned int index;
        bool visited{};

    protected:
        double value{0};
        double activationValue{0};
        virtual void updateActivationValue();

    public:
        /**
         * Create a nev vertex object.
         * @param index vertex index
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
        [[nodiscard]] double getValue() const;

        [[nodiscard]] double getActivationValue() const;

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
         * Sets the indexArg of this vertex.
         * @param indexArg indexArg to set
         */
        void setIndex(unsigned int indexArg);

        /**
         * Has this vertex already been visited during traversal?
         * @return true or false
         */
        [[nodiscard]] bool isVisited() const;

        /**
         * Sets the visitedArg property.
         * @param visitedArg has this vertex been visitedArg
         */
        void setVisited(bool visitedArg);

        /**
         * Get a string with useful information about this object.
         * @param technical format for https://csacademy.com/app/graph_editor/
         * @return this object's information
         */
        virtual std::string toString(bool technical);

        /**
         * Vertex type getter virtual function. Do not use this as it is.
         * @return 0
         */
        virtual enums::VertexType getType() = 0;
    };
}

#endif //NEUROEVOLUTION_VERTEX_H
