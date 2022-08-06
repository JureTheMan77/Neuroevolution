//
// Created by jure on 24. 01. 22.
//

#ifndef NEUROEVOLUTION_ICROSSOVERABLE_H
#define NEUROEVOLUTION_ICROSSOVERABLE_H

namespace data_structures {
    class ICrossoverable {
    private:
        double mutationChance{};
        bool dominant{};
        //double chanceToGetDominated{};
        unsigned int maxChildren{};
    public:

        explicit ICrossoverable(double mutationChance, bool dominant,
                                unsigned int maxChildren) : mutationChance(mutationChance), dominant(dominant),
                                                            maxChildren(maxChildren) {}

        /**
         * Get the mutation chance.
         * @return
         */
        double getMutationChance() const;

        /**
         * Is this object dominant?
         * @return true or false
         */
        bool isDominant() const;

        /**
         * Set the object as dominant or recessive.
         * @param dominant dominant or recessive
         */
        void setDominant(bool dominant);

        /**
         * Max number of children this object can produce.
         * @return value
         */
        unsigned int getMaxChildren() const;
    };
}


#endif //NEUROEVOLUTION_ICROSSOVERABLE_H
