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
        double chanceToGetDominated{};
        unsigned int maxChildren{};
    public:

        explicit ICrossoverable(double mutationChance, bool dominant, double chanceToGetDominated,
                                unsigned int maxChildren) : mutationChance(mutationChance), dominant(dominant),
                                                            chanceToGetDominated(chanceToGetDominated),
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
         * What's the chance of this object to be dominated?
         * @return value
         */
        double getChanceToGetDominated() const;

        /**
         * Max number of children this object can produce.
         * @return value
         */
        unsigned int getMaxChildren() const;
    };
}


#endif //NEUROEVOLUTION_ICROSSOVERABLE_H
