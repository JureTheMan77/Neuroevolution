//
// Created by jure on 24. 01. 22.
//

#ifndef NEUROEVOLUTION_ICROSSOVERABLE_H
#define NEUROEVOLUTION_ICROSSOVERABLE_H

namespace data_structures {
    class ICrossoverable {
    private:
        bool dominant{};
    public:

        /**
         * Create a new ICrossoverable object with dominance of your choice.
         * @param dominant influences crossover step
         */
        explicit ICrossoverable(bool dominant) : dominant(dominant) {}

        /**
         * Is this object dominant?
         * @return true or false
         */
        [[nodiscard]] bool isDominant() const;

        /**
         * Set the object as dominantArg or recessive.
         * @param dominantArg dominantArg or recessive
         */
        void setDominant(bool dominantArg);
    };
}

#endif //NEUROEVOLUTION_ICROSSOVERABLE_H
