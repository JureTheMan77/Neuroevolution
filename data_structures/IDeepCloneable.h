//
// Created by jure on 12. 11. 21.
//

#ifndef NEUROEVOLUTION_IDEEPCLONEABLE_H
#define NEUROEVOLUTION_IDEEPCLONEABLE_H

#include <type_traits>
#include <memory>

namespace data_structures {
    /**
     * Interface that adds deep cloning functionality.
     * @tparam T target class type
     */
    template<typename T>
    class IDeepCloneable {
    public:
        /**
         * Deep clones this abject.
         * @return new instance, wrapped in a shared pointer
         */
        virtual std::shared_ptr<T> deepClone() = 0;
    };
}


#endif //NEUROEVOLUTION_IDEEPCLONEABLE_H
