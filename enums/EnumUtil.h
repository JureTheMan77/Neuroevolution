//
// Created by jure on 27. 10. 20.
//

#ifndef NEUROEVOLUTION_ENUMUTIL_H
#define NEUROEVOLUTION_ENUMUTIL_H

#include <string>
#include "VertexType.h"

namespace enums {
    class EnumUtil {
    public:
        /**
         * Returns the name of a VertexType enum.
         * @param vertexType enum type
         * @return enum name
         */
        static std::string VertexTypeToString(enums::VertexType vertexType);
    };
}

#endif //NEUROEVOLUTION_ENUMUTIL_H
