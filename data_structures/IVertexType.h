//
// Created by jure on 17. 11. 21.
//

#ifndef NEUROEVOLUTION_IVERTEXTYPE_H
#define NEUROEVOLUTION_IVERTEXTYPE_H

#include "../enums/VertexType.h"

namespace data_structures {
    class IVertexType {
    public:
        virtual enums::VertexType getType() = 0;
    };
}


#endif //NEUROEVOLUTION_IVERTEXTYPE_H
