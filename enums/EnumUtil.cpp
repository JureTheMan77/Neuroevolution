//
// Created by jure on 27. 10. 20.
//

#include "EnumUtil.h"

std::string enums::EnumUtil::VertexTypeToString(enums::VertexType vertexType) {
    switch (vertexType) {
        case enums::VertexType::Input:
            return "Input";
        case enums::VertexType::Output:
            return "Output";
        case enums::VertexType::Deep:
            return "Deep";
        default:
            return "";
    }
}
