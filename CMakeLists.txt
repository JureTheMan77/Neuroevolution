cmake_minimum_required(VERSION 3.15)
project(Neuroevolution)

set(CMAKE_CXX_STANDARD 17)

add_executable(Neuroevolution main.cpp data_structures/Vertex.cpp data_structures/Vertex.h data_structures/Edge.cpp data_structures/Edge.h data_structures/Graph.cpp data_structures/Graph.h data_structures/Edge.fwd.h data_structures/Vertex.fwd.h enums/VertexType.h logging/logging.h enums/EnumUtil.cpp enums/EnumUtil.h evolution/Population.cpp evolution/Population.h data_structures/InputVertex.cpp data_structures/InputVertex.h data_structures/DeepVertex.cpp data_structures/DeepVertex.h data_structures/OutputVertex.cpp data_structures/OutputVertex.h util/util.h data_structures/DataInstance.cpp data_structures/DataInstance.h data_structures/UniqueVertexQueue.cpp data_structures/UniqueVertexQueue.h evolution/Agent.cpp evolution/Agent.h enums/SelectionType.h data_structures/InputEdges.cpp data_structures/InputEdges.h data_structures/EdgeGroup.cpp data_structures/EdgeGroup.h data_structures/OutputEdges.cpp data_structures/OutputEdges.h data_structures/IDeepCloneable.h data_structures/ICrossoverable.cpp data_structures/ICrossoverable.h data_structures/MulticlassConfusionMatrix.cpp data_structures/MulticlassConfusionMatrix.h enums/FitnessMetric.h nlohmann/json.hpp evolution/Metrics.cpp evolution/Metrics.h)
#add_executable(Neuroevolution main.cpp data_structures/Vertex.cpp data_structures/Vertex.h data_structures/Edge.cpp data_structures/Edge.h data_structures/Graph.cpp data_structures/Graph.h data_structures/Edge.fwd.h data_structures/Vertex.fwd.h enums/VertexType.h logging/logging.h enums/EnumUtil.cpp enums/EnumUtil.h data_structures/InputVertex.cpp data_structures/InputVertex.h data_structures/DeepVertex.cpp data_structures/DeepVertex.h data_structures/OutputVertex.cpp data_structures/OutputVertex.h util/util.h data_structures/DataInstance.cpp data_structures/DataInstance.h data_structures/UniqueVertexQueue.cpp data_structures/UniqueVertexQueue.h evolution/Agent.cpp evolution/Agent.h enums/SelectionType.h data_structures/InputEdges.cpp data_structures/InputEdges.h data_structures/EdgeGroup.cpp data_structures/EdgeGroup.h data_structures/OutputEdges.cpp data_structures/OutputEdges.h)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
target_link_libraries(Neuroevolution PRIVATE Threads::Threads)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_CXX_FLAGS "-Wextra -Wall -march='native' -mtune='native'")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -ggdb3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")