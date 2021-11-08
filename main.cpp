#include "data_structures/Graph.h"
#include "logging//logging.h"
#include "evolution/Population.h"

int main() {
    logging::logs("Pozdravljen, svet!");

    auto inputLabels = std::vector<std::string>();
    inputLabels.insert(inputLabels.end(), {"sepal length", "sepal width", "petal length", "petal width"});
    auto outputLabels = std::vector<std::string>();
    outputLabels.insert(outputLabels.end(), {"Iris Setosa", "Iris Versicolour", "Iris Virginica"});

    auto graph = data_structures::Graph::createGraph(4,inputLabels,5,3,outputLabels);
    graph->addEdge(enums::VertexType::Input,0,enums::VertexType::Deep,1,2,1);
    graph->addEdge(enums::VertexType::Deep,1,enums::VertexType::Output,2,1,1);
    logging::logs(graph->toString());

    auto agent = evolution::Agent::create(graph, 0);
    logging::logs(agent->toString());

//    auto *pop = new evolution::Population(10, 4, inputLabels, 15, 3, outputLabels, 10, true);


    //logging::logs(pop->toString());

    //delete pop;

//    auto *pop = new evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/iris/iris.data");
//    pop->initialisePopulation(50,20,50,2,true);
//    logging::logs("Population initialised.");
//    pop->calculateFitness();
//    logging::logs("Fitness calculated.");
//    pop->sample(enums::SelectionType::StochasticUniversalSampling, 25);
//    logging::logs("Population sampled.");
//    //logging::logs(pop->toString());
//    delete pop;


    logging::logs("done");

    return 0;
}