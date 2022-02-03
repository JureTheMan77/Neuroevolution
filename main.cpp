#include "data_structures/Graph.h"
#include "logging//logging.h"
#include "evolution/Population.h"
#include "data_structures/MulticlassConfusionMatrix.h"

int main() {
    logging::logs("Pozdravljen, svet!");

    auto inputLabels = std::vector<std::string>();
    inputLabels.insert(inputLabels.end(),
                       {"sepal length in cm", "sepal width in cm", "petal length in cm", "petal width in cm"});
    auto outputLabels = std::vector<std::string>();
    outputLabels.insert(outputLabels.end(), {"Iris Setosa", "Iris Versicolour", "Iris Virginica"});

    auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/iris/iris.data");

//    auto graph = data_structures::Graph::createGraph(4,inputLabels,1,3,outputLabels);
    //graph->addEdge(enums::VertexType::Input,0,enums::VertexType::Deep,0,2,1);
    //graph->addEdge(enums::VertexType::Deep,0,enums::VertexType::Output,0,0.75,1);
    //graph->addEdge(enums::VertexType::Input,1,enums::VertexType::Output,1,0.45,1);
    // recursive
    //graph->addEdge(enums::VertexType::Input,2,enums::VertexType::Deep,1,1,1);
    //graph->addEdge(enums::VertexType::Deep,1,enums::VertexType::Deep,1,2,1);
    //graph->addEdge(enums::VertexType::Deep,1,enums::VertexType::Output,2,1,1);

//    graph->addEdge(enums::VertexType::Input,0,enums::VertexType::Deep,0,1,1);
//    graph->addEdge(enums::VertexType::Deep,0,enums::VertexType::Deep,0,2,1);
//    graph->addEdge(enums::VertexType::Deep,0,enums::VertexType::Output,0,1,1);
    //logging::logs(graph->toString());
//    auto agent = evolution::Agent::create(graph, 0);
//
//    std::shared_ptr<data_structures::DataInstance> di = data_structures::DataInstance::createDataInstance({5.1,3.5,1.4,0.2,0});
//    agent->getGraph()->traverse(di);
//    logging::logs(agent->toString());
//
//    logging::logs(agent->toString());
//    auto clonedAgent = agent->deepClone();
//
//    logging::logs(agent->toString());
//    logging::logs(clonedAgent->toString());

    //auto graph = data_structures::Graph::createGraph(4,inputLabels,0,3,outputLabels);
    //graph->addEdge(enums::VertexType::Input,1,enums::VertexType::Output,0,0.467950335030075,1);
    //graph->addEdge(enums::VertexType::Input,3,enums::VertexType::Output,1,1.31808422766154,1);

    //pop.addAgent(evolution::Agent::create(graph));
    //logging::logs(pop.toString());
    //pop.calculateFitness(-1, -0.5);
    //logging::logs(pop.getFittestAgent()->toString());
    logging::logs("Start");
    pop.initialisePopulation(1000, 20, 40, 2, true,0.2);
    logging::logs("Population initialised.");
    //logging::logs(pop.toString());

    std::shared_ptr<evolution::Agent> fittestAgent;
    for (int i = 0; i < 2000; i++) {
        // std::cout << "\033[2J" << "\033[1;1H" << std::endl;
        logging::logs("Generation " + std::to_string(i));
        pop.calculateFitness(enums::FitnessMetric::MatthewsCorrelationCoefficient,-0.002, -0.002);
        //logging::logs("Fitness calculated.");
        //logging::logs(pop.toString());

        logging::logs("Average fitness: " + std::to_string(pop.getAverageFitness()));
        logging::logs("Fittest agent: " + std::to_string(pop.getFittestAgent()->getFitness()));
        //logging::logs(pop.getFittestAgent()->toString());
        if (i == 1999) {
            fittestAgent = pop.getFittestAgent();
        } else {
            pop.sample(enums::SelectionType::StochasticUniversalSampling, 800);
            //logging::logs("Population sampled.");

            pop.crossover();
            //logging::logs("Population crossovered.");
        }
    }

    logging::logs(fittestAgent->toString());

    auto mcm = data_structures::MulticlassConfusionMatrix(fittestAgent, pop.getTestingValues(),
                                                          pop.getNumberOfOutputs());
    logging::logs(mcm.toString(pop.getOutputLabels()));

    //logging::logs(pop->toString());
    //delete pop;




    logging::logs("done");

    return 0;
}