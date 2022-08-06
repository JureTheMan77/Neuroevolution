#include <filesystem>
#include <fstream>
#include "data_structures/Graph.h"
#include "logging//logging.h"
#include "evolution/Population.h"
#include "data_structures/MulticlassConfusionMatrix.h"

int main() {
    // iris
    auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/iris/iris.data");
    // wine
    //auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/wine/wine.data");
    //car
    //auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/car/car.data");
    //statlog
    //auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/statlog/shuttle.data");
    //transfusion
    //auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/transfusion/transfusion.data");

    //auto inputLabels = pop.getInputLabels();
    //auto outputLabels = pop.getOutputLabels();
    //auto graph = data_structures::Graph::createGraph(4, inputLabels, 16, 3, outputLabels, 1);
    //graph->addEdge(enums::VertexType::Input, 0, enums::VertexType::Output, 2, 0, 0.499, 1, 1);
    //graph->addEdge(enums::VertexType::Input, 0, enums::VertexType::Output, 0, 0, 0.0133, 1, 1);
    //graph->addEdge(enums::VertexType::Input, 1, enums::VertexType::Output, 0, 0, 0.414, 1, 1);
    //graph->addEdge(enums::VertexType::Input, 1, enums::VertexType::Deep, 6, 0, 0.0211, 1, 1);
    //graph->addEdge(enums::VertexType::Input, 1, enums::VertexType::Output, 2, 0, 0.72, 2, 1);
    //graph->addEdge(enums::VertexType::Input, 1, enums::VertexType::Deep, 10, 0, 0.147, 1, 1);
    //graph->addEdge(enums::VertexType::Input, 2, enums::VertexType::Output, 1, 0, 0.17, 2, 1);
    //graph->addEdge(enums::VertexType::Input, 2, enums::VertexType::Output, 2, 0, 0.0341, 2, 1);
    //graph->addEdge(enums::VertexType::Input, 2, enums::VertexType::Deep, 15, 0, 0.0505, 1, 1);
    //graph->addEdge(enums::VertexType::Input, 2, enums::VertexType::Output, 0, 0, 0.00873, 2, 1);
    //graph->addEdge(enums::VertexType::Input, 3, enums::VertexType::Output, 2, 0, 0.147, 2, 1);
    //graph->addEdge(enums::VertexType::Input, 3, enums::VertexType::Output, 0, 0, 0.134, 1, 1);
    //graph->addEdge(enums::VertexType::Input, 3, enums::VertexType::Output, 1, 0, 0.698, 1, 1);
    //graph->addEdge(enums::VertexType::Deep, 6, enums::VertexType::Output, 1, 0, 0.415, 2, 1);
    //graph->addEdge(enums::VertexType::Deep, 6, enums::VertexType::Deep, 6, 0, 0.0291, 1, 1);
    //graph->addEdge(enums::VertexType::Deep, 10, enums::VertexType::Output, 0, 0, 0.355, 1, 1);
    //graph->addEdge(enums::VertexType::Deep, 15, enums::VertexType::Output, 1, 0, 1, 1, 1);


    //auto agent = evolution::Agent::create(graph, 0);
    //auto minimizedAgent = pop.minimizeAgent(agent);
    //logging::logs(agent->toString(true));
    //auto mcmtest = data_structures::MulticlassConfusionMatrix(agent, pop.getTestingValues(), pop.getNumberOfOutputs());
    //logging::logs(mcmtest.toString(pop.getOutputLabels()));
//
    //mcmtest = data_structures::MulticlassConfusionMatrix(minimizedAgent, pop.getTestingValues(), pop.getNumberOfOutputs());
    //logging::logs(mcmtest.toString(pop.getOutputLabels()));
    //logging::logs(minimizedAgent->toString(true));

int a = 0;
    //graph->addEdge(enums::VertexType::Input, 1, enums::VertexType::Output, 1, 0.45, 1);
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


    std::ofstream fitnessFile("fitness.csv");


    logging::logs("Start");
    pop.initialisePopulation(1000, 50, 100, 4, true, 0.1);
    logging::logs("Population initialised.");
    //logging::logs(pop.toString());

    for (int i = 0; i < 1000; i++) {
        // std::cout << "\033[2J" << "\033[1;1H" << std::endl;
        logging::logs("Generation " + std::to_string(i));
        pop.calculateFitness(enums::FitnessMetric::Accuracy, 0.000, -0.001);
        //logging::logs("Fitness calculated.");
        //logging::logs(pop.toString());

        logging::logs("Average fitness: " + std::to_string(pop.getAverageFitness()));
        logging::logs("Fittest agent: " + std::to_string(pop.getFittestAgent()->getFitness()));

        // write pop info to files
        fitnessFile << pop.fitnessToCSVString(';', i) << std::endl;

        //logging::logs(pop.getFittestAgent()->toString());
        if (i == 999) {
            fitnessFile.close();
        } else {
            pop.sample(enums::SelectionType::StochasticUniversalSampling, 750);
            //logging::logs("Population sampled.");

            pop.crossoverAndMutate();
            //logging::logs("Population crossovered.");
        }
    }

    //fittestAgent->minimize();
    std::shared_ptr<evolution::Agent> bestAgent = pop.getPopulation().at(0);
    auto bestMcm = data_structures::MulticlassConfusionMatrix(bestAgent,pop.getTestingValues(), pop.getNumberOfOutputs());
    for (const auto &agent: pop.getPopulation()) {
        auto mcm = data_structures::MulticlassConfusionMatrix(agent, pop.getTestingValues(), pop.getNumberOfOutputs());
        if (bestMcm.getAccuracy() < mcm.getAccuracy()) {
            bestAgent = agent;
            bestMcm = mcm;
        }
        //if(bestAgent == nullptr || bestMcm.getMatthewsCorrelationCoefficient() < mcm.getMatthewsCorrelationCoefficient()){
        //    bestAgent = agent;
        //    bestMcm = mcm;
        //}
    }
    logging::logs(bestAgent->toString(true));
    logging::logs(bestMcm.toString(pop.getOutputLabels()));
    bestAgent->getGraph()->reset();
    auto minimizedBestAgent = pop.minimizeAgent(bestAgent);
    logging::logs(minimizedBestAgent->toString(true));
    auto mcm = data_structures::MulticlassConfusionMatrix(minimizedBestAgent, pop.getTestingValues(),
                                                          pop.getNumberOfOutputs());
    logging::logs(mcm.toString(pop.getOutputLabels()));

    //logging::logs(pop->toString());
    //delete pop;




    logging::logs("done");


    return 0;
}