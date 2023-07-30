#include <filesystem>
#include <fstream>
#include "data_structures/Graph.h"
#include "logging//logging.h"
#include "evolution/Population.h"
#include "data_structures/MulticlassConfusionMatrix.h"
#include "util/util.h"

int main(int argc, char *argv[]) {

    // process command line arguments
    // source https://stackoverflow.com/a/6034177
    std::vector<std::string> cmdArgs(&argv[0], &argv[0 + argc]);

    // first argument: path to dataset
    std::string pathToDataset = cmdArgs.at(1);
    // second argument: population size
    unsigned int populationSize = std::stoi(cmdArgs.at(2));
    // third argument: maximum number of deep vertices
    unsigned int maxDeepVertexCount = std::stoi(cmdArgs.at(3));
    // fourth argument: maximum number od edges
    unsigned int maxEdgeCount = std::stoi(cmdArgs.at(4));
    // fifth argument: maximum number of times edges can be traversed
    unsigned int edgeTraverseLimit = std::stoi(cmdArgs.at(5));
    // sixth argument: should dormant vertices and edges be kept (preform minimization step after every crossover + mutation step)
    bool keepDormant = cmdArgs.at(6) == "true";
    // seventh argument: percentage of children to be mutated
    double mutationChance = std::stod(cmdArgs.at(7));
    // eighth argument: agents to keep after sampling
    unsigned int agentsToKeep = std::stoi(cmdArgs.at(8));
    // ninth argument: keep fittest agent
    bool keepFittestAgent = cmdArgs.at(9) == "true";
    // tenth argument: edge contribution
    double edgeContribution = std::stod(cmdArgs.at(10));
    // eleventh argument: vertex contribution
    double vertexContribution = std::stod(cmdArgs.at(11));
    // twelfth argument: number of iterations
    unsigned int numIterations = std::stoi(cmdArgs.at(12));
    // thirteenth argument: fitness metric, default is Accuracy
    std::string fitnessMetricTemp = cmdArgs.at(13);
    enums::FitnessMetric fitnessMetric = enums::FitnessMetric::Accuracy;
    if (fitnessMetricTemp == "MatthewsCorrelationCoefficient") {
        fitnessMetric = enums::FitnessMetric::MatthewsCorrelationCoefficient;
    }


    // iris
    //auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/iris/iris.data");
    // wine
    auto pop = evolution::Population(pathToDataset);
    //car
    //auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/car/car.data");
    //statlog
    //auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/statlog/shuttle.data");
    //transfusion
    //auto pop = evolution::Population("/home/jure/CLionProjects/Neuroevolution/datasets/transfusion/transfusion.data");

//    auto inputLabels = pop.getInputLabels();
//    auto outputLabels = pop.getOutputLabels();
//    auto graph = data_structures::Graph::createGraph(13, inputLabels, 3, 3, outputLabels);
//
//    graph->addEdge(enums::VertexType::Input, 0, enums::VertexType::Deep, 0, 0, 1, 1);
//    graph->addEdge(enums::VertexType::Deep, 0, enums::VertexType::Deep, 1, 0, 0.5, 1);
//    graph->addEdge(enums::VertexType::Deep, 1, enums::VertexType::Deep, 2, 0, 0.5, 1);
//    graph->addEdge(enums::VertexType::Deep, 2, enums::VertexType::Deep, 0, 0, 0.5, 1);
//    graph->addEdge(enums::VertexType::Deep, 0, enums::VertexType::Output, 0, 0, 0.5, 2);


    //graph->addEdge(enums::VertexType::Input, 1, enums::VertexType::Output, 0, 1, -0.00963, 3);
    //graph->addEdge(enums::VertexType::Input, 2, enums::VertexType::Output, 0, 2, -0.712, 1);
    //graph->addEdge(enums::VertexType::Input, 6, enums::VertexType::Output, 0, 3, 1, 2);
    //graph->addEdge(enums::VertexType::Input, 8, enums::VertexType::Output, 0, 4, -0.0214, 2);
    //graph->addEdge(enums::VertexType::Input, 8, enums::VertexType::Output, 1, 5, 0.289, 2);
    //graph->addEdge(enums::VertexType::Input, 9, enums::VertexType::Output, 1, 6, -0.719, 1);
    //graph->addEdge(enums::VertexType::Deep, 0, enums::VertexType::Output, 1, 7, 1, 5);
    //graph->addEdge(enums::VertexType::Deep, 1, enums::VertexType::Deep, 1, 8, -0.0459, 5);
    //graph->addEdge(enums::VertexType::Deep, 1, enums::VertexType::Deep, 2, 9, 0.193, 3);
    //auto agent = evolution::Agent::create(graph);
    //pop.addAgent(agent);
    //auto minimized = pop.minimizeAgent(agent);
    //pop.calculateFitness(enums::FitnessMetric::Accuracy,-0.001,-0.001);

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

    //auto graph = data_structures::Graph::createGraph(4,inputLabels,2,3,outputLabels);
    //graph->addEdge(enums::VertexType::Input,1,enums::VertexType::Output,2,0,-1.01,2);
    //graph->addEdge(enums::VertexType::Input,1,enums::VertexType::Output,1,1,-0.49755895629207125,2);
    //graph->addEdge(enums::VertexType::Input,1,enums::VertexType::Deep,0,2,0.13420725976756206,2);
    //graph->addEdge(enums::VertexType::Input,2,enums::VertexType::Output,0,3,-0.6040223462636878,2);
    //graph->addEdge(enums::VertexType::Input,2,enums::VertexType::Output,2,4,0.23177392817540032,2);
    //graph->addEdge(enums::VertexType::Deep,0,enums::VertexType::Output,1,5,-1.0,2);
    //graph->addEdge(enums::VertexType::Deep,0,enums::VertexType::Deep,1,6,1,2);
    //graph->addEdge(enums::VertexType::Deep,1,enums::VertexType::Deep,0,7,1,2);
    //auto agent = evolution::Agent::create(graph);
    //std::shared_ptr<data_structures::DataInstance> di = data_structures::DataInstance::createDataInstance({5.1,3.5,1.4,0.2,0});
    //agent->getGraph()->traverse(di);


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


    std::string csvHeader = "Generation;Worst;0-10%;10-20%;20-30%;30-40%;40-50%;50-60%;60-70%;70-80%;80-90%;90-100%;Best;Average";
    std::ofstream fitnessFile("fitness.csv");
    fitnessFile << csvHeader << std::endl;
    std::ofstream accuracyFile("accuracy.csv");
    accuracyFile << csvHeader << std::endl;
    std::ofstream mccFile("mcc.csv");
    mccFile << csvHeader << std::endl;

    logging::logs("Start");
    pop.initialisePopulation(populationSize, maxDeepVertexCount, maxEdgeCount, edgeTraverseLimit, keepDormant,
                             mutationChance);
    logging::logs("Population initialised.");
    //logging::logs(pop.toString());

    for (int i = 0; i < numIterations; i++) {
        // std::cout << "\033[2J" << "\033[1;1H" << std::endl;
        logging::logs("Generation " + std::to_string(i));
        evolution::Metrics m = pop.calculateFitness(fitnessMetric, vertexContribution,
                                                    edgeContribution);


        logging::logs("AVERAGES - fitness: " + std::to_string(m.getBottomFitnessPercentile(100)) +
                      ", accuracy: " + std::to_string(m.getBottomAccuracyPercentile(100)) +
                      ", Matthews correlation coefficient: " + std::to_string(m.getBottomMccPercentile(100)));
        auto fittestAgent = pop.getFittestAgent();
        logging::logs("FITTEST AGENT - fitness: " + std::to_string(fittestAgent->getFitness()) +
                      ", accuracy: " + std::to_string(fittestAgent->getAccuracy()) +
                      ", Matthews correlation coefficient: " +
                      std::to_string(fittestAgent->getMatthewsCorrelationCoefficient()));

        // write pop info to files
        fitnessFile << util::fitnessToCsv(i, m) << std::endl;
        accuracyFile << util::accuracyToCsv(i, m) << std::endl;
        mccFile << util::mccToCsv(i, m) << std::endl;
        // write entire population fitness into a file
        //std::ofstream fitnessGenFile("fitnessGen.csv");
        //int ccc = 0;
        //for(const auto fitness : m.getFitnessList()){
        //    fitnessGenFile << std::to_string(ccc) + ";" + util::doubleToCsv(fitness) << std::endl;
        //    ccc+=1;
        //}

        //logging::logs(pop.getFittestAgent()->toString());
        if (i == numIterations - 1) {
            fitnessFile.close();
            accuracyFile.close();
            mccFile.close();
        } else {
            pop.sample(enums::SelectionType::StochasticUniversalSampling, agentsToKeep, keepFittestAgent,
                       fitnessMetric);
            //logging::logs("Population sampled.");

            pop.crossoverAndMutate();
            //logging::logs("Population crossovered.");
        }
    }

    //fittestAgent->minimize();
    std::shared_ptr<evolution::Agent> bestAgent = pop.getFittestAgent();
    logging::logs(bestAgent->toString(true));
    auto bestMcm = data_structures::MulticlassConfusionMatrix(bestAgent, pop.getTestingValues(),
                                                              pop.getNumberOfOutputs());
    logging::logs(bestMcm.toString(pop.getOutputLabels()));
    bestAgent->getGraph()->reset();
    std::ofstream jsonFullFile("/home/jure/CLionProjects/Neuroevolution/visualization/graphFull.json");
    jsonFullFile << bestAgent->getGraph()->toForceGraphJson();
    jsonFullFile.close();

    auto minimizedBestAgent = pop.minimizeAgent(bestAgent);
    logging::logs(minimizedBestAgent->toString(true));
    auto minimizedMcm = data_structures::MulticlassConfusionMatrix(minimizedBestAgent, pop.getTestingValues(),
                                                                   pop.getNumberOfOutputs());
    logging::logs(minimizedMcm.toString(pop.getOutputLabels()));

    //logging::logs(pop->toString());
    //delete pop;


    auto json = minimizedBestAgent->getGraph()->toForceGraphJson();
    std::ofstream jsonFile("/home/jure/CLionProjects/Neuroevolution/visualization/graph.json");
    jsonFile << json;
    jsonFile.close();

    fitnessFile.close();
    accuracyFile.close();
    mccFile.close();

    logging::logs("done");


    return 0;
}