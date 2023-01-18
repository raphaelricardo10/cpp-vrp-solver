#include <iostream>
#include "routes.cpp"

void print_cromossome(int chromossome)
{
    std::cout << chromossome << ' ';
}

void print_individual(ga::Individual *individual)
{
    individual->chromossome.map(print_cromossome);
    std::cout << '\n';
}

int main()
{
    int popSize = 200;
    int qtyLocations = 15;
    int qtyRoutes = 1;
    int maxGenerations = 2000;
    int selectionK = 3;
    float mutationRate = 0.05;
    float optRate = 0.02;

    ga::RoutingGA ga(maxGenerations, popSize, qtyLocations, qtyRoutes, selectionK, mutationRate, optRate);

    ga.run();

    std::cout << ga.population.best->fitness << '\n';

    return 0;
}