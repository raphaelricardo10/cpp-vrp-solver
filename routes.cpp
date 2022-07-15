#include <functional>
#include <tuple>
#include <deque>
#include <unordered_set>
#include <iostream>
#include <string>
#include "genetic.cpp"

namespace ga
{
    class Interval
    {
    public:
        int startIndex;
        int endIndex;
        std::vector<int> *v;

        Interval(std::vector<int> &v, int index1, int index2)
        {
            this->v = &v;
            this->startIndex = std::min(index1, index2);
            this->endIndex = std::max(index1, index2);
        }

        Interval(std::vector<int> &v, RandomizerInt &randomizer)
        {
            randomizer.set_range(v);
            int index1 = randomizer.get_number();
            int index2 = randomizer.get_number(index1);

            *this = Interval(v, index1, index2);
        }

        Interval(std::vector<int> &v, int bpIndex, RandomizerInt &randomizer)
        {
            std::vector<int> breakpoints;
            breakpoints.reserve(v.size() - bpIndex + 1);
            breakpoints.push_back(0);
            breakpoints.insert(breakpoints.begin() + 1, v.begin() + bpIndex, v.end());

            randomizer.set_range(breakpoints);

            int start = randomizer.get_number();
            int end = start < breakpoints.size() - 1 ? start + 1 : start - 1;

            *this = Interval(v, breakpoints[start], breakpoints[end]);
        }

        Interval(Interval &interval, std::vector<int> &v, RandomizerInt &randomizer)
        {
            v.insert(v.begin(), interval.begin(), interval.end());
            *this = Interval(v, randomizer);
        }

        std::vector<int>::iterator begin()
        {
            return this->v->begin() + startIndex;
        }

        std::vector<int>::iterator end()
        {
            return this->v->begin() + endIndex + 1;
        }

        std::vector<int>::iterator at(int pos)
        {
            return this->v->begin() + startIndex + pos;
        }

        int size()
        {
            return this->endIndex - this->startIndex + 1;
        }

        void rotate_left(int n)
        {
            for (int i = 0; i < n; i++)
            {
                this->v->insert(this->end() + 1, *this->begin());
                this->v->erase(this->begin());
            }
        }
    };

    class Crossover
    {
    private:
        bool initialized;
        int numberOfTrials;

        template <class _ContainerType, class _ElementType>
        bool isInContainer(_ContainerType container, _ElementType elem)
        {
            return container.find(elem) != container.end();
        }

        bool isInRange(int base, int start, int end, int value)
        {
            return value >= base + start && value <= base + end;
        }

    public:
        int maxOfTrials;
        Individual *parent1;
        Individual *parent2;
        Individual offspring;

        Crossover()
        {
            this->maxOfTrials = 0;
            this->parent1 = 0;
            this->parent2 = 0;
            this->numberOfTrials = 0;
            this->initialized = false;
        }

        Crossover(Individual &parent1, Individual &parent2, int maxOfTrials)
        {
            this->parent1 = &parent1;
            this->parent2 = &parent2;
            this->maxOfTrials = maxOfTrials;
            this->initialized = true;
            this->numberOfTrials = 0;
        }

        bool is_acceptable()
        {
            if (!this->initialized)
            {
                return false;
            }

            if (this->offspring.fitness == 0)
            {
                return false;
            }

            if (this->numberOfTrials > this->maxOfTrials)
            {
                return false;
            }

            if (this->offspring.fitness > this->parent1->fitness)
            {
                return false;
            }

            if (this->offspring.fitness > this->parent2->fitness)
            {
                return false;
            }

            return true;
        }

        void make_offspring(int bpIndex, RandomizerInt &randomizer)
        {
            Interval p1Interval(this->parent1->chromossome.genes, bpIndex, randomizer);
            Interval p2Interval(this->parent2->chromossome.genes, bpIndex, randomizer);

            std::vector<int> p1Part;
            Interval crossoverInterval(p1Interval, p1Part, randomizer);

            std::deque<int> p2Part(p2Interval.begin(), p2Interval.end());

            std::unordered_set<int> crossoverMap(crossoverInterval.begin(), crossoverInterval.end());

            int rotationOffset = crossoverInterval.size() / 2;
            rotate_deq(p2Part, rotationOffset);

            std::deque<int> offspring(this->parent2->chromossome.genes.begin(), this->parent2->chromossome.genes.begin() + bpIndex);

            int insertionPoint = std::min((int)p2Part.size() - 1, crossoverInterval.startIndex);

            p2Part.insert(p2Part.begin() + insertionPoint, crossoverInterval.begin(), crossoverInterval.end());
            offspring.erase(offspring.begin() + p2Interval.startIndex, offspring.begin() + p2Interval.endIndex + 1);
            offspring.insert(offspring.begin() + p2Interval.startIndex, p2Part.begin(), p2Part.end());

            int i = -1;
            auto it = std::remove_if(offspring.begin(), offspring.end(), [&crossoverMap, &p2Interval, &crossoverInterval, &i, insertionPoint, this](int elem)
                                     {
                i++;

                if(!this->isInContainer(crossoverMap, elem)){
                    return false;
                }

                if(!this->isInRange(p2Interval.startIndex, insertionPoint, insertionPoint + crossoverInterval.size() - 1, i)){
                    return true;
                }

                return false; });

            offspring.erase(it, offspring.end());
            offspring.insert(offspring.end(), this->parent2->chromossome.genes.begin() + bpIndex, this->parent2->chromossome.genes.end());

            this->numberOfTrials++;
            this->offspring = Individual(std::vector<int>(offspring.begin(), offspring.end()));
        }
    };

    class RoutingGA : public GeneticBase
    {
    private:
        bool should_update_best(int fitness)
        {
            if (!this->population.best)
            {
                return true;
            }

            return fitness < this->population.best->fitness;
        }

        int get_distance(int location1, int location2)
        {
            int i = std::max(location1, location2);
            int j = std::min(location1, location2);

            return this->distances[i][j];
        }

        std::vector<std::vector<int>> vector_from_pointer(int *ptr, int width, int height)
        {
            std::vector<std::vector<int>> v(height, std::vector<int>(width));

            for(int i = 0; i < height; i++)
            {
                for(int j = 0; j < width; j++)
                {
                    v[i][j] = ptr[width * i + j];
                }
            }

            return v;
        }

        void map_routes(Individual &individual, std::function<void(int begin, int end)> func){
            for (int i = this->numberOfLocations - 1; i < individual.chromossome.genes.size(); i++)
            {
                int firstLocationOfRoute = i == this->numberOfLocations - 1 ? 0 : individual.chromossome.genes[i];
                int LastLocationOfRoute = i + 1 >= individual.chromossome.genes.size() ? this->numberOfLocations : individual.chromossome.genes[i + 1];

                func(firstLocationOfRoute, LastLocationOfRoute);
            }
        }

    public:
        int numberOfRoutes;
        int numberOfLocations;
        float optRate;
        std::vector<std::vector<int>> distances;

        RoutingGA(int maxGenerations, int populationSize, int numLocations, int numRoutes, int selectionK, float mutationRate, float optRate)
        {
            this->maxGenerations = maxGenerations;
            this->numberOfRoutes = numRoutes;
            this->selectionK = selectionK;
            this->mutationRate = mutationRate;
            this->numberOfLocations = numLocations;
            this->population = Population(populationSize, numLocations, numRoutes);
            this->generate_google_distances();
        }

        RoutingGA(int maxGenerations, int populationSize, int numLocations, int numRoutes, int selectionK, float mutationRate, float optRate,int *v_ptr)
        {
            this->maxGenerations = maxGenerations;
            this->numberOfRoutes = numRoutes;
            this->selectionK = selectionK;
            this->mutationRate = mutationRate;
            this->optRate = optRate;
            this->numberOfLocations = numLocations;
            this->population = Population(populationSize, numLocations, numRoutes);
            this->distances = this->vector_from_pointer(v_ptr, numLocations + 1, numLocations + 1);
        }

        void generate_distances(int individualSize)
        {
            RandomizerInt randomizer(9000, 10000);
            for (int i = 0; i <= individualSize; i++)
            {
                std::vector<int> v;
                for (int j = 0; j < i; j++)
                {
                    v.push_back(randomizer.get_number());
                }
                this->distances.push_back(v);
            }
        }

        void generate_google_distances()
        {
            this->distances = {
                {0, 548, 776, 696, 582, 274, 502, 194, 308, 194, 536, 502, 388, 354, 468,
                 776, 662},
                {548, 0, 684, 308, 194, 502, 730, 354, 696, 742, 1084, 594, 480, 674,
                 1016, 868, 1210},
                {776, 684, 0, 992, 878, 502, 274, 810, 468, 742, 400, 1278, 1164, 1130,
                 788, 1552, 754},
                {696, 308, 992, 0, 114, 650, 878, 502, 844, 890, 1232, 514, 628, 822,
                 1164, 560, 1358},
                {582, 194, 878, 114, 0, 536, 764, 388, 730, 776, 1118, 400, 514, 708,
                 1050, 674, 1244},
                {274, 502, 502, 650, 536, 0, 228, 308, 194, 240, 582, 776, 662, 628, 514,
                 1050, 708},
                {502, 730, 274, 878, 764, 228, 0, 536, 194, 468, 354, 1004, 890, 856, 514,
                 1278, 480},
                {194, 354, 810, 502, 388, 308, 536, 0, 342, 388, 730, 468, 354, 320, 662,
                 742, 856},
                {308, 696, 468, 844, 730, 194, 194, 342, 0, 274, 388, 810, 696, 662, 320,
                 1084, 514},
                {194, 742, 742, 890, 776, 240, 468, 388, 274, 0, 342, 536, 422, 388, 274,
                 810, 468},
                {536, 1084, 400, 1232, 1118, 582, 354, 730, 388, 342, 0, 878, 764, 730,
                 388, 1152, 354},
                {502, 594, 1278, 514, 400, 776, 1004, 468, 810, 536, 878, 0, 114, 308,
                 650, 274, 844},
                {388, 480, 1164, 628, 514, 662, 890, 354, 696, 422, 764, 114, 0, 194, 536,
                 388, 730},
                {354, 674, 1130, 822, 708, 628, 856, 320, 662, 388, 730, 308, 194, 0, 342,
                 422, 536},
                {468, 1016, 788, 1164, 1050, 514, 514, 662, 320, 274, 388, 650, 536, 342,
                 0, 764, 194},
                {776, 868, 1552, 560, 674, 1050, 1278, 742, 1084, 810, 1152, 274, 388,
                 422, 764, 0, 798},
                {662, 1210, 754, 1358, 1244, 708, 480, 856, 514, 468, 354, 844, 730, 536,
                 194, 798, 0},
            };
        }

        void calculate_fitness(Individual &individual)
        {
            int totalDistance = 0;

            this->map_routes(individual, [this, &totalDistance, individual] (int firstLocation, int lastLocation) {
                for (int i = firstLocation; i <= lastLocation; i++)
                {
                    int currLocation = i == lastLocation ? 0 : individual.chromossome.genes[i];
                    int prevLocation = i == firstLocation ? 0 : individual.chromossome.genes[i - 1];

                    totalDistance += this->get_distance(currLocation, prevLocation);
                }
            });

            individual.fitness = totalDistance;
        }

        void print_routes(Individual &individual){
            int routeNum = 1;
            this->map_routes(individual, [this, individual, &routeNum] (int firstLocation, int lastLocation) {
                std::string message("Route: " + std::to_string(routeNum) + ": 0 >> ");

                for(int i = firstLocation; i < lastLocation; i++){
                    if(i != firstLocation){
                        message += " >> ";
                    }

                    message += std::to_string(individual.chromossome.genes[i]);
                }

                routeNum++;
                std::cout << message << '\n';
            });
        }

        void two_opt()
        {
            Randomizer<std::uniform_real_distribution<float>, float> random_float(0, 1);

            this->population.map([this, &random_float] (Individual &individual){

                if(random_float.get_number() < this->optRate){
                    Interval randomRoute(individual.chromossome.genes, this->numberOfLocations,this->randomizer);
                    this->randomizer.set_range(randomRoute.startIndex, randomRoute.endIndex);

                    int best_j = 0;
                    int best_distance = 0;

                    int best_i = randomRoute.startIndex;

                    for(int i = randomRoute.startIndex + 1; i < randomRoute.endIndex - 1; i++){
                        int i_neigh = this->get_distance(i, i-1) + this->get_distance(i, i+1);

                        for(int j = i + 1; j < randomRoute.endIndex - 1; j++){
                            int j_neigh = this->get_distance(j, j-1) + this->get_distance(j, j+1);

                            int j_in_i = this->get_distance(j, i-1) + this->get_distance(j, i+1);
                            int i_in_j = this->get_distance(i, j-1) + this->get_distance(i, j+1);

                            int currTotalDistance = i_neigh + j_neigh;
                            int newTotalDistance = i_in_j + j_in_i;

                            if(best_j == 0 && best_i == randomRoute.startIndex){
                                best_j = j;
                                best_distance = std::min(currTotalDistance, newTotalDistance);
                            }

                            if(newTotalDistance < currTotalDistance && newTotalDistance < best_distance){
                                best_j = j;
                            }
                        }
                    }

                    int aux = individual.chromossome.genes[best_i];
                    individual.chromossome.genes[best_i] = individual.chromossome.genes[best_j];
                    individual.chromossome.genes[best_j] = aux;

                }
            });
        }

        void make_mutation()
        {
            Randomizer<std::uniform_real_distribution<float>, float> random_float(0, 1);

            this->population.map([this, &random_float] (Individual &individual){
                random_float.get_number();

                if(random_float.get_number() < this->mutationRate){
                    this->randomizer.set_range(0, this->numberOfLocations - 1);
                    int index1 = this->randomizer.get_number();
                    int index2 = this->randomizer.get_number(index1);

                    int aux = individual.chromossome.genes[index1];

                    individual.chromossome.genes[index1] = individual.chromossome.genes[index2];
                    individual.chromossome.genes[index2] = aux;

                    Interval routesInterval(individual.chromossome.genes, this->numberOfLocations, individual.chromossome.genes.size());
                    this->randomizer.set_range(1, this->numberOfLocations - 1);
                    std::unordered_set<int> bp_map(routesInterval.begin(), routesInterval.end() - 1);

                    int new_bp = this->randomizer.get_number(bp_map);
                    this->randomizer.set_range(routesInterval.startIndex, routesInterval.endIndex - 1);
                    int pos = this->randomizer.get_number();

                    individual.chromossome.genes[pos] = new_bp;
                }
            });
        }

        void run()
        {
            this->population.map([this](Individual &individual)
                                 {
                                     this->calculate_fitness(individual);

                                     if (this->should_update_best(individual.fitness))
                                     {
                                         this->population.best = &individual;
                                     } });

            for ( ; this->population.generation < this->maxGenerations; this->population.generation++)
            {

                int p1, p2;
                Crossover crossover1;
                Crossover crossover2;
                int maxOfTries = 15;
                int tries = 1;

                while ((!crossover1.is_acceptable() && !crossover2.is_acceptable()) || tries <= maxOfTries)
                {
                    std::tie(p1, p2) = this->make_selection();

                    crossover1 = Crossover(this->population.individuals[p1], this->population.individuals[p2], 5);
                    crossover2 = Crossover(this->population.individuals[p2], this->population.individuals[p1], 5);

                    for (int i = 0; i < 5; i++)
                    {
                        if (!crossover1.is_acceptable())
                        {
                            crossover1.make_offspring(this->numberOfLocations, this->randomizer);
                            this->calculate_fitness(crossover1.offspring);
                        }
                        if (!crossover2.is_acceptable())
                        {
                            crossover2.make_offspring(this->numberOfLocations, this->randomizer);
                            this->calculate_fitness(crossover2.offspring);
                        }
                    }

                    tries++;
                }

                this->population.individuals[p1] = crossover1.offspring;
                this->population.individuals[p2] = crossover2.offspring;

                if(this->should_update_best(this->population.individuals[p1].fitness)){
                    this->population.best = &this->population.individuals[p1];
                }

                if(this->should_update_best(this->population.individuals[p2].fitness)){
                    this->population.best = &this->population.individuals[p2];
                }

                this->two_opt();
                this->make_mutation();
                
                // std::cout << this->population.generation << '\t' << this->population.best->fitness << '\n';
            }

            // this->print_routes(*this->population.best);
        }
    };
    
    extern "C"
    void ga_interface(int popSize, int qtyLocations, int qtyRoutes, int maxGenerations, int selectionK, float mutationRate, float optRate, int *distances, int *result)
    {
        ga::RoutingGA ga(maxGenerations, popSize, qtyLocations, qtyRoutes, selectionK, mutationRate, optRate, distances);

        ga.run();

        // std::cout << ga.population.best->fitness << '\n';

        for(int i = 0; i < ga.population.best->chromossome.genes.size(); i++){
            result[i] = ga.population.best->chromossome.genes[i];
        }
    }
}
