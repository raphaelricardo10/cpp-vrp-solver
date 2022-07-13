#include <vector>
#include <algorithm>
#include <tuple>
#include <random>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <deque>

namespace ga
{
    template <class _DistributionType, class _DataType>
    class Randomizer
    {
    public:
        std::random_device rd;
        std::default_random_engine generator;
        _DistributionType distribution;

        Randomizer()
        {
            this->generator = std::default_random_engine(this->rd());
        }

        Randomizer(_DataType minValue, _DataType maxValue)
        {
            this->generator = std::default_random_engine(this->rd());
            this->distribution = _DistributionType(minValue, maxValue);
        }

        _DataType get_number()
        {
            return this->distribution(this->generator);
        }

        _DataType get_number(_DataType uniqueValue)
        {
            _DataType number;
            do
            {
                number = this->get_number();
            } while (number == uniqueValue);

            return number;
        }

        _DataType get_number(std::unordered_set<_DataType> &set)
        {
            _DataType number;
            do
            {
                number = this->get_number();
            } while (set.find(number) != set.end());

            return number;
        }

        template <class _Container>
        void set_range(_Container &container)
        {
            this->distribution = _DistributionType(0, container.size() - 1);
        }

        void set_range(_DataType minValue, _DataType maxValue)
        {
            this->distribution = _DistributionType(minValue, maxValue);
        }
    };

    typedef std::uniform_int_distribution<int> UniformIntDistribution;
    typedef Randomizer<UniformIntDistribution, int> RandomizerInt;

    void rotate_deq(std::deque<int> &deq, int qty)
    {
        for (int i = 0; i < qty; i++)
        {
            int front = deq.front();
            deq.pop_front();
            deq.push_back(front);
        }
    }

    class BreakpointSet
    {
    private:
        std::unordered_set<int> unordered_values;

    public:
        std::set<int> values;

        BreakpointSet(std::vector<int> &v, int qty, RandomizerInt &randomizer)
        {
            int numParts = v.size() / qty;
            for (int i = 0; i < qty; i++)
            {
                randomizer.set_range(numParts * i, std::min(numParts * (i + 1), (int) v.size() - 1));
                int breakpoint = randomizer.get_number(this->unordered_values);

                this->unordered_values.insert(breakpoint);
                this->values.insert(breakpoint);
            }
        }
    };

    class Permutator
    {
    private:
        std::default_random_engine rng;
        std::vector<int> generate_vector(int n)
        {
            std::vector<int> v(n);

            for (int i = 0; i < n; i++)
                v[i] = i + 1;

            return v;
        }

    public:
        int n;
        std::vector<int> vector;

        Permutator() {}

        Permutator(int n)
        {
            this->vector = this->generate_vector(n);
        }

        void shuffle()
        {
            std::shuffle(this->vector.begin(), this->vector.end(), this->rng);
        }
    };

    class Chromossome
    {
    public:
        std::vector<int> genes;

        void map(std::function<void(int)> func)
        {
            for (auto i = this->genes.begin(); i != this->genes.end(); ++i)
            {
                func(*i);
            }
        }
    };

    class Individual
    {
    private:
        void generate_breakpoints(int qty, RandomizerInt &randomizer)
        {
            BreakpointSet breakpoints(this->chromossome.genes, qty - 1, randomizer);

            this->chromossome.genes.reserve(breakpoints.values.size());
            this->chromossome.genes.insert(this->chromossome.genes.end(), breakpoints.values.begin(), breakpoints.values.end());
        }

    public:
        int fitness;
        Chromossome chromossome;

        Individual(std::vector<int> genes, int qtyBreaks, RandomizerInt &randomizer)
        {
            this->fitness = 0;
            this->chromossome.genes = genes;
            this->generate_breakpoints(qtyBreaks, randomizer);
        }

        Individual(std::vector<int> genes){
            this->fitness = 0;
            this->chromossome.genes = genes;
        }

        Individual()
        {
            this->fitness = 0;
        }
    };

    class Population
    {
    private:
        Permutator permutator;

    public:
        int generation;
        int size;
        Individual *best;
        std::vector<Individual> individuals;

        Population() {}

        Population(int size, int n, int qtyBreaks)
        {
            this->best = 0;
            this->generation = 0;
            this->size = size;
            this->permutator = Permutator(n);
            this->generate_individuals(qtyBreaks);
        }

        void map(std::function<void(Individual &)> func)
        {
            for (auto i = this->individuals.begin(); i != this->individuals.end(); ++i)
            {
                func(*i);
            }
        }

        void generate_individuals(int qtyBreaks)
        {
            this->permutator.shuffle();
            RandomizerInt randomizer(2, this->permutator.vector.size() - 1);

            for (int i = 0; i < this->size; i++)
            {
                this->permutator.shuffle();
                Individual individual(this->permutator.vector, qtyBreaks, randomizer);
                this->individuals.push_back(individual);
            }
        }
    };

    class GeneticBase
    {
    private:
        // Tournament selection
        int select_parent(std::unordered_set<int> &parents)
        {
            int winner = this->randomizer.get_number(parents);
            parents.insert(winner);

            for (int i = 1; i < this->selectionK; i++)
            {
                int chosen = this->randomizer.get_number(parents);
                parents.insert(chosen);

                if (this->population.individuals[i].fitness < this->population.individuals[winner].fitness)
                {
                    winner = chosen;
                }

                i++;
            }

            return winner;
        };

        int select_parent(std::unordered_set<int> &parents, int parent)
        {
            parents.insert(parent);
            return this->select_parent(parents);
        }

    public:
        int maxGenerations;
        int selectionK;
        float mutationRate;
        Population population;
        RandomizerInt randomizer;

        virtual void calculate_fitness(Individual &individual) = 0;

        auto make_selection()
        {
            this->randomizer.set_range(this->population.individuals);
            std::unordered_set<int> parents;
            int p1 = this->select_parent(parents);
            int p2 = this->select_parent(parents, p1);

            return std::make_tuple(p1, p2);
        }
    };
}
