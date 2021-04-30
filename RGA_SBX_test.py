
import numpy as np
import random
import matplotlib.pyplot as plt
from itertools import count
import math
plt.style.use('fivethirtyeight')
POPULATION_SIZE=1000
TOURNAMENT_SELECTION_SIZE=2
MUTATION_RATE=0.25
NUMBER_OF_ELITE_CHROMOSOMES=0
sizeList=[10,10]
xUpperList=[100,100]
xLowerList=[-100,-100]
a=0.5 # alpha for crossover
nc=0
sigma=[0.01,0.01]
class Chromosome:
    def __init__(self):
        self._genes=[]
        self._fitness=0
        for i in range(xUpperList.__len__()):
            a=np.random.uniform(low=xLowerList[i], high=xUpperList[i])
            self._genes.append(a)

    def get_fitness(self):
        # self._fitness=100*(self._genes[1]-self._genes[0]**2)**2+(1-self._genes[0])**2 
        OF=0
        for i,X in enumerate(self._genes):
            OF+=(X**2)-10*math.cos(2*math.pi)+10
        self._fitness=OF
        return self._fitness

    def get_genes(self):
        return self._genes
    
    def __str__(self):
        return self._genes.__str__()

class Population:
    def __init__(self,size):
        self._chromosomes=[]
        i=0
        while i< size:
            self._chromosomes.append(Chromosome())
            i+=1
    def get_chromosomes(self):return self._chromosomes

class GeneticAlgorithm:
    @staticmethod
    def evolve(pop):
        # return GeneticAlgorithm._crossover_population(pop)
        newPop=GeneticAlgorithm._mutate_population(GeneticAlgorithm._crossover_population(pop))
        return GeneticAlgorithm._survivorStage(pop,newPop)
    
    #creating the Crossover Population:
    @staticmethod
    def _crossover_population(pop):
        crossover_pop=Population(0) # HERE I DEFINED the empty chromosomes when put the zero
        for i in range(NUMBER_OF_ELITE_CHROMOSOMES):
            crossover_pop.get_chromosomes().append(pop.get_chromosomes()[i]) ## we will move this chromose to the next gen
        i=NUMBER_OF_ELITE_CHROMOSOMES
        ## here we will exclude the elite chromsome
        while i < POPULATION_SIZE//2:
            # crossover population will have population after the selection.Then we select the best 2 Chromosome
            chromosome1=GeneticAlgorithm._select_tournament_population(pop).get_chromosomes()[0] # here we call to the _tournament_population
            chromosome2=GeneticAlgorithm._select_tournament_population(pop).get_chromosomes()[0]
            offspring1,offspring2 =GeneticAlgorithm.sbx_crossover_chromosomes(chromosome1, chromosome2)
            crossover_pop.get_chromosomes().append(offspring1)
            crossover_pop.get_chromosomes().append(offspring2)
            i+=1
        return crossover_pop
    #creating the Mutation Population:
    @staticmethod
    def _mutate_population(pop):
        for i in range(NUMBER_OF_ELITE_CHROMOSOMES,POPULATION_SIZE):# EXCLUDE THE ELITE CHROMOSOME
            GeneticAlgorithm._Normaldistribution_mutate_chromosome(pop.get_chromosomes()[i])
        return pop
    
    @staticmethod
    def sbx_crossover_chromosomes(chromosome1,chromosome2): # this method does the random gen selection from each one of the parent 
        ## chromosomes
        crossover_chrom1 = Chromosome()
        crossover_chrom2 = Chromosome()
        # Loop over the Number-variable
        for i in range(chromosome1.get_genes().__len__()):
            ### have to select the two parent and assumption is the p1 < p2 always
            x1=chromosome1.get_genes()[i] 
            x2=chromosome2.get_genes()[i]
            if x2 > x1:
                p2=x2
                p1=x1
            else:
                p2=x1
                p1=x2
            # dt= a*(p2-p1) # define dt 
            u=random.random() # randomly picked the 
            # calculation of beta 
            if u<=0.5:
                beta=(2*u)**(1 / (nc+1))
            else:
                beta=(1/ (2* (1-u)) )**(1 / (nc+1))
            o1 = 0.5* (  (p1+p2)-beta*(p2-p1) )
            o2 = 0.5* (  (p1+p2)+beta*(p2-p1) )
            crossover_chrom1.get_genes()[i]=o1
            crossover_chrom2.get_genes()[i]=o2
            # crossover_chrom.get_genes()[i]=p1*(1-gamma)+gamma*p2
            # print(crossover_chrom.get_genes()[i])
        return crossover_chrom1,crossover_chrom2
            

    
    @staticmethod
    def _blx_alpha_crossover(chromosome1,chromosome2):
        crossover_chrom = Chromosome()
        # Loop over the Number-variable
        for i in range(chromosome1.get_genes().__len__()):
            ### have to select the two parent and assumption is the p1 < p2 always
            x1=chromosome1.get_genes()[i] 
            x2=chromosome2.get_genes()[i]
            if x2 > x1:
                p2=x2
                p1=x1
            else:
                p2=x1
                p1=x2
            # dt= a*(p2-p1) # define dt 
            u=random.random() # randomly picked the 
            gamma=(1+2*a)*u-a
            crossover_chrom.get_genes()[i]=p1*(1-gamma)+gamma*p2
            # print(crossover_chrom.get_genes()[i])
        return crossover_chrom
            

    @staticmethod
    def _Normaldistribution_mutate_chromosome(chromosome):
        for i in range(chromosome.get_genes().__len__()):
            chromosome.get_genes()[i]=chromosome.get_genes()[i]+np.random.normal(0, sigma[i])
            if chromosome.get_genes()[i] > xUpperList[i]:
                chromosome.get_genes()[i]=xUpperList[i]
            if chromosome.get_genes()[i] < xLowerList[i]:
                chromosome.get_genes()[i]=xLowerList[i]


    @staticmethod
    def _select_tournament_population(pop):
        tournament_pop=Population(0)
        i=0
        while i < TOURNAMENT_SELECTION_SIZE:
            tournament_pop.get_chromosomes().append(pop.get_chromosomes()[random.randrange(0,POPULATION_SIZE)])
            i+=1
        tournament_pop.get_chromosomes().sort(key=lambda x:x.get_fitness(), reverse=False)
        return tournament_pop
    @staticmethod
    
    def _survivorStage(oldpop,newpop):
        survivalPop=Population(0)
        oldpop.get_chromosomes().sort(key=lambda x:x.get_fitness(), reverse=False)
        newpop.get_chromosomes().sort(key=lambda x:x.get_fitness(), reverse=False)
        n=POPULATION_SIZE//2
        i=0
        while i < n:
            survivalPop.get_chromosomes().append(oldpop.get_chromosomes()[i])
            i+=1
        j=0
        while j < n:
            survivalPop.get_chromosomes().append(newpop.get_chromosomes()[j])
            j+=1
        return survivalPop



def _print_population(pop, gen_number):
    print('\n------------------------------------------------')
    print("Generation #", gen_number, "|Fittest Chromosome fitness :", pop.get_chromosomes()[0].get_fitness())
    print("---------------------------------------------------")
    i=0
    for x in pop.get_chromosomes():
        print("Chromosome #",i,":", x, "| fitness :", x.get_fitness())
        i+=1
population=Population(POPULATION_SIZE)
population.get_chromosomes().sort(key=lambda x:x.get_fitness(), reverse=False)
_print_population(population,0)
# newpopulation=GeneticAlgorithm().evolve(population)
# newpopulation.get_chromosomes().sort(key=lambda x:x.get_fitness(), reverse=False)
# _print_population(newpopulation,1)

generation_number=1
plt.figure()
while generation_number< 100:
    population=GeneticAlgorithm().evolve(population)
    population.get_chromosomes().sort(key=lambda x:x.get_fitness(), reverse=False)
    plt.scatter(generation_number,population.get_chromosomes()[0].get_fitness() )
    _print_population(population,generation_number)
    plt.pause(0.01)
    generation_number+=1
print('Minmization -- at the function :',population.get_chromosomes()[0])
plt.show()



