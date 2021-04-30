
import numpy as np
import random
import matplotlib.pyplot as plt
from itertools import count
plt.style.use('fivethirtyeight')
POPULATION_SIZE=3000
TARGET_CHROMOSOME=[1,1,0,1,0,0,1,1,1,0,1,1,0,1,0,0,1,1,1,0]
TOURNAMENT_SELECTION_SIZE=2
MUTATION_RATE=0.25
NUMBER_OF_ELITE_CHROMOSOMES=2
sizeList=[10,10]
xUpperList=[5,5]
xLowerList=[-5,-5]
class Chromosome:
    def __init__(self):
        self._genes=[]
        self._chunks=0
        self._dcodedValue=0
        self._xRealValues=0
        self._fitness=0
        i=0
        while i < TARGET_CHROMOSOME.__len__():
            if random.random() >=0.5:
                self._genes.append(1)
            else:
                self._genes.append(0)
            i+=1
    def get_genes(self):
        return self._genes
    
    def get_fitness(self):
        self._fitness=0
        self._chunks=Chromosome.chunkIt(self._genes)
        self._dcodedValue=Chromosome.DecodedValue(self._chunks)
        self._xRealValues=Chromosome.real_values(self._dcodedValue,sizeList,xUpperList,xLowerList)
        self._fitness=100*(self._xRealValues[1]-self._xRealValues[0]**2)**2+(1-self._xRealValues[0])**2 
        return self._fitness
    
    @staticmethod
    def chunkIt(dataList):
        ## Size is in linst
        it = iter(dataList)
        return [[next(it) for _ in range(size)] for size in sizeList]
    @staticmethod
    def DecodedValue(dataList):
        DV=[]
        for i, eachVar_gens in enumerate(dataList):
            temp=0
            for j,value in enumerate(eachVar_gens[::-1]):# resverd for the fourmula ---b1*2^0+b2*2^1+....bn*2^(n-1)+
                temp=temp+(value*2**j)
            DV.append(temp)
        return DV
    
    @staticmethod
    def real_values(Dv,sizeList,xUpperList,xLowerList):
        real_value=[]
        for i, value in enumerate(Dv):
            pre=(xUpperList[i]-xLowerList[i])/(2**sizeList[i]-1)
            xi = xLowerList[i]+pre*value
            real_value.append(xi)
        return real_value

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
        newPop=GeneticAlgorithm._mutate_population(GeneticAlgorithm._crossover_population(pop))
        return GeneticAlgorithm._survivorStage(pop,newPop)
        
    
    @staticmethod
    def _crossover_population(pop):
        crossover_pop=Population(0) # HERE I DEFINED the empty chromosomes when put the zero
        for i in range(NUMBER_OF_ELITE_CHROMOSOMES):
            crossover_pop.get_chromosomes().append(pop.get_chromosomes()[i]) ## we will move the chromose to the next gen
        i=NUMBER_OF_ELITE_CHROMOSOMES
        ## here we will exclude the elite chromsome
        while i < POPULATION_SIZE:
            chromosome1=GeneticAlgorithm._select_tournament_population(pop).get_chromosomes()[0]
            chromosome2=GeneticAlgorithm._select_tournament_population(pop).get_chromosomes()[0]
            crossover_pop.get_chromosomes().append(GeneticAlgorithm._crossover_chromosomes(chromosome1, chromosome2))
            i+=1
        return crossover_pop


    @staticmethod
    def _mutate_population(pop):
        for i in range(NUMBER_OF_ELITE_CHROMOSOMES,POPULATION_SIZE):# EXCLUDE THE ELITE CHROMOSOME
            GeneticAlgorithm._mutate_chromosome(pop.get_chromosomes()[i])
        return pop
    @staticmethod
    def _crossover_chromosomes(chromosome1,chromosome2): # this method does the random gen selection from each one of the parent 
        ## chromosomes
        crossover_chrom = Chromosome()
        for i in range(TARGET_CHROMOSOME.__len__()):
            if random.random() < 0.5:
                crossover_chrom.get_genes()[i]=chromosome1.get_genes()[i]
            else:
                crossover_chrom.get_genes()[i]=chromosome2.get_genes()[i]
        return crossover_chrom
    
    @staticmethod
    def _mutate_chromosome(chromosome):
        for i in range(TARGET_CHROMOSOME.__len__()):
            if random.random()< MUTATION_RATE:
                if random.random() < 0.5:
                    chromosome.get_genes()[i]=1
                else:
                    chromosome.get_genes()[i]=0    
        
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
    print("Target Chromosome:", TARGET_CHROMOSOME)
    print("---------------------------------------------------")
    i=0
    for x in pop.get_chromosomes():
        print("Chromosome #",i,":", x, "| fitness :", x.get_fitness())

population=Population(POPULATION_SIZE)#--gen zero
# print(population.get_chromosomes()[0].get_fitness())
# print(population.get_chromosomes()[0]._chunks)
# print(population.get_chromosomes()[0]._dcodedValue)
# print(population.get_chromosomes()[0]._xRealValues)

population.get_chromosomes().sort(key=lambda x:x.get_fitness(), reverse=False)

_print_population(population,0)
generation_number=1
plt.figure()
while generation_number< 50:
    population=GeneticAlgorithm().evolve(population)
    population.get_chromosomes().sort(key=lambda x:x.get_fitness(), reverse=False)
    plt.scatter(generation_number,population.get_chromosomes()[0].get_fitness() )
    _print_population(population,generation_number)
    plt.pause(0.01)
    generation_number+=1
print('Minmization -- at the function :',population.get_chromosomes()[0]._xRealValues)
plt.show()


    

