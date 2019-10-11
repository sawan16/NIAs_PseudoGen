


import random
from scipy.stats import truncnorm
import os

DNA_SIZE    = 10
POP_SIZE    = 20
GENERATIONS = 30



def weighted_choice(items):

  weight_total = sum((item[1] for item in items))
  #print(weight_total)
  #print('-----')
  n = random.uniform(0, weight_total)
  #print(n)
  for item, weight in items:
    if n < weight:
      return item
    n = n - weight
  return item

def random_data():

  return chr(int(random.randrange(32, 126, 1)))

def random_population():
  
  pop = []
  '''
  for i in range(POP_SIZE):
    dna = ""
    for c in range(DNA_SIZE):
      dna += random_data()
    pop.append(dna)
  return pop
  '''
  for i in range(POP_SIZE):
    pop.append(truncnorm(a=-.01, b=.01).rvs(size=10))
  return pop

def fitness(dna):

  f  = open('tune/travatar.ini','w')
  #print('in fitness:',dna)
  data = ['[tm_file]', 'travatar-model/model/rule-table.gz', '', '[lm_file]', 'lm/lm.blm', '', '[binarize]', 'right', '', '[weight_vals]', ]
  weight_name = ['lfreq', 'fgel','egfl','egfp','lm','p', 'w','fgep','parse','unk']
  
  for i in range(10):
   data.append(weight_name[i]+'='+str(float(dna[i])))

  for i in range(len(data)):
   f.write(data[i]+'\n')

  f.flush()
  f.close()


  os.system("./run-pseudogen.sh -f tune/travatar.ini<./dev.reducedtree 2>>hyp_dev.anno")
  os.system("./test-pseudogen.sh -r dev.entok -h hyp_dev.anno >result.txt")


  f = open('result.txt')
  score = f.readline()
  f.close()
  return float(score.split('\t')[0].split('=')[1])


def mutate(dna):

  dna_out = []
  mutation_chance = 100
  for c in range(DNA_SIZE):
    if random.random() <0.1 :
      dna_out.append(float(truncnorm(a=-.01, b=.01).rvs(size=1)))
    else:
      dna_out.append(float(dna[c]))
  return dna_out

def crossover(dna1, dna2):

  dna1 = list(dna1)
  dna2 = list(dna2)
  pos = int(random.random()*DNA_SIZE)
  return (dna1[:pos]+dna2[pos:], dna2[:pos]+dna1[pos:])



if __name__ == "__main__":
 
 population = random_population()
 #print(type(population[4]))
 
 
 
 for generation in range(GENERATIONS):
  print(generation)
  weighted_population = []
  
  for individual in population:

   fitness_val = fitness(individual)
   pair = (individual, fitness_val)
   weighted_population.append(pair)
   
  population = []
  print(weighted_population,'\n')
  for _ in range(int(POP_SIZE/2)):
   # Selection
   ind1 = weighted_choice(weighted_population)
   ind2 = weighted_choice(weighted_population)
   ind1, ind2 = crossover(ind1, ind2)
   
   #print('ind1 and ind2:',ind1,ind2,'\n')
   population.append(mutate(ind1))
   population.append(mutate(ind2))

 fittest_string = population[0]
 minimum_fitness = fitness(population[0])
 
 for individual in population:
  ind_fitness = fitness(individual)
  if ind_fitness > minimum_fitness:
   fittest_string = individual
   minimum_fitness = ind_fitness
 
 print ("Fittest String: %s" % fittest_string)
 exit(0)
  