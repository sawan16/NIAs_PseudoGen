
import random
import os


def func3(dna):
    f  = open('tune/travatar.ini','w')
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
    b = float(score.split('\t')[0].split('=')[1])
    r = float(score.split('\t')[1].split('=')[1])
    os.system("java -Xmx2G -jar tools/meteor-1.5/meteor-*.jar hyp_dev.anno dev.entok -norm -writeAlignments >result.txt")
    f = open('result.txt')
    score = f.readlines()
    #print(score[-1].split('            ')[1][:-1])
    f.close()
    m= float(score[-1].split('            ')[1][:-1])
    print(b,r,m)
    return b*0.5+r*0.3+m*0.2
    



def ensure_bounds(vec, bounds):

    vec_new = []
    
    for i in range(len(vec)):

        
        if vec[i] < bounds[i][0]:
            vec_new.append(bounds[i][0])

        
        if vec[i] > bounds[i][1]:
            vec_new.append(bounds[i][1])

      
        if bounds[i][0] <= vec[i] <= bounds[i][1]:
            vec_new.append(vec[i])
        
    return vec_new




def main(cost_func, bounds, popsize, mutate, recombination, maxiter):
    best_fl = open('best_bleu.txt','w')
    
    
    population = []
    for i in range(0,popsize):
        indv = []
        for j in range(len(bounds)):
            indv.append(random.uniform(bounds[j][0],bounds[j][1]))
        population.append(indv)
    print(population)

    for i in range(1,maxiter+1):
        print ('GENERATION:',i)

        gen_scores = [] # score keeping

        
        for j in range(0, popsize):

          
            canidates = list(range(0,popsize))
            canidates.remove(j)
            random_index = random.sample(canidates, 3)

            x_1 = population[random_index[0]]
            x_2 = population[random_index[1]]
            x_3 = population[random_index[2]]
            x_t = population[j]    

           
            x_diff = [x_2_i - x_3_i for x_2_i, x_3_i in zip(x_2, x_3)]

            
            v_donor = [x_1_i + mutate * x_diff_i for x_1_i, x_diff_i in zip(x_1, x_diff)]
            v_donor = ensure_bounds(v_donor, bounds)

            
            v_trial = []
            for k in range(len(x_t)):
                crossover = random.random()
                if crossover <= recombination:
                    v_trial.append(v_donor[k])

                else:
                    v_trial.append(x_t[k])
                    
           

            score_trial  = cost_func(v_trial)
            score_target = cost_func(x_t)

            if score_trial > score_target:
                population[j] = v_trial
                gen_scores.append(score_trial)
                print ('   >',score_trial, v_trial)

            else:
                print ('   >',score_target, x_t)
                gen_scores.append(score_target)

       
        gen_avg = sum(gen_scores) / popsize                         # current generation avg. fitness
        gen_best = max(gen_scores)                                  # fitness of best individual
        gen_sol = population[gen_scores.index(max(gen_scores))]     # solution of best individual
        best_fl.write(str(gen_best))
        best_fl.write('\n')
        print ('      > GENERATION AVERAGE:',gen_avg)
        print ('      > GENERATION BEST:',gen_best)
        print ('         > BEST SOLUTION:',gen_sol,'\n')
    best_fl.flush()
    best_fl.close()    
    return gen_sol



cost_func = func3                  # Cost function
bounds = [(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01)]            # Bounds [(x1_min, x1_max), (x2_min, x2_max),...]
popsize = 20                       # Population size, must be >= 4
mutate = 0.5                        # Mutation factor [0,2]
recombination = 0.95                # Recombination rate [0,1]
maxiter = 30                        # Max number of generations (maxiter)


main(cost_func, bounds, popsize, mutate, recombination, maxiter)

