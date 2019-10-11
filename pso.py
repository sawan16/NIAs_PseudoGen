from __future__ import division
import random
from scipy.stats import truncnorm
import os

def func2(dna):
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
    os.system("./test-pseudogen.sh -r dev.entok -h hyp_dev.anno >result4.txt")
    f = open('result4.txt')
    score = f.readline()
    f.close()
    b = float(score.split('\t')[0].split('=')[1])
    r = float(score.split('\t')[1].split('=')[1])
    os.system("java -Xmx2G -jar tools/meteor-1.5/meteor-*.jar hyp_dev.anno dev.entok -norm -writeAlignments >result4.txt")
    f = open('result4.txt')
    score = f.readlines()
    #print(score[-1].split('            ')[1][:-1])
    f.close()
    m= float(score[-1].split('            ')[1][:-1])
    print(b,r,m)
    return b*0.5+r*0.3+m*0.2


class Particle:
    def __init__(self,x0):
        self.position_i=[]          # particle position
        self.velocity_i=[]          # particle velocity
        self.pos_best_i=[]          # best position individual
        self.err_best_i=0.0          # best error individual
        self.err_i=0.0             # error individual

        for i in range(0,num_dimensions):
            self.velocity_i.append(random.uniform(-0.01,0.01))
            self.position_i.append(x0[i])


    def evaluate(self,costFunc):
        self.err_i=costFunc(self.position_i)
        print('position:',self.position_i,'velocity:',self.velocity_i,'score:',self.err_i)
        
        if self.err_i > self.err_best_i or self.err_best_i==1.0:
            self.pos_best_i=self.position_i
            self.err_best_i=self.err_i

    
    def update_velocity(self,pos_best_g):
        w=0.65      # constant inertia weight (how much to weigh the previous velocity)
        c1=1.63        # cognative constant
        c2=0.62        # social constant

        for i in range(0,num_dimensions):
            r1=random.random()
            r2=random.random()

            vel_cognitive=c1*r1*(self.pos_best_i[i]-self.position_i[i])
            vel_social=c2*r2*(pos_best_g[i]-self.position_i[i])
            self.velocity_i[i]=w*self.velocity_i[i]+vel_cognitive+vel_social

   
    def update_position(self,bounds):
        for i in range(0,num_dimensions):
            self.position_i[i]=self.position_i[i]+self.velocity_i[i]

            
            if self.position_i[i]>bounds[i][1]:
                self.position_i[i]=bounds[i][1]

            
            if self.position_i[i] < bounds[i][0]:
                self.position_i[i]=bounds[i][0]
                
class PSO():
    def __init__(self,costFunc,x0,bounds,num_particles,maxiter):
        global num_dimensions

        num_dimensions=len(x0)
        err_best_g=0.0                  # best error for group
        pos_best_g=[]                   # best position for group

       
        swarm=[]
        for i in range(0,num_particles):
            swarm.append(Particle(x0))

      
        i=0
        while i < maxiter:
            print(i)
            
            for j in range(0,num_particles):
                swarm[j].evaluate(costFunc)
                
                if swarm[j].err_i > err_best_g or err_best_g == 1.0:
                    pos_best_g=list(swarm[j].position_i)
                    err_best_g=float(swarm[j].err_i)

            
            for j in range(0,num_particles):
                swarm[j].update_velocity(pos_best_g)
                swarm[j].update_position(bounds)
            i+=1

        
        print ('FINAL:')
        print (pos_best_g)
        print (err_best_g)



initial=truncnorm(a=-.01, b=.01).rvs(size=10) # initial starting location [x1,x2...]
bounds = [(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01),(-0.01,0.01)]  # input bounds [(x1_min,x1_max),(x2_min,x2_max)...]
PSO(func2,initial,bounds,num_particles=20,maxiter=30)

