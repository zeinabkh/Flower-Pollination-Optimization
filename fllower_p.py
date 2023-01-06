import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gamma,levy
import math
def Rosenbrock(X):
    # print(X,"*")
    '''INPUTS
    X: arguments of the function Rosenbrock
    OUTPUTS
    f : evaluation of the Rosenbrock function given the inputs

    DOMAIN         : Xi is within [-5,10] although can be [-2.048,2.048]
    DIMENSIONS     : any
    GLOBAL MINIMUM : f(x)=0 x=[1,...,1]
'''
    f = sum(100.0 * (X[i + 1] - X[i] ** 2) ** 2 + (1 - X[i]) ** 2 for i in range(0, len(X) - 1))
    # print(X, f)
    return f


def StyblinskiTang(X):
    '''INPUTS
    X: arguments of the Styblinski-Tang Function
    OUTPUTS
    f : evaluation of the Styblinski-Tang function given the inputs

    DOMAIN         : [-5,5]
    DIMENSIONS     : any
    GLOBAL MINIMUM : f(x)=(-39.166*d) x=[-2.9035,...,-2.9035]
    '''
    f_sum = sum((X[i] ** 4) - (16 * X[i] ** 2) + (5 * X[i]) for i in range(len(X)))
    # print(f_sum)
    return f_sum / 2


def Ackley(x):
    '''
    INPUTS
    x : arguments of the function Ackley
    Output
    f : evaluation of the Ackley function given the inputs

    DOMAIN           : [-32,32]
    DIMENSIONS       : any
    GLOBAL MINIMUM   : f(x)=0 x=[0...0]
    '''
    d = len(x)
    a = 20
    b = 0.2
    c = np.pi * 2
    sum1 = sum(x[i] ** 2 for i in range(d))
    sum1 = (-a) * np.exp(((-b) * np.sqrt(sum1 / d)))
    sum2 = sum(np.cos(c * x[i]) for i in range(d))
    sum2 = np.exp((sum2 / d))
    return sum1 - sum2 + a + np.exp(1)


class Flower_Pllination:
    def __init__(self, nop, nog, bound, dimension, p, optimum_dir, fittness_fun, l):
        self.number_of_population = nop
        self.number_of_generation = nog
        self.dimension = dimension
        self.dimension_bound = bound
        self.population = []
        self.min = optimum_dir
        self.switch_p = p
        self.fittness_fun = fittness_fun
        self.objective_values_generation = []
        self.lambda_v = l
        self.best_chrom = None

    def levy(self, d):
        lamda = 1.5
        sigma = (math.gamma(1 + lamda) * math.sin(math.pi * lamda / 2) / (
                math.gamma((1 + lamda) / 2) * lamda * (2 ** ((lamda - 1) / 2)))) ** (1 / lamda)
        # sigma = 0.6965745025576968
        u = np.random.randn(1, d) * sigma
        v = np.random.randn(1, d)
        step = u / abs(v) ** (1 / lamda)
        return 0.01 * step

    def initializing(self):
        pop = [np.array([np.random.uniform(self.dimension_bound[i][0], self.dimension_bound[i][1])\
               for i in range(self.dimension)]).squeeze() for j in range(self.number_of_population)]

        self.population = [{"xt": x,
                            "f": self.fittness_fun(x)} for x in pop]

        self.population = sorted(self.population, key=lambda x: x['f'], reverse=False)
        self.best_chrom = self.population[0]['xt']

        self.objective_values_generation.append(self.population[0]['f'])

    def flower_pollination(self):

        g = 0
        while g < self.number_of_generation:
            print(g)
            for chrom in self.population:
                p = np.random.uniform(0, 1)
                if p < self.switch_p:
                   chrom['xt1'] = chrom['xt'] + self.levy(self.dimension) * (self.best_chrom - chrom['xt'])
                else:
                    j, k = np.random.randint(0, self.number_of_population, 2)
                    chrom['xt1'] = chrom['xt'] + 0.01 * (self.population[j]['xt'] - self.population[k]['xt'])
                chrom['xt1'] = chrom['xt1'].squeeze()
            print(self.population)
            for chrom in self.population:
               if self.fittness_fun(chrom['xt']) > self.fittness_fun(chrom['xt1']):
                   chrom['xt'] = chrom['xt1']
                   chrom['f'] = self.fittness_fun(chrom['xt1'])
            self.population = sorted(self.population, key=lambda x: x['f'], reverse=False)
            self.best_chrom = self.population[0]['xt']

            self.objective_values_generation.append(self.population[0]['f'])

            g += 1

    def plot_objective_value(self):
        print(self.objective_values_generation,self.best_chrom)
        plt.plot([i for i in range(self.number_of_generation+1)], self.objective_values_generation)
        plt.show()


if __name__ == "__main__":
    p = [0.05, 0.2, 0.6]
    e = [0.5, 0.1, 0.01]
    d = 2
    flower_p = Flower_Pllination(1000, 100, [[-5, 10] for i in range(d)], d, 0.2, "min", Rosenbrock, 0.5)
    flower_p.initializing()
    flower_p.flower_pollination()
    flower_p.plot_objective_value()
    # flower_p = Flower_Pllination(1000, 1000, [[-5, 5] for i in range(d)], d, 0.5, "min", StyblinskiTang, 0.5)
    # flower_p.initializing()
    # flower_p.flower_pollination()
    # flower_p.plot_objective_value()
    # flower_p = Flower_Pllination(1000, 1000, [[-32, 32] for i in range(d)], d, 0.5, "min", Ackley, 0.5)
    # flower_p.initializing()
    # flower_p.flower_pollination()
    # flower_p.plot_objective_value()
