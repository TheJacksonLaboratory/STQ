"""
Created on Feb 14 2017 by Tony Szedlak
Modified on Dec 06 2018 by Sergii Domanskyi
"""

import os
import copy
import numpy as np
import matplotlib
if os.name != 'nt': matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocess as mp
import platform
import time

from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
from pandas.errors import PerformanceWarning
simplefilter("ignore", category=ConvergenceWarning)
simplefilter("ignore", category=PerformanceWarning)

class partial:

    """Return "partially evaluated" version of given function/arguments.
    This is similar to functools partial"""

    def __init__(self, func, *args, **kwargs):

        self.func, self.args, self.kwargs = func, args, kwargs

    def __call__(self, *more_args, **more_kwargs):

        all_kwargs = {**self.kwargs, **more_kwargs}

        return self.func(*self.args, *more_args, **all_kwargs)

class Agent:
    
    def __init__(self, pool, numTargets, mutationProb, genome=None, ObjFunction=None):

        self.mutationProb = mutationProb
        self.pool = pool
        self.numTargets = numTargets
        self.ObjFunction = ObjFunction

        if genome is not None:
            self.genome = copy.deepcopy(genome)
        else:
            self.RandomAgent()

        self.genome = np.sort(self.genome).tolist()

    def Mutate(self):

        for i in range(len(self.genome)):
            if np.random.rand() < self.mutationProb:
                availablePool = np.setdiff1d(self.pool, self.genome[:i] + self.genome[i + 1:]).tolist()
                self.genome[i] = np.random.choice(availablePool).tolist()
        
        self.genome = np.sort(self.genome).tolist()

    def RandomAgent(self):

        def GetOneRandom():
            return np.random.choice(self.pool, size=self.numTargets, replace=False).tolist()

        if not self.ObjFunction is None:
            bestRandomGenome = GetOneRandom()
            bestRandomScore = self.ObjFunction(bestRandomGenome)
            for i in range(1): # 10**2
                genome = GetOneRandom()
                score = self.ObjFunction(genome)
                if score > bestRandomScore:
                    bestRandomGenome, bestRandomScore = genome, score
            self.genome = bestRandomGenome
        else:
            self.genome = GetOneRandom()

        return

def MateAgents(agentList, scoreList, fractionToMate):

    def Mate(agent1, agent2):

        availableGenes = agent1.genome + agent2.genome
        uniqueGenes = np.unique(availableGenes).tolist()
        genesAndScoresDict = {i:0.0 for i in uniqueGenes}

        for i in availableGenes:
            genesAndScoresDict[i] = genesAndScoresDict[i] + 1.0

        numHits = [i[1] for i in genesAndScoresDict.items()]
        probs = numHits / np.sum(numHits)
        outputGenome = np.random.choice(uniqueGenes, size=len(agent1.genome), p=probs, replace=False).tolist()

        return outputGenome

    bestToWorst = np.argsort(scoreList)[::-1]
    whichToMate = bestToWorst[:int(np.ceil(len(bestToWorst) * fractionToMate))]
    parentPool = [agentList[i] for i in whichToMate]

    newAgentList = copy.deepcopy(parentPool)

    def chooseParents(parentPool):

        np.random.seed()

        return np.random.choice(parentPool, size=2, replace=False)

    while len(newAgentList) < len(agentList):
        bestParents, bestSimilarity = chooseParents(parentPool), 1
        for i in range(1):
            parents = chooseParents(parentPool)
            similarity = len(np.intersect1d(parents[0].genome, parents[1].genome)) / len(parents[0].genome)
            if similarity < bestSimilarity:
                bestSimilarity, bestParents = similarity, parents

        parents = bestParents
        newAgent = Agent(parents[0].pool,
                        parents[0].numTargets,
                        mutationProb = parents[0].mutationProb,
                        genome = Mate(parents[0],parents[1]),)

        newAgent.Mutate()

        newAgentList.append(newAgent)

    return newAgentList

def Maximize(ObjFunction, pool, numTargets, numIterations, numAgents, startingGenome=None, fractionToMate=0.3, mutationProb=0.3, saveFigure=True, parallel=True,
             trackScores=False, saveFreq=10, verbose=False, saveDir='../genAlgPlots/', NumberOfAvailableCPUs=4, keepOnlyBestLiveResults=False, saveGenomeAtEachIteration=False):

    agentList = [Agent(pool, numTargets, mutationProb, genome=startingGenome, ObjFunction=None) for i in range(numAgents)]
    print('Set of %s agents ready' % numAgents, flush=True)

    if parallel:
        pool = mp.Pool(processes=NumberOfAvailableCPUs)
    print('Pool of %s workers ready' % NumberOfAvailableCPUs, flush=True)

    if trackScores:
        allScores = []
        fig, ax = plt.subplots(2, 1, figsize=(6, 6))

        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
    
    iterable = range(numIterations)
    allBestResults = [[None,None]]
    
    for iteration in iterable:
        print('Iteration:', iteration, end='\t', flush=True)
        if parallel:
            scoreList = np.vstack(pool.map(ObjFunction, [agent.genome for agent in agentList])).squeeze().tolist()
        else:
            scoreList = [ObjFunction(agent.genome) for agent in agentList]

        ## Add any agents to best results, even if such agents are already in the best results.
        #allBestResults += list(zip([agent.genome for agent in agentList], scoreList))

        # Do not add agents to best results, if such agents are already in the best results.
        allBestResults += [i for i in list(zip([agent.genome for agent in agentList], scoreList)) if i[0] not in list(zip(*allBestResults))[0]]

        if [None,None] in allBestResults:
            allBestResults.remove([None,None])

        allBestResults.sort(key=lambda x:float(x[-1]), reverse=True)

        try:
            temp = np.unique(np.array([item[0] + [item[1]] for item in allBestResults]), axis=0, return_index=True)
            temp = temp[0][np.argsort(temp[1])]
            np.savetxt('%s/bestGenomes.csv' % (saveDir), temp[:100, :], delimiter=',', fmt='%s')
        except:
            print('Best results file is unavailable')

        bestResults = allBestResults[:numAgents]

        print('Overall best: %s' % (bestResults[0][1]), '\tCurrent best: %s' % (max(scoreList)))

        if keepOnlyBestLiveResults:
            allBestResults = bestResults
        
        if trackScores:
            allScores.append(scoreList)
    
        if trackScores and iteration != 0 and (iteration == numIterations - 1 or iteration % saveFreq == 0):

            allScoresForPlot = np.vstack(allScores)
            allScoresForPlot_sorted = np.sort(allScoresForPlot, axis=1)
            meanAllScoresForPlot = np.mean(allScoresForPlot, axis=1)

            if saveFigure:
                ax[0].cla()
                ax[1].cla()

                for i in range(numAgents):
                    ax[0].plot(allScoresForPlot[:,i], '-o', linewidth=0.5, markersize=2)
                    ax[1].plot(allScoresForPlot_sorted[:,i], '-o', linewidth=0.5, markersize=2)

                ax[0].plot(meanAllScoresForPlot, 'k-', linewidth=3)
                ax[1].plot(meanAllScoresForPlot, 'k-', linewidth=3)
                ax[0].plot(meanAllScoresForPlot, 'w-o', linewidth=1, markersize=6, markeredgecolor='k', markeredgewidth=1)
                ax[1].plot(meanAllScoresForPlot, 'w-o', linewidth=1, markersize=6, markeredgecolor='k', markeredgewidth=1)

                ax[0].set_xlabel('Generation index', fontsize=14)
                ax[1].set_xlabel('Generation index', fontsize=14)

                ax[0].set_ylabel('Agent scores', fontsize=14)
                ax[1].set_ylabel('Agent ranked scores', fontsize=14)

                fig.tight_layout()
            
                try:
                    fig.savefig('%s/geneticAlgorithm.png' % (saveDir), dpi=300)
                except:
                    print('Figure file is unavailable')

            if saveGenomeAtEachIteration:
                try:
                    np.savetxt('%s/genomeAtIteration%s.csv'%(saveDir, iteration), [agent.genome + [score] for agent, score in zip(agentList, scoreList)], delimiter=',', fmt='%s')
                except:
                    print('CSV file is unavailable')

        if iteration != numIterations - 1:
            agentList = MateAgents(agentList, scoreList, fractionToMate)

    if parallel:
        pool.close()
        pool.join()

    try:
        temp = np.unique(np.array([item[0] + [item[1]] for item in allBestResults]), axis=0, return_index=True)
        temp = temp[0][np.argsort(temp[1])]
        np.savetxt('%s/allBestGenomes.csv' % (saveDir), temp[:, :], delimiter=',', fmt='%s')
    except:
        print('All best results file is unavailable')

    return bestResults

if __name__ == '__main__':
    
    #=================================================
    #     Test to demonstrate Genetic Algorithm #
    #=================================================

    np.random.seed(42)

    def objFunction(x): 
        return np.sum(x)

    pool = range(10000)
    numTargets = 20
    
    bestResults = Maximize(objFunction, 
                           pool = pool, 
                           numTargets = numTargets, 
                           numIterations = 100, 
                           numAgents = 40, 
                           fractionToMate = 0.3, 
                           mutationProb = 0.3, 
                           trackScores=True, 
                           saveDir='genAlgPlots/', 
                           saveFreq=100)
    
    for agent_genome,score in bestResults:
        print(np.sort(agent_genome),' || ',score)
   
    scoreList = list(tuple(zip(*bestResults))[1])
    print('Mean rel. score = ', np.mean(scoreList) / objFunction(pool[-numTargets:]))
