from supplementary import *
from GeneticAlgorithm import *

if __name__ == '__main__':

    # Fix the random sequence starting seed for reproducibility
    np.random.seed(142)

    drug_case = 'JAX PDX Docetaxel logTMM'

    # genes expressed in less than specified fraction of training samples (e.g. 50%) are filtered out
    minForGene = 0.95 

    # Prepare dir for results
    wdir = 'PDX search 01 07 2022 ttest-selected ridge %s %s/' % (drug_case, minForGene)
    mdir(wdir)

    print(wdir)

    if not os.path.isfile(wdir + 'df_train.h5'):
        # Load the data: The loaded df_expr should be a pandas.DataFrame:
        # Genes identifiers are in rows.
        # Values are normalized gene expression.
        # Columns should be pandas.MultiIndex and have at least 2 levels. There should be level named "response", containing True/False/NA.
        # One could designate False for response PD/SD, and True for PR/CR.
        #
        # Data table may look like:
        #sample_id    241        242        249        252        ...
        #response     NA         True       NA         False      ...
        #SCYL3        16.230947  15.589821  16.713638  18.227494  ...
        #FGR          12.703414  10.696765  12.109998  11.872076  ...
        #RAD52        15.147076  16.381272  15.348719  15.806286  ...
        #...          ...        ...        ...        ...        ...
        #
        # Note: if one has a CSV spreadsheet in this format. Then to load it do:
        # df_expr = pd.read_csv('path/to/data.csv', index_col=0, header=[0, 1])
        try:
            df_expr = pd.read_csv('data/%s.csv.gz' % drug_case, index_col=0, header=[0, 1, 2, 3])
        except:
            df_expr = pd.read_csv('data/%s.csv.gz' % drug_case, index_col=0, header=[0, 1, 2])

        dfi = df_expr.columns.to_frame()
        dfi['response'] = dfi['response'].replace({'True': True, 'False': False}).astype(bool)
        df_expr.columns = pd.MultiIndex.from_frame(dfi)
    
        # Remove any samples with unknown (not measured) response.
        df_expr = df_expr[df_expr.columns[df_expr.columns.get_level_values('response').isin([True, False])]]

        # Save 20 non-responders and 3 responders for validation later
        df_validate = pd.concat([df_expr.xs(False, level='response', axis=1, drop_level=False).sample(20, axis=1),
                                 df_expr.xs(True,  level='response', axis=1, drop_level=False).sample(3,  axis=1)], axis=1, sort=False)
        df_train = df_expr[df_expr.columns.difference(df_validate.columns)]

        # Determine highly variable genes (hvg)
        wgenes = df_train.index[((df_train>0).mean(axis=1)>=minForGene)]
        print('Genes expressed in at at least %s%% of the training samples:' % (100*minForGene), len(wgenes))

        dfs = [df_train.loc[wgenes].xs(v, level='response', axis=1, drop_level=False) for v in [True, False]]
        hvg = ((dfs[0].mean(axis=1) - dfs[1].mean(axis=1)).abs() / np.sqrt(dfs[0].var(axis=1)**2 + dfs[1].var(axis=1)**2)/len(dfs[0])).sort_values(ascending=False)[:10000].index
        #hvg = (df_train.loc[wgenes].var(axis=1)/df_train.loc[wgenes].mean(axis=1)).sort_values(ascending=False)[:10000].index

        # Retain hvg in the data
        df_train = df_train.loc[hvg]
        df_validate = df_validate.loc[hvg]
        print('Train/test:', df_train.shape, '\tValidate:', df_validate.shape)

        # Save data to use later in the validation step
        df_train.to_hdf(wdir + 'df_train.h5', key='df', mode='a', complevel=4, complib='zlib')
        df_validate.to_hdf(wdir + 'df_validate.h5', key='df', mode='a', complevel=4, complib='zlib')
    else:
        # Load data
        df_train = pd.read_hdf(wdir + 'df_train.h5', key='df')

    # NumberOfAvailableCPUs < numAgents, or the excess resources will be idle
    if platform.system()=="Linux": # For sumner HPC runs
        numAgents = 36
        NumberOfAvailableCPUs = 36
    else: # For local machine tests
        numAgents = 128
        NumberOfAvailableCPUs = 4

    #df_train.iloc[:] = scale(df_train.T).T

    print(partial(scoreCorr, df=df_train)(df_train.index.values[100:120]))
    exit()


    # Do the search: it may take from minutes to hours depending on the data and objective function
    bestResults = Maximize(partial(scoreCorr, df=df_train), pool=df_train.index.values, 
    #bestResults = Maximize(partial(scoreCorr, df=df_train, method='ridge', scoring='f1score'), pool=df_train.index.values, 
                           fractionToMate=0.3, mutationProb=0.3, trackScores=True, 
                           numTargets=20, numIterations=1000, numAgents=numAgents, 
                           parallel=True, NumberOfAvailableCPUs=NumberOfAvailableCPUs,
                           saveDir=wdir, saveFreq=10)

    # partial(regressionSamplingObjective, df=df_train, N=100, sample=4)
