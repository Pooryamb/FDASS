import pandas as pd
def AddingPredProbAndLabel2DF(X, model):
    """This is for adding the predictions and the confidence score to each alignment record"""
    PredProb  = model.predict_proba(X[:,:9])[:, 1]
    PredLabel = model.predict(X[:,:9])
    df = pd.DataFrame(data = X, columns = ['fident', 'alnlen', 'MismatchRatio', 'GapOpenRatio', 'tlen', 'bits',
       'alntmscore','lddt', "FracOfPf","p_evalue", 'query', 'target','qstart',
       'qend', 'tstart', 'tend', 'qlen','PredPF', 'PF', 'PFstart', 'PFend', 'evalue', 'Status'])
    df["PredLabel"] = PredLabel
    df["RF_prob"] = PredProb
    df = df[df["PredLabel"]==1]
    df.sort_values(by=["query", "RF_prob"], ascending=[True, False], ignore_index = True, inplace=True)
    dfNewOrder = df[['query','PredPF','qstart', 'qend',"RF_prob", "PredLabel",  
        'fident', 'alnlen', 'MismatchRatio', 'GapOpenRatio', 'qlen', 'bits',
        'alntmscore', "FracOfPf","p_evalue", 'PF', 'PFstart', 'PFend', 'evalue', 'Status']]
    return df