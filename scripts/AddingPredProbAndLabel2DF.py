import pandas as pd
def AddingPredProbAndLabel2DF(X, model):
    """This is for adding the predictions and the confidence score to each alignment record"""
    PredProb  = model.predict_proba(X[:,:9])[:, 1]
    PredLabel = model.predict(X[:,:9])
    df = pd.DataFrame(data = X, columns = ['fident', 'alnlen', 'MismatchRatio', 'GapOpenRatio', 'tlen', 'bits',
       'alntmscore', 'lddt', 'FracOfPf', 'p_evalue','query','qstart', 'qend','PredPF', 'PF',"Status"])
    df["PredLabel"] = PredLabel
    df["RF_prob"] = PredProb
    df = df[df["PredLabel"]==1]
    df.sort_values(by=["query", "RF_prob"], ascending=[True, False], ignore_index = True, inplace=True)
    df = df[['query','PredPF','qstart', 'qend',"RF_prob", "PredLabel",  
        'fident', 'alnlen', 'MismatchRatio', 'GapOpenRatio', 'tlen', 'bits',
        'alntmscore', "FracOfPf","p_evalue", 'PF',"Status"]]
    return df