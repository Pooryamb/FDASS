from numpy import log

def PreprocessingLabeledData(df):
    """This is for calculating the pfam coverage, and replacing 0 e-values, it will also change
    the order of the columns. Before running this, the input must have been labeled. """
    df["FracOfPf"] = (df["tend"] -df["tstart"] + 1)/df["tlen"]
    df["p_evalue"] = -log(df["evalue"] + 1e-300)
    df["MismatchRatio"] = df["mismatch"]/df['alnlen'] *100
    df["GapOpenRatio"] = df["gapopen"]/df['alnlen'] *100
    ColumnNames = ['fident', 'alnlen', 'MismatchRatio', 'GapOpenRatio', 'tlen', 'bits',
       'alntmscore','lddt', "FracOfPf","p_evalue", 'query', 'target','qstart',
       'qend', 'tstart', 'tend', 'qlen','PredPF', 'PF', 'PFstart', 'PFend', 'evalue', 'Status']
    df = df[ColumnNames]
    return df