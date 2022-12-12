import pandas as pd


Coverage_Of_PfamChoppedStruct = 0.8
Coverage_Of_PfamOnQueryPosThresh = 0.5
FalsePositiveThreshold = 0.25

##### columns in FS search
# query 0 # target 1 # fident 2 # alnlen 3 # mismatch 4 # gapopen 5 # qstart 6 # qend 7 # tstart 8 # tend 9 # qlen 10 # tlen 11 # evalue 12 # bits 13 # alntmscore 14

def PreparingAliDF(AliRawDf):
    AliRawDf.columns = "query target fident alnlen mismatch gapopen qstart qend tstart tend qlen tlen evalue bits alntmscore".split()
    AliRawDf["query"] = AliRawDf["query"].str.replace(".pdb.gz","", regex=False)
    AliRawDf["target"] = AliRawDf["target"].str.replace("-F1-model_v3.pdb.gz","", regex=False).str.replace("-F1-model_v3.cif.gz","", regex=False).str.replace("AF-",'', regex=False)
    AliRawDf["PredPF"] = AliRawDf["query"].str.split("_", expand=True)[3]
    return AliRawDf

def PreparingPfamDF(PF_map):
    PF_map.columns = ["GeneID","PF","PFstart","PFend"]
    return PF_map
    


def FindingLabels(FoldSeekRawOutput, PfamFile, FoldSeekLabeledOutput): #it takes the address of the pfam as input   
    AliRawDf =  pd.read_csv(FoldSeekRawOutput,sep="\t", header=None)
    PF_map = pd.read_csv(PfamFile, sep="\t", header=None)
    PF_map = PreparingPfamDF(PF_map)
   
    FormattedAli = PreparingAliDF(AliRawDf)
    del AliRawDf 
    merged_df = pd.merge(FormattedAli, PF_map, left_on = ["target"], right_on = ["GeneID"], how='left')
    merged_df["max_start"] = merged_df[["tstart","PFstart"]].max(axis=1)
    merged_df["min_end"] = merged_df[["tend","PFend"]].min(axis=1)
    merged_df["Status"] = -1
    
    TruePosInds = (merged_df["PredPF"] == merged_df["PF"]) & ((merged_df["min_end"] - merged_df["max_start"] +1 )/(merged_df['PFend'] - merged_df['PFstart'] +1) > Coverage_Of_PfamOnQueryPosThresh)    
    merged_df.loc[TruePosInds, "Status"] = 1

    FalsPosInds = (merged_df["PredPF"] != merged_df["PF"]) & ((merged_df["min_end"] - merged_df["max_start"] +1 )/(merged_df['tend'] - merged_df['tstart'] +1 ) > FalsePositiveThreshold)
    merged_df.loc[FalsPosInds, "Status"] = 0
    
    merged_df.loc[pd.isnull(merged_df["PF"]),"Status" ] = -1  
    #As I am doing left merge, I have to set the status of cases where the protein has no 
    # pfam annotation to -1. If I didn't include this step, they would have been labeled as 0
    
    merged_df2 = merged_df.sort_values("Status", ascending=False).drop_duplicates(['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart',
       'qend', 'tstart', 'tend', 'qlen', 'tlen', 'evalue', 'bits','alntmscore'])
    merged_df2 = merged_df2[['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'qlen', 'tlen', 'evalue', 'bits', 
    'alntmscore', 'PredPF', 'PF', 'PFstart', 'PFend', 'Status']]
    #LablelledOutputAdd = FoldSeekOutput + "_labeled.tsv"
    merged_df2.to_csv(FoldSeekLabeledOutput , sep="\t", index=None)
    return 0
    
    
RawFSoutputFormat =  "../rawinput/{}/result_aln_pf_{}_tm"
PfamAddFormat = "../rawinput/{}/Pfam{}.txt"
LabeledData = "../intermediates/{}/result_aln_pf_{}_tm_labeled.tsv"

Orgs = ["Tb", "Mj", "Ec", "Sc"]

for org in Orgs:
    FindingLabels(RawFSoutputFormat.format(org, org), PfamAddFormat.format(org, org), LabeledData.format(org, org))
