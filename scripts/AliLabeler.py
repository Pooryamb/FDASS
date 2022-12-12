import pandas as pd


Coverage_Of_PfamChoppedStruct = 0.8
Coverage_Of_PfamOnQueryPosThresh = 0.5
FalsePositiveThreshold = 0.25

##### columns in FS search
# query 0 # target 1 # fident 2 # alnlen 3 # mismatch 4 # gapopen 5 # qstart 6 # qend 7 # tstart 8 # tend 9 # qlen 10 # tlen 11 # evalue 12 # bits 13 # alntmscore 14

def FindingLabels(FoldSeekOutput, PfamFile): #it takes the address of the pfam as input   
    PF_df =  pd.read_csv(FoldSeekOutput,sep="\t", header=None)
    PF_map = pd.read_csv(PfamFile, sep="\t", header=None)
    PF_map.columns = ["GeneID","PF","PFstart","PFend"]
    
    PF_df.columns = "query target fident alnlen mismatch gapopen qstart qend tstart tend qlen tlen evalue bits alntmscore".split()
    PF_df["query"] = PF_df["query"].str.replace(".pdb.gz","")
    PF_df["target"] = PF_df["target"].str.replace("-F1-model_v3.pdb.gz","").str.replace("-F1-model_v3.cif.gz","").str.replace("AF-",'')
    PF_df["PredPF"] = PF_df["query"].str.split("_", expand=True)[3]
    
    merged_df = pd.merge(PF_df, PF_map, left_on = ["target"], right_on = ["GeneID"], how='left')
    merged_df["max_start"] = merged_df[["tstart","PFstart"]].max(axis=1)
    merged_df["min_end"] = merged_df[["tend","PFend"]].min(axis=1)
    merged_df["Status"] = -1
    
    TruePosInds = (merged_df["PredPF"] == merged_df["PF"]) & ((merged_df["min_end"] - merged_df["max_start"] +1 )/(merged_df['PFend'] - merged_df['PFstart'] +1) > Coverage_Of_PfamOnQueryPosThresh)    
    merged_df.loc[TruePosInds, "Status"] = 1

    FalsPosInds = (merged_df["PredPF"] != merged_df["PF"]) & ((merged_df["min_end"] - merged_df["max_start"] +1 )/(merged_df['tend'] - merged_df['tstart'] +1 ) > FalsePositiveThreshold)
    merged_df.loc[FalsPosInds, "Status"] = 0
    
    merged_df.loc[pd.isnull(merged_df["PF"]),"Status" ] = -1

    merged_df2 = merged_df.sort_values("Status", ascending=False).drop_duplicates(['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart',
       'qend', 'tstart', 'tend', 'qlen', 'tlen', 'evalue', 'bits','alntmscore'])
    merged_df2 = merged_df2[['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'qlen', 'tlen', 'evalue', 'bits', 
    'alntmscore', 'PredPF', 'PF', 'PFstart', 'PFend', 'Status']]
    LablelledOutputAdd = FoldSeekOutput + "_labeled.tsv"
    merged_df2.to_csv(LablelledOutputAdd , sep="\t", index=None)
    return 0
    
    
RawFSoutputFormat =  "D:/McGillThesis/WholeGenomeStructuralBasedAnnotation/{}AgPfam/result_aln_pf_{}_tm"
PfamAddFormat = "D:/McGillThesis/WholeGenomeStructuralBasedAnnotation/{}AgPfam/Pfam{}.txt"

Orgs = ["Sm"] 
#["Tb", "Mj", "Ec", "Sc"]

for org in Orgs:
    FindingLabels(RawFSoutputFormat.format(org, org), PfamAddFormat.format(org, org))
