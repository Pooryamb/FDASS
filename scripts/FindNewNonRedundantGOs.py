from goatools.base import download_go_basic_obo
import pandas as pd
from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag

fin_dag = download_go_basic_obo("go-basic.obo")
GAF = GafReader('../rawinput/Tb/TriTrypDB-60_TbruceiTREU927_GO.gaf')
InDictFormat = GAF.get_ns2assc()

godag = GODag(fin_dag, optional_attrs={'consider', 'replaced_by'}, load_obsolete=True)

PfAnns = pd.read_csv("../predictions/NewPredictionsTb.txt", sep="\t")
PF2Go = pd.read_csv("../rawinput/pfam2go.txt", sep="\t")
mappings = pd.read_csv("../rawinput/Tb/MappingFromUnipToTriTrypID.txt", sep="\t", header=None)

PfAnns = PfAnns[["Query", "PredPF"]]

PfAnnsTri = PfAnns.merge(mappings, left_on="Query", right_on=1)
PfAnnsTri = PfAnnsTri[["Query", "PredPF", 0]]

PredGo = PfAnnsTri.merge(PF2Go, left_on = "PredPF", right_on="PF")

PredGo.columns = ['UnipID', 'PredPF', "TriTrypID", 'PF', 'GO', 'GO_desc']
PredGo = PredGo[["TriTrypID", 'UnipID', 'PredPF', 'GO', 'GO_desc']]

Associations = {}
for ns in InDictFormat.keys():
    for gene in InDictFormat[ns].keys(): 
        #SetOfGOsOfGene will store the GOs and all their parents
        SetOfGOsOfGene = Associations.get(gene, set())
        PrevGOs = InDictFormat[ns][gene]
        for GO in PrevGOs:
            if godag[GO].is_obsolete and godag[GO].replaced_by != "":
                GO = godag[GO].replaced_by
            GOs2Add = godag[GO].get_all_parents()
            GOs2Add.add(GO)
            SetOfGOsOfGene = SetOfGOsOfGene.union(GOs2Add)
        Associations[gene] = SetOfGOsOfGene


NewAnnotations = []

for i in range(PredGo.shape[0]):
    line = PredGo.iloc[i]
    if line["GO"] not in Associations.get(line["TriTrypID"],set()):
        NewAnnotations.append(line)

NewGO_ann = pd.DataFrame(NewAnnotations)

NewGO_ann = NewGO_ann[["TriTrypID", "UnipID", "GO",'GO_desc']].drop_duplicates()
NewGO_ann.to_csv("../predictions/NewlyPredictedGOs.tsv", sep="\t", index=None)
