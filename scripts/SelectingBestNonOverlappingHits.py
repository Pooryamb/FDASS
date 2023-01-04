import pandas as pd
from OverLapReturner import *
def SelectingBestNonOverlappingHits(df,PF_Anns = {}):
    """This function takes the dataframe containing all preprocessed alignments as input, and reports 
    the highest scoring non-overlapping alignments on each target protein. Optionally, you can specify 
    the current Pfam annotations of the proteins in a dictionary whose keys are gene names, and values
    are tuples of the intervals of the pfams"""
    ListOfRows = []
    OverlapThreshold = 10
    for (i, CurID,start, end)  in zip(df.index, df["query"], df['qstart'],df['qend']):
        row = df.loc[i,:]
        ListOfIntervals = PF_Anns.get(CurID, [])
        SignificantOverLap = False
        for Interval in ListOfIntervals:
            if OverLapReturner(Interval, (start,end)) >= OverlapThreshold:
                SignificantOverLap = True
                break
        if (SignificantOverLap==False):
            ListOfRows.append(row)
            ListOfIntervals.append((start,end))
            PF_Anns[CurID] = ListOfIntervals
    OutDF = pd.DataFrame(data= ListOfRows, columns = df.columns)
    return OutDF