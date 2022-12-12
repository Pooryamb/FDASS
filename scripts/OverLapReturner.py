def OverLapReturner(Interval1, Interval2):
    """It takes two intervals as input and returns the overlap between them as output"""
    MaxStart = max(Interval1[0], Interval2[0])
    MinEnd  =  min(Interval1[1], Interval2[1])
    return max(MinEnd-MaxStart+1,0)