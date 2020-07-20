import numpy as np
class TestResult():



    def __init__(self,rho0,rhoF,fVal,fid,inten,count,cTalk,errorMes = ""):
        self.startingRho = rho0
        self.myDensitie = rhoF
        self.myfVal = fVal
        self.myFidel = fid
        self.intensity = inten
        self.totalCountActual = count
        self.crossTalkMatrix = cTalk
        self.errorTraceBack = errorMes
    # def __eq__(self, other):
    #     if not isinstance(other, Test):
    #         # don't attempt to compare against unrelated types
    #         return NotImplemented
    #     return self.foo == other.foo and self.bar == other.bar

    def __str__(self):

        # Actual
        FORREPLACE = '<tr>'
        FORREPLACE += '<td colspan="2">'
        FORREPLACE += printMatrix(self.startingRho)
        FORREPLACE += '</td>'
        FORREPLACE += '<td>'
        FORREPLACE += '<div class="title">Counts: </div>' + str(int(self.totalCountActual))
        FORREPLACE += '</td>'

        if(self.myfVal == -1):
            FORREPLACE += '<td colspan="3">'
            FORREPLACE += '<b style="color:red;">Error</b>'
            FORREPLACE += '<p>'+self.errorTraceBack+'</p>'
            FORREPLACE += '<a href="">Click for eval file</a>'
            FORREPLACE += '</td>'
        else:
            FORREPLACE += '<td>'
            FORREPLACE += printfValsFidelity(self.myfVal, self.myFidel,self.intensity)
            FORREPLACE += '</td>'
            FORREPLACE += '<td colspan="2">'
            FORREPLACE += printMatrix(self.myDensitie)
            FORREPLACE += '</td>'

            # # Calculated
            # FORREPLACE += '<td>'
            # FORREPLACE += printfValsFidelity(self.myfVal, self.myFidel,self.intensity)
            # FORREPLACE += '</td>'
            # FORREPLACE += '<td colspan="2">'
            # FORREPLACE += printMatrix(self.myDensitie)
            # FORREPLACE += '</td>'

        # Ctalk
        if not isinstance(self.crossTalkMatrix,int):
            FORREPLACE += '<td colspan="2">'
            FORREPLACE += printMatrix(self.crossTalkMatrix)
            FORREPLACE += '</td>'


        FORREPLACE += '</tr>'
        return FORREPLACE
def printMatrix(M):
    s = np.shape(M)
    res = str('<table class="densityMat">')
    for i in range(s[0]):
        res += ' <tr>'
        for j in range(s[1]):
            res += '<td>' + str(round(M[i][j].real, 6)) + '<div class="purp">+</div><BR>' + str(
                round(M[i][j].imag, 6))
            res += '<div class="purp">j</div></td>'
        res += '</tr>'
    res += '</table>'
    d, v = np.linalg.eig(M)
    colorR = 0
    if (any(d < 0) or any(d > 1)):
        colorR = 255
    eigenVals = '<h5 class="eigVals" style="color:rgb(' + str(
        colorR) + ',0,0);">Eigen Values: </h5><p class="eigVals">'
    for x in range(0, len(d)):
        eigenVals += str(round(d[x].real, 5))
        if (abs(d[x].imag) > .00001):
            eigenVals += '<div class="purp">+</div>'
            eigenVals += str(round(d[x].imag, 5))
            eigenVals += '<div class="purp">j</div>'
        eigenVals += " , "
    eigenVals = str(eigenVals)[0:len(str(eigenVals)) - 2]
    eigenVals += "<p>"
    res += eigenVals
    return res


def printfValsFidelity(fVal, fid, i):
    colorR = mapFloat(fid, 1, .9, 0, 255)  # red
    res = '<ul class="fVals">'

    # fval
    res += '<li>'
    if (fVal == -1):
        res += "ERROR THIS SHOW NOT SHOW"
    else:
        res += '<div class="title">fVal: </div>'
        res += str(round(fVal, 6))
    res += '</li>'

    # fidelity
    res += '<li style="color:rgb(' + str(colorR) + ',0,0);">'
    res += '<div class="title">Fidelity: </div><div class="inline">'
    res += str(round(fid, 6))
    res += '</div></li>'

    # intensity
    res += '<li>'
    res += 'Intensity: ' + str(i)
    res += '</li>'

    res += str('</ul>')
    return res
# Maps a range of numbers to another range of numbers
def mapFloat(value, x1, x2, y1, y2, f=1):
    newVal = value - x1

    newVal = (y2 - y1) / (x2 - x1) * newVal
    newVal += y1
    return newVal