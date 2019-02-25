# -*- coding: utf-8 -*-
"""
Project:    Verkefni I - Kinematics of the Stewart Platform
            Numerical Analysis, 2nd ed., Sauer, Chapter 1
Course:     Toluleg Greining STAE405
Institute:  Haskoli Islands
Professor:  Sigurdur Freyr Hafstein
Authors:    Erling Oskar Kristjansson   eok4@hi.is
            David Freyr Bjornsson       dfb2@hi.is
"""

import numpy as np
import math
import scipy.optimize as scOpt
import matplotlib.pyplot as plt

""" 

Ytri breytur (stikar) sem við skilgreinum sem
global eða viðværar breytur 

"""
# global x1, x2, y2, L1, L2, L3, gamma, p1, p2, p3 
x1=y2   = 4     # tekid af fig 1.15
x2      = 0     # tekid af fig 1.15
L1      = 2     
L2 = L3 = np.sqrt(2)
gamma   = np.pi/2
p1=p2=p3= np.sqrt(5)

''' Suggested Activity 1
Params L1, L2, L3, γ, x1, x2, y2 are fixed constants,
strut lengths p1, p2, p3 will be known for a given pose.
'''
print("\n Suggested Activity 1 \n")
# Fallið f skilar f(theta, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2)
# fyrir gefin gildi á theta, p1, p2, p3, L1, L2, L3, gamma, x1, x2 og y2
def f(theta, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2):
    A2      = L3*np.cos(theta)-x1
    B2      = L3*np.sin(theta)
    A3      = L2*np.cos(theta+gamma)-x2
    B3      = L2*np.sin(theta+gamma)-y2
    aNum    = p2**2-p1**2-A2**2-B2**2
    bNum    = p3**2-p1**2-A3**2-B3**2
    D       = 2*(A2*B3-B2*A3)
    N1      = B3*aNum-B2*bNum
    N2      =-A3*aNum+A2*bNum 
    return N1**2+N2**2-p1**2*D**2

print("Látum theta = -pi/4 og theta = pi/4")
print("Þá ætti f(theta) = 0")
print("Nú fæst að f(theta) er:")
print(f(np.pi/4, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2))
print("og")
print(f(-np.pi/4, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2))     

'''
Suggested Activity 2
Plots f(θ) on [−π,π].
Her sest nokkud merkilegt nefnilega hvad grafid er nalaegt nulli lengi.
Hefur sannarlega nullpunkta i +/- 0.25*pi eins og lysingin segir.
Maetti teikna thad staerra eda skyrar...
'''
print("\n Suggested Activity 2 \n")
def fPlot(thetaL, thetaH, fileName, graphTitle):
        # Setjum bilið á theta sem við viljum
        # teikna f(theta) fyrir
        theRange = np.arange(thetaL, thetaH, 0.1)
        # Upphafsstillum myndina og ásana
        figure, axis = plt.subplots()
        # Teiknum grafið
        axis.plot(theRange, f(theRange, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2))
        # Setjum titil á grafið og ásana
        axis.set(xlabel=r"${\Theta}$ [radians]", 
                ylabel="$f({\Theta})$",
                title = graphTitle)
        axis.set_ylim([-50000,50000])
        axis.grid()
        # Vistum myndina
        fileString = fileName + ".png"
        figure.savefig(fileString)
        # Setjum minnsta sýnilega bilið á x - ásnum
        plt.xticks(np.arange(thetaL, thetaH, 0.5))
        plt.yticks(np.arange(-50000,50000,5000))
        plt.show()

# Teiknum f(theta) fyrir theta frá -pi upp í pi
# fPlot(-np.pi, np.pi, "sa2", 'Suggested Activity 2')

# Finnum núllstöðvar f(theta) með útreikningi
# Þær ættu að vera theta = pi/4 og theta = -pi/4

# fSolver skilar núllstöð f. Þarf upphafságiskun sem byggir á
# því að skoða graf af f(theta)
def fSolver(f, thetaGuess): 
        func = lambda theta : f(theta, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2)
        # Aðferðin sem fsolve notar til að finna núllstöð: 
        # https://www.math.utah.edu/software/minpack/minpack/hybrd.html
        thetaSol = scOpt.fsolve(func, thetaGuess)
        return thetaSol

# Skoðum graf af f(theta) og sjáum þá fyrir hvaða gildi
# á theta, f(theta) er nálægt núlli
thetaCalculated = fSolver(f, np.pi/5)
thetaGiven = np.pi/4
print("Núllstöð f(theta) sem fékkst með ítrun er:", thetaCalculated)
print("Gefin núllstöð f(theta) er:", thetaGiven)

thetaCalculated = fSolver(f, -np.pi/5)
thetaGiven = -np.pi/4
print("Núllstöð f(theta) sem fékkst með ítrun er:", thetaCalculated)
print("Gefin núllstöð f(theta) er:", thetaGiven)

'''
Suggested Activity 3
Reproduce Figure 1.15
'''
print("\n Suggested Activity 3 \n")
# Latum xL2 := x + L2*cos(theta + gamma)
#       yL2 := y + L2*sin(theta + gamma)
#       xL3 := x + L3*cos(theta + gamma)
#       yL3 := y + L3*sin(theta + gamma)

# Teiknar Stewart Platform
def plotStewartPlatform(x, y, x1, x2, y2, xL2, yL2, xL3, yL3, plotTitle="Vantar titil!", plotName="Vantar"):
        # Teiknar þrihyrning. Byrjum lengst til 
        # vinstri og forum rettsaelis
        # Her er (x, y) = (1,2), 
        #        (xL2, yL2) = (2,3)
        #        (xL3, yL3) = (2,1)
        (triangleX, triangleY)  = [x, xL2, xL3, x],[y, yL2, yL3, y]
        fig, plotObject = plt.subplots()
        plotObject.plot(triangleX, triangleY, 'b', linewidth=2.5)

        # Baetum vid fyrsta akkerinu
        # sem hefur hnit fra (0,0) i (x,y)
        firstAnchorX, firstAnchorY  = [0,x],[0,y]
        plotObject.plot(firstAnchorX, firstAnchorY, 'b')

        # Baetum vid odru akkerinu
        # sem hefur hnit fra (x2, y2)
        # til (xL2, yL2) 
        secondAnchorX, secondAnchorY  = [x2,xL2],[y2,yL2]
        plotObject.plot(secondAnchorX, secondAnchorY, 'b')

        # Baetum vid thridja akkerinu
        # sem hefur hnit fra (xL3, yL3)
        # til (x1, 0)
        thirdAnchorX, thirdAnchorY  = [x1,xL3],[0,yL3]
        plotObject.plot(thirdAnchorX, thirdAnchorY, 'b')

        # Baetum vid punktum fyrir 
        # akkaeri eitt, tvo og thrju 
        # sem hafa hnitin
        # (0,0), (x2, y2) og (x1, 0)
        anchorPointsX, anchorPointsY  = [0, x1, x2], [0, 0, y2]
        plotObject.plot(anchorPointsX, anchorPointsY, 'bo')

        # Baetum vid punktum fyrir
        # hvert horn thrihyrningsins
        plotObject.plot(triangleX, triangleY, 'bo')

        plotObject.set(xlabel="x", 
                ylabel="y",
                title = plotTitle)

        # Vistum myndina
        plotSaveFile = plotName + ".png"
        fig.savefig(plotSaveFile)
        plt.show()

''' Plot SA3(a) '''
x, y = 1, 2
xL2, yL2 = 2, 3
xL3, yL3 = 2, 1
plotName = "sa3a"
plotTitle = 'Suggested Activity 3 (a)'
# plotStewartPlatform(x, y, x1, x2, y2, xL2, yL2, xL3, yL3, plotTitle, plotName)

''' Plot SA3(b) '''
x, y = 2, 1
xL2, yL2 = 1, 2
xL3, yL3 = 3, 2
plotName = "sa3b"
plotTitle = 'Suggested Activity 3 (b)'
# plotStewartPlatform(x, y, x1, x2, y2, xL2, yL2, xL3, yL3, plotTitle, plotName)

'''
Suggested Activity 4
Solve the forward kinematics problem for the planar Stewart platform.
Plot f(theta) and then solve f(theta) = 0
'''
print("\n Suggested Activity 4 \n")
# Stikar
x1 = 5 
x2, y2 = 0, 6
L1 = L3 = 3 
L2 = 3*np.sqrt(2)
gamma = np.pi/4
p1 = p2 = 5
p3 = 3

# Teiknum f(theta) fyrir theta frá -pi upp í pi
fPlot(-np.pi, np.pi, "sa4-0", 'Suggested Activity 4')

# Finnum núna fjórar núllstöðvar f(theta)
thetas = np.zeros(4, dtype=float)

# Útfrá mynd sést að ágætis upphafsgisk fyrir
# þessar fjórar núllstöðvar f(theta) eru
# -0.7, -0.3, 1.1, 2.1
thetasGuess = [-0.7, -0.3, 1.1, 2.1]

# Hér fáum við út réttar núllstöðvar
print("Ef theta er raunveruleg núllstöð þá ætti")
print("f(theta) að vera mjög nálægt núlli \n")
for i in range(len(thetas)):
        thetas[i] = fSolver(f, thetasGuess[i])
        print("Nú fæst að theta er:", thetas[i])
        print("Svo er f(theta):", f(thetas[i], p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2), "\n")

# Geymir x og y gildin fyrir útreiknaðar núllstöðvar f(theta)
xy = np.zeros([4,2], dtype=float)

# Fallið xyCalc reiknar út gildin á x og y miðað við gefin
# gildi á theta, p1, p2, p3, L1, L2, L3, gamma, x1, x2 og y2.
def xyCalc(theta, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2):
    A2      = L3*np.cos(theta)-x1
    B2      = L3*np.sin(theta)
    A3      = L2*np.cos(theta+gamma)-x2
    B3      = L2*np.sin(theta+gamma)-y2
    aNum    = p2**2-p1**2-A2**2-B2**2
    bNum    = p3**2-p1**2-A3**2-B3**2
    D       = 2*(A2*B3-B2*A3)
    N1      = B3*aNum-B2*bNum
    N2      =-A3*aNum+A2*bNum 
    x = N1/D
    y = N2/D
    return [x, y]

for i in range(0, 4): 
        xy[i] = xyCalc(thetas[i], p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2)

# Útfrá útreiknuðum theta, x og y gildum getum við síðan reiknað xL2, yL2, xL3, yL3
xyL23 = np.zeros([4,4], dtype=float)

def xyLCalc(theta, x, y, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2): 
        xL2 = x + L2*np.cos(theta + gamma)
        yL2 = y + L2*np.sin(theta + gamma)
        xL3 = x + L3*np.cos(theta)
        yL3 = y + L3*np.sin(theta)
        return [xL2, yL2, xL3, yL3]

for i in range(0, 4):
        xyL23[i] = xyLCalc(thetas[i], xy[i,0], xy[i,1], p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2)

# Teiknum Stewart platform fyrir þessar fjórar núllstöðvar f(theta)
# og könnum í leiðinni hvort lengdirnar á struts séu p1, p2 og p3
for i in range(0, 4):
        p1Calc = math.hypot(xy[i,0] - 0, xy[i,1] - 0)
        p2Calc = math.hypot(xyL23[i,2] - x1, xyL23[i,3] - 0)
        p3Calc = math.hypot(xyL23[i,0] - x2, xyL23[i,1] - y2)
        print("Fyrir núllstöðina", thetas[i], "þá fást eftirfarandi niðurstöður")
        print("Nú er p1 = 5 en útreiknað gildi er:", p1Calc)
        print("Nú er p2 = 5 en útreiknað gildi er:", p2Calc)
        print("Nú er p3 = 3 en útreiknað gildi er:", p3Calc)
        print("\n")
        pName = "sa4-" + str(i+1)
        plotStewartPlatform(xy[i,0], xy[i,1], x1, x2, y2, xyL23[i, 0], xyL23[i, 1], xyL23[i, 2], xyL23[i, 3], "Suggested Activity 4", pName)

'''
Suggested Activity 5
Change strut length to p2 = 7 and re-solve the problem. 
For these parameters, there are six poses
'''

print("Suggested Activity 4 \n")
# Change the value of p2
p2 = 7

# Teiknum f(theta) fyrir theta frá -pi upp í pi
fPlot(-np.pi, np.pi, "sa5-0", 'Suggested Activity 5')

# Finnum núna sex núllstöðvar f(theta)
# Hér er ferillinn mjög flatur í kringum fimm
# af sex núllstöðvunum
thetas2 = np.zeros(6, dtype=float)

# Útfrá mynd sést að ágætis upphafsgisk fyrir
# þessar sex núllstöðvar f(theta) eru
# -0.7, -0.4, 0, 0.42, 1, 2.5
thetasGuess2 = [-0.7, -0.4, 0, 0.42, 1, 2.5]

# Hér fáum við út réttar núllstöðvar
print("Ef theta er raunveruleg núllstöð þá ætti")
print("f(theta) að vera mjög nálægt núlli \n")
for i in range(0,len(thetas2)):
        thetas2[i] = fSolver(f, thetasGuess2[i])
        print("Nú fæst að theta er:", thetas2[i])
        print("Svo er f(theta):", f(thetas2[i], p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2), "\n")

# Geymir x og y gildin fyrir útreiknaðar núllstöðvar f(theta)
xy2 = np.zeros([6,2], dtype=float)

for i in range(0, 6): 
        xy2[i] = xyCalc(thetas2[i], p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2)

# Útfrá útreiknuðum theta, x og y gildum getum við síðan reiknað xL2, yL2, xL3, yL3
xy2L23 = np.zeros([6,4], dtype=float)

for i in range(0, 6):
        xy2L23[i] = xyLCalc(thetas2[i], xy2[i,0], xy2[i,1], p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2)

# Teiknum Stewart platform fyrir þessar sex núllstöðvar f(theta)
# og könnum í leiðinni hvort lengdirnar á struts séu p1, p2 og p3
for i in range(0, 6):
        p1Calc2 = math.hypot(xy2[i,0] - 0, xy2[i,1] - 0)
        p2Calc2 = math.hypot(xy2L23[i,2] - x1, xy2L23[i,3] - 0)
        p3Calc2 = math.hypot(xy2L23[i,0] - x2, xy2L23[i,1] - y2)
        print("Fyrir núllstöðina", thetas2[i], "þá fást eftirfarandi niðurstöður")
        print("Nú er p1 = 5 en útreiknað gildi er:", p1Calc2)
        print("Nú er p2 = 7 en útreiknað gildi er:", p2Calc2)
        print("Nú er p3 = 3 en útreiknað gildi er:", p3Calc2)
        print("\n")
        pName = "sa5-"+str(i+1)
        plotStewartPlatform(xy2[i,0], xy2[i,1], x1, x2, y2, xy2L23[i, 0], xy2L23[i, 1], xy2L23[i, 2], xy2L23[i, 3], "Suggested Activity 5", pName)

'''
Suggested Activity 6
Find a strut length p2, with the rest of the parameters as in Step 4, 
for which there are only two poses.
'''
# Stikar úr Suggested Activity 4
x1 = 5 
x2, y2 = 0, 6
L1 = L3 = 3 
L2 = 3*np.sqrt(2)
gamma = np.pi/4
p1 = 5
p3 = 3

# Fall sem telur fjölda róta fallsins f(theta)
# fyrir gefna stika
def countRoots(f, p2): 
        deltaTheta = 0.001
        thetaOld = -np.pi
        oldSign = f(thetaOld, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2) > 0
        count = 0
        for theta in np.arange(-np.pi, np.pi+deltaTheta, deltaTheta):
                thetaNew = theta
                newSign = f(thetaNew, p1, p2, p3, L1, L2, L3, gamma, x1, x2, y2) > 0
                if(newSign != oldSign):
                        count += 1
                oldSign = newSign
        return count

# Höfum séð að það eru 4 rætur og þar með 4 stöður þegar p2 = 5
# Höfum séð að það eru 6 rætur og þar með 6 stöður þegar p2 = 7
# Skoðum gildi þar í kring
# Vitum að það eru (líklegast) formerkjaskipti í rótum
#   svo leitum að formerkjaskiptum á bilinu [0,12]

# Prófum mismunandi gildi á p2
p2Array = np.arange(0, 12.5, 0.5)
p2Roots = list()
for i in range(len(p2Array)):
        if(countRoots(f, i) == 2):
                p2Roots.append(i)
print(p2Roots)

'''
Suggested Activity 7
Calculate the intervals in p2, with the rest of the parameters
as in Step 4, for which there are 0, 2, 4 and 6 poses,
respectively 
'''
# Hér voru skoðuð gildi á p frá 0 upp í 30

# Stikar úr Suggested Activity 4
x1 = 5 
x2, y2 = 0, 6
L1 = L3 = 3 
L2 = 3*np.sqrt(2)
gamma = np.pi/4
p1 = 5
p3 = 3

no = 100
noSolutions = []
print(noSolutions)
i = 0.0
for p2 in range(no):
    noRoots = countRoots(f, i)
    noSolutions.append(countRoots(f, i))
    i+=0.2

noSolutions = np.array(noSolutions)

theRange = np.arange(0, 100, 1)

figure, axis = plt.subplots()

# Teiknum grafið
axis.plot(theRange, noSolutions[theRange])

# Setjum titil á grafið og ásana
axis.set(xlabel=r"${p_2}$", 
                ylabel="Fjöldi staða (núllstöðva)",
                title = "")
axis.set_ylim([0,6])
axis.grid()

# Setjum minnsta sýnilega bilið á x - ásnum
plt.xticks(np.arange(0, 100, 5))
plt.yticks(np.arange(0, 6, 1))
# Vistum myndina
fileString = "sa7" + ".png"
figure.savefig(fileString)
plt.show()