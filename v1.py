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
import matplotlib.pyplot as plt


''' Vidvaerar breytur fyrir Suggested Activities 1-3 allavega '''
x1=y2   = 4     # tekid af fig 1.15
x2      = 0     # tekid af fig 1.15
L1      = 2     # skritid ad thessi se ekki notadur, en hann kemur ekki fram i jofnunni i bokinni
L2 = L3 = np.sqrt(2)
gamma   = np.pi/2
p1=p2=p3= np.sqrt(5)
error = 1e-10   # spyr að þessu í dæmatíma á eftir


''' Suggested Activity 1
Params L1,L2,L3,γ,x1,x2,y2 are fixed constants,
strut lengths p1,p2,p3 will be known for a given pose.
'''
def f(theta):
    ### Derived:
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
    
theta0 = np.pi/4
#print(f(theta0))     # -4.547473508864641e-13 simeq 0... Eigum vid ad meta ovissu?
#print(f(-theta0))
#print(f(0))
#print(f(np.pi))


'''
Suggested Activity 2
Plots f(θ) on [−π,π].
Her sest nokkud merkilegt nefnilega hvad grafid er nalaegt nulli lengi.
Hefur sannarlega nullpunkta i +/- 0.25*pi eins og lysingin segir.
Maetti teikna thad staerra eda skyrar...
'''
t2 = np.arange(-np.pi, np.pi, 0.1)

#fig2, ax2 = plt.subplots()
#ax2.plot(t2, f(t2))

#ax2.set(xlabel=r"${\Theta}$ [deg]", 
#        ylabel="$f({\Theta})$",
#        title ='Suggested Activity 2')
#ax2.grid()

#fig2.savefig("sa2.png")
#plt.show()



'''
Suggested Activity 3
Reproduce Figure 1.15
'''
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
        # plotObject.ylabel.set_rotation('0')     Hef reynt ad snua thessu y
        # plotObject.grid()

        plotSaveFile = plotName + ".png"
        fig.savefig(plotSaveFile)
        plt.show()

''' Plot (a) '''
x, y = 1, 2
xL2, yL2 = 2, 3
xL3, yL3 = 2, 1
plotName = "sa3a"
plotTitle = 'Suggested Activity 3 (a)'
plotStewartPlatform(x, y, x1, x2, y2, xL2, yL2, xL3, yL3, plotTitle, plotName)

''' Plot (b) '''
x, y = 2, 1
xL2, yL2 = 1, 2
xL3, yL3 = 3, 2
plotName = "sa3b"
plotTitle = 'Suggested Activity 3 (b)'
plotStewartPlatform(x, y, x1, x2, y2, xL2, yL2, xL3, yL3, plotTitle, plotName)