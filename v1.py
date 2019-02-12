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
print(f(theta0))     # -4.547473508864641e-13 simeq 0... Eigum vid ad meta ovissu?
print(f(0))
print(f(np.pi))


'''
Suggested Activity 2
Plots f(θ) on [−π,π].
Her sest nokkud merkilegt nefnilega hvad grafid er nalaegt nulli lengi.
Hefur sannarlega nullpunkta i +/- 0.25*pi eins og lysingin segir.
Maetti teikna thad staerra eda skyrar...
'''
t2 = np.arange(-np.pi, np.pi, 0.1)

fig2, ax2 = plt.subplots()
ax2.plot(t2, f(t2))

ax2.set(xlabel=r"${\Theta}$ [deg]", 
        ylabel="$f({\Theta})$",
        title ='Suggested Activity 2')
ax2.grid()

fig2.savefig("sa2.png")
plt.show()



'''
Suggested Activity 3
Reproduce Figure 1.15
'''

''' Plot (a) '''
t3a1, t3a2  = [1,2,2,1],[2,3,1,2]
t3a3, t3a4  = [0, x1, x2], [0, 0, y2]
t3a5, t3a6  = [0,1],[0,2]
t3a7, t3a8  = [0,2],[y2,3]
t3a9, t3a0  = [x1,2],[0,1]
        
fig3a, ax3a = plt.subplots()
ax3a.plot(t3a1, t3a2, 'b', linewidth=2.5)
ax3a.plot(t3a5, t3a6, 'b')
ax3a.plot(t3a7, t3a8, 'b')
ax3a.plot(t3a9, t3a0, 'b')
ax3a.plot(t3a1, t3a2, 'bo')
ax3a.plot(t3a3, t3a4, 'bo')

ax3a.set(xlabel="x", 
        ylabel="y",
        title ='Suggested Activity 3 (a)')
# ax3a.ylabel.set_rotation('0')     Hef reynt ad snua thessu y
# ax3a.grid()

fig3a.savefig("sa3a.png")
plt.show()


''' Plot (b) '''
t3b1, t3b2  = [1,3,2,1],[2,2,1,2]
t3b3, t3b4  = [0, x1, x2], [0, 0, y2]
t3b5, t3b6  = [0,2],[0,1]
t3b7, t3b8  = [0,1],[y2,2]
t3b9, t3b0  = [x1,3],[0,2]


fig3b, ax3b = plt.subplots()
ax3b.plot(t3b1, t3b2, 'b', linewidth=2.5)
ax3b.plot(t3b5, t3b6, 'b')
ax3b.plot(t3b7, t3b8, 'b')
ax3b.plot(t3b9, t3b0, 'b')
ax3b.plot(t3b1, t3b2, 'bo')
ax3b.plot(t3b3, t3b4, 'bo')


ax3b.plot(t3b1, t3b2, 'bo')
ax3b.plot(t3b3, t3b4, 'bo')

ax3b.set(xlabel="x", 
        ylabel="y",
        title ='Suggested Activity 3 (b)')
# ax3b.ylabel.set_rotation('0')     Hef reynt ad snua thessu y
# ax3b.grid()

fig3b.savefig("sa3b.png")
plt.show()
