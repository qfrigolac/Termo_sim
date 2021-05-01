from vpython import *
import numpy as np
#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood

win = 500

Natoms = 400  # change this to have more or fewer atoms

# Typical values
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container
mass = 4E-3/6E23 # helium mass
Ratom = 0.05 # wildly exaggerated size of helium atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature
dt = 1E-5

mpress = 0
m_glob_press = 0


animation = canvas( width=win, height=win, align='left')
animation.range = L
animation.title = 'A "hard-sphere" gas'
s = """  Theoretical and averaged speed distributions (meters/sec).
  Initially all atoms have the same speed, but collisions
  change the speeds of the colliding atoms. One of the atoms is
  marked and leaves a trail so you can follow its path.
  Codi de Bruce Sherwood. Modificat per Joaquim Frigola i Pau Sànchez.
  
"""
animation.caption = s

d = L/2+Ratom
r = 0.005
boxbottom = curve(color=gray, radius=r)
boxbottom.append([vector(-d,-d,-d), vector(-d,-d,d), vector(d,-d,d), vector(d,-d,-d), vector(-d,-d,-d)])
boxtop = curve(color=gray, radius=r)
boxtop.append([vector(-d,d,-d), vector(-d,d,d), vector(d,d,d), vector(d,d,-d), vector(-d,d,-d)])
vert1 = curve(color=gray, radius=r)
vert2 = curve(color=gray, radius=r)
vert3 = curve(color=gray, radius=r)
vert4 = curve(color=gray, radius=r)
vert1.append([vector(-d,-d,-d), vector(-d,d,-d)])
vert2.append([vector(-d,-d,d), vector(-d,d,d)])
vert3.append([vector(d,-d,d), vector(d,d,d)])
vert4.append([vector(d,-d,-d), vector(d,d,-d)])

wt = wtext(text='Pressió mitjana: {:1.2f}\n'.format(mpress))
wt2 = wtext(text='Temperatura mitjana: {:1.2f}'.format(mpress))
wt3 = wtext(text='PV: {:1.2f}'.format(mpress))
wt4 = wtext(text='nkT: {:1.2f}'.format(mpress))

def setpress(s):
    wt.text = 'Pressió mitjana: {:e} Pa\n'.format(s)
def settemp(s):
    wt2.text = 'Temperatura mitjana: {:.2f} K\n'.format(s)
def setPV(s):
    wt3.text = 'pV= {:e}\t'.format(s)
def setnkT(s):
    wt4.text = 'nkT: {:e} '.format(s)

Atoms = []
p = []
apos = []

for i in range(Natoms):
    x = L*random()-L/2
    y = L*random()-L/2
    z = L*random()-L/2
    u1 = random()
    u2 = random()
    u3 = random()
    u4 = random()
    u5 = random()
    u6 = random()
    px = sqrt(mass*k*T)*sqrt(-2*log(u1))*cos(2*np.pi*u2)
    py = sqrt(mass*k*T)*sqrt(-2*log(u3))*cos(2*np.pi*u4)
    pz = sqrt(mass*k*T)*sqrt(-2*log(u5))*cos(2*np.pi*u6)
    if i == 0:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z))
    p.append(vector(px,py,pz))

deltav = 100 # binning for v histogram

def barx(v):
    return int(v/deltav) # index into bars array

nhisto = int(4500/deltav)
histo = []
for i in range(nhisto): histo.append(0.0)
for i in range(Natoms): histo[barx(p[i].mag/mass)]+=1


gg = graph( width=win, height=0.4*win, xmax=3000, align='left',
    xtitle='speed, m/s', ytitle='Number of atoms', ymax=Natoms*deltav/1000)
theory = gcurve( color=color.cyan )    

dv = 10
for v in range(0,3001+dv,dv):  # theoretical prediction
    theory.plot( v, (deltav/dv)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) *exp(-0.5*mass*(v**2)/(k*T))*(v**2)*dv )

accum = []
for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0])
vdist = gvbars(color=color.red, delta=deltav )

t=0
gg2 = graph( width=win, height=0.4*win, align='left',
                xtitle='temps (t)', ytitle='p (Pa)')
press_graf = gcurve( color=color.cyan )
mpress_graf = gcurve( color=color.yellow )
m_glob_press_graf = gcurve( color=color.orange)

gg3 = graph( width=win, height=0.4*win, align='left',
                xtitle='temps (t)', ytitle='T (K)')
temp_graf = gcurve( color=color.cyan )
mtemp_graf = gcurve( color=color.orange )
tempx_graf = gcurve( color=color.yellow )

def interchange(v1, v2):  # remove from v1 bar, add to v2 bar
    barx1 = barx(v1)
    barx2 = barx(v2)
    if barx1 == barx2:  return
    if barx1 >= len(histo) or barx2 >= len(histo): return
    histo[barx1] -= 1
    histo[barx2] += 1

def checkCollisions():
    hitlist = []
    r2 = 2*Ratom
    r2 *= r2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i) :
            aj = apos[j]
            dr = ai - aj
            if mag2(dr) < r2: hitlist.append([i,j])
    return hitlist

nhisto = 0 # number of histogram snapshots to average


while True:
    rate(300)
    # Accumulate and average histogram snapshots
    for i in range(len(accum)): accum[i][1] = (nhisto*accum[i][1] + histo[i])/(nhisto+1)
    if nhisto % 10 == 0:
        vdist.data = accum
    nhisto += 1
    t+=dt
    # Update all positions
    for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt
    
    # Check for collisions
    hitlist = checkCollisions()

    # If any collisions took place, update momenta of the two atoms
    for ij in hitlist:
        i = ij[0]
        j = ij[1]
        ptot = p[i]+p[j]
        posi = apos[i]
        posj = apos[j]
        vi = p[i]/mass
        vj = p[j]/mass
        vrel = vj-vi
        a = vrel.mag2
        if a == 0: continue;  # exactly same velocities
        rrel = posi-posj
        if rrel.mag > Ratom: continue # one atom went all the way through another
    
        # theta is the angle between vrel and rrel:
        dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta)
        dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
        # alpha is the angle of the triangle composed of rrel, path of atom j, and a line
        #   from the center of atom i to the center of atom j where atome j hits atom i:
        alpha = asin(dy/(2*Ratom)) 
        d = (2*Ratom)*cos(alpha)-dx # distance traveled into the atom from first contact
        deltat = d/vrel.mag         # time spent moving from first contact to position inside atom
        
        posi = posi-vi*deltat # back up to contact configuration
        posj = posj-vj*deltat
        mtot = 2*mass
        pcmi = p[i]-ptot*mass/mtot # transform momenta to cm frame
        pcmj = p[j]-ptot*mass/mtot
        rrel = norm(rrel)
        pcmi = pcmi-2*pcmi.dot(rrel)*rrel # bounce in cm frame
        pcmj = pcmj-2*pcmj.dot(rrel)*rrel
        p[i] = pcmi+ptot*mass/mtot # transform momenta back to lab frame
        p[j] = pcmj+ptot*mass/mtot
        apos[i] = posi+(p[i]/mass)*deltat # move forward deltat in time
        apos[j] = posj+(p[j]/mass)*deltat
        interchange(vi.mag, p[i].mag/mass)
        interchange(vj.mag, p[j].mag/mass)
    
    mom_tan=0
    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L/2:
            mom_tan+=abs(p[i].x) #Per a mesurar la pressió
            if loc.x < 0: p[i].x =  abs(p[i].x)
            else: p[i].x =  -abs(p[i].x)
        
        if abs(loc.y) > L/2:
            mom_tan+=abs(p[i].y) #Per a mesurar la pressió
            if loc.y < 0: p[i].y = abs(p[i].y)
            else: p[i].y =  -abs(p[i].y)
        
        if abs(loc.z) > L/2:
            mom_tan+=abs(p[i].z) #Per a mesurar la pressió
            if loc.z < 0: p[i].z =  abs(p[i].z)
            else: p[i].z =  -abs(p[i].z)
            
            
    # Càlcul de la pressió
    press=2/(6*L*L*dt)*mom_tan #Noteu que al treballar amb moments ens estalviem posar la massa a tot arreu :)
    
    if nhisto % 10 == 0: 
        mpress_graf.plot((t,mpress/10))
        setpress(m_glob_press)
        setPV(m_glob_press)
        mpress=0
        m_glob_press_graf.plot((t,m_glob_press))
    
    mpress+=press
    press_graf.plot((t,press))
    m_glob_press = (nhisto*m_glob_press+press)/(nhisto+1)
    
    
    #Mesura de la temperatura
    suma=0
    suma2=0
    for i in range(Natoms): 
        suma2+=((p[i].mag)**2)/3
        suma+=(p[i].x)**2
    Temp=1/k*suma2/(Natoms*mass)
    Tempx=1/k*suma/(Natoms*mass)
    temp_graf.plot((t,Temp))
    tempx_graf.plot((t,Tempx))
    setnkT(Natoms*k*Temp)
    settemp(Temp)
    
    #if nhisto % 10 == 0: 
    #   print('P =',press, 'Pa T=',Temp,'K') 
