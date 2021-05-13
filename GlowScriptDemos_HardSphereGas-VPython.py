#Termodinàmica i Mecànica Estadística - Treball de Simulació
#Pau Sánchez Pascual i Joaquim Frigola Casals


from vpython import *
import numpy as np

#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood. Modificat per Pau Sànchez i Joaquim Frigola

win = 500


###### VARIABLES GENERALS I INICIALS ######
Natoms = 400  # change this to have more or fewer atoms


L = 1 # container is a cube L on a side

gray = color.gray(0.7) # color of edges of container
mass = 4E-3/6E23 # helium mass
Ratom = 0.05 # wildly exaggerated size of helium atom
k = 1.4E-23 # Boltzmann constant
T0 = 300 # around room temperature
dt = 1E-5

nderiv=50 #Iter que farà en cada cas (per assolir l'eq) al càlcul de les deriv.

deltaL=0.1
dq4 = 1
#############################################


#Inicialitzem el canvas:
animation = canvas( width=win, height=win, align='left') #Creació del canvas
animation.range = L
animation.title = 'Gas de "boles dures"'
s = """  Distribució teòrica i pràctica de les velocitats (m/sec).
  Al inici apliquem la distribució de Maxwell-Boltzmann perquè
  les partícules es trobin en equilibri tèrmic.
  Codi de Bruce Sherwood. Modificat per Joaquim Frigola i Pau Sànchez.
"""
animation.caption = s

###### TEXT PER A IMPRIMIR INFORMACIÓ ######
val=-1 #Valor per defecte fins que el text sigui actualizat per primer cop
wt = wtext(text=' Pressió mitjana: {:1.2f}\n'.format(val))
wt2 = wtext(text=' Temperatura mitjana: {:1.2f}\n'.format(val))
wt3 = wtext(text=' pV= {:1.2f}\t'.format(val))
wt4 = wtext(text=' NkT= {:1.2f}\n'.format(val))
wt5 = wtext(text=' n= {:.0f}'.format(val))
wt6 = wtext(text=' \n')
wt7 = wtext(text=' \n')


def setpress(s): #Funcions per a actualitzar el text
    wt.text = ' Pressió mitjana: {:e} Pa\n'.format(s)
def settemp(s):
    wt2.text = ' Temperatura mitjana: {:.2f} K\n'.format(s)
def setPV(s):
    wt3.text = ' pV= {:e}\t'.format(s)
def setnkT(s):
    wt4.text = ' NkT= {:e}\n '.format(s)
def setn(s):
    wt5.text = ' n= {:.0f}\t'.format(s)
def setflag(s):
    wt6.text = ' Iteració {:.0f} de 28\n'.format(s)
def seta(s):
    wt7.text = ' Coef. dilatació= {:.5f}\n'.format(s)



###### VARIABLES DE FUNCIONAMENT ######
Lx = L/2  #L'usarem per allargar la caixa més endavant


Atoms = []  #Llista d'àtoms
p = []  #Moment dels àtoms
apos = []  #Posició dels àtoms
histo = []   #Particions de l'histograma

T=T0

deltav = 100 # binning for v histogram
deltaa = 0.0005
nhisto = int(4500/deltav)
nhisto2 = int(2*45000*deltaa)

t=0  #Per controlar el temps que ha passat (gràfics)
n=0  #Nombre de passos (t=n*dt) que serà més útil per fer mitj que la var ant.

m_glob_press=0 #Mitjana de la pressió des de t=0 fins a t
mpress=0

#Derivades:
press_vol = []
temp_vol = []
press_temp = []
vol_temp = []
temp_temp = []

coef_dil=[] #Per si després volem fer estadístiques del coef de dilatació
p_exp=[] #Per si després volem fer estadístiques de p vs 1/V
inv_v_exp=[] #Per si després volem fer estadístiques de p vs 1/V


flag=0 #Per comptar les derivades

###### FUNCIONS DE INICI I RESET ######

def barx(v):
    return int(v/deltav) # index into bars array

def bara(a):
    return int(a/deltaa)


def inicialitzacio(Atoms,p,apos,histo,nhisto,Ratom):
    d = L/2+Ratom
    dx = Lx + Ratom
    r = 0.005
    boxbottom = curve(color=gray, radius=r)
    boxbottom.append([vector(-dx,-d,-d), vector(-dx,-d,d), 
                      vector(dx,-d,d), vector(dx,-d,-d), vector(-dx,-d,-d)])
    boxtop = curve(color=gray, radius=r)
    boxtop.append([vector(-dx,d,-d), vector(-dx,d,d),
                   vector(dx,d,d), vector(dx,d,-d), vector(-dx,d,-d)])
    vert1 = curve(color=gray, radius=r)
    vert2 = curve(color=gray, radius=r)
    vert3 = curve(color=gray, radius=r)
    vert4 = curve(color=gray, radius=r)
    vert1.append([vector(-dx,-d,-d), vector(-dx,d,-d)])
    vert2.append([vector(-dx,-d,d), vector(-dx,d,d)])
    vert3.append([vector(dx,-d,d), vector(dx,d,d)])
    vert4.append([vector(dx,-d,-d), vector(dx,d,-d)])
    for i in range(Natoms):
        #Posicions aleatòries
        x = L*random()-Lx #Durant l'expansió només a un extrem de la caixa
        y = L*random()-L/2
        z = L*random()-L/2
        #Mètode de Box-Muller per calcular el moment:
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
            Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, 
                                color=color.cyan))
        else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
        apos.append(vec(x,y,z)) #Guardem la posició
        p.append(vector(px,py,pz)) #Guardem el moment
    
    for i in range(nhisto): histo.append(0.0)
    for i in range(Natoms): histo[barx(p[i].mag/mass)]+=1

def resetvars(animation,Atoms,p,apos,histo,vdist,press_graf,
              mpress_graf,m_glob_press_graf,temp_graf,tempx_graf):#Netejar tot
    global n
    global t
    global n2
    
    Atoms.clear()
    p.clear()
    apos.clear()
    histo.clear()
    
    vdist.delete()
    press_graf.delete()
    mpress_graf.delete()
    m_glob_press_graf.delete()
    temp_graf.delete() 
    tempx_graf.delete()
    for obj in animation.objects:
        obj.visible=False
        del obj
    t=0
    n=0
    n2=0


###### FUNCIONS DEL BOTÓ #######
running = False
reset = False
def Run(b):
    global running
    global reset
    running = not running
    if running: 
        b.text = "Parar derivades"
        reset = True
    else: b.text = "Fer derivades programades"
    
button(text="Fer derivades programades",bind=Run)

restart = False
def Rest(b):
    global restart
    restart = True
    
button(text="Restart",bind=Rest)

q4contl = 0


quest4 = False
def Q4(b):
    global quest4
    quest4 = True
    
button(text="Expansió lliure (pas)",bind=Q4)

inst = False
def Q4_inst(b):
    global inst
    inst = True
    
button(text="Instantània pV",bind=Q4_inst)

reset_press = False
def Res_press(b):
    global reset_press
    reset_press = True
    
button(text="R. Pressió Mitjana",bind=Res_press)


###### INICI I CREACIÓ DE GRÀFICS ######
inicialitzacio(Atoms,p,apos,histo,nhisto,Ratom)



gg = graph( width=win, height=0.4*win, xmax=3000, align='left',
    xtitle='speed, m/s', ytitle='Number of atoms', ymax=Natoms*deltav/1000)

theory = gcurve( color=color.cyan )
dv = 10
for v in range(0,3001+dv,dv):  # theoretical prediction
    theory.plot( v, (deltav/dv)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) *
                exp(-0.5*mass*(v**2)/(k*T))*(v**2)*dv )

accum = []
for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0])
vdist = gvbars(color=color.red, delta=deltav )

gg2 = graph( width=win, height=0.4*win, align='left',
                xtitle='temps (s)', ytitle='p (Pa)')
press_graf = gcurve( color=color.cyan )
mpress_graf = gcurve( color=color.yellow )
m_glob_press_graf = gcurve( color=color.orange)

gg3 = graph( width=win, height=0.4*win, align='left',
                xtitle='temps (s)', ytitle='T (K)')
temp_graf = gcurve( color=color.cyan )
tempx_graf = gcurve( color=color.yellow )

gg4 = graph( width=win, height=0.4*win, align='left',
                xtitle='1/V (m^-3)', ytitle='p (Pa)')
pV_graf = gcurve( color=color.cyan )


histo2=[]
for i in range(nhisto2): histo2.append([deltaa*(i+0.5),0.0])
gg5= graph( width=win, height=0.4*win, align='left',
                xtitle='valor', ytitle='coef. dilatació')

vdist2 = gvbars(color=color.red, delta=deltaa )



###### FUNCIONS DE FUNCIONAMENT ######

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


###### PROGRAMA ######

n = 0 # number of iterations
n2 = 0 # per a fer la mitjana global de la p i poder-ho resetejar
while True:
    rate(300)
    t += dt
    
    # Accumulate and average histogram snapshots
    for i in range(len(accum)): accum[i][1] = (n*accum[i][1] + histo[i])/(n+1)
    if n % 10 == 0:
        vdist.data = accum


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
        if abs(loc.x) > Lx:
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
    press=2/(6*L*L*dt)*mom_tan 
    #Noteu que al treballar amb moments estalviem posar la massa a tot arreu 
    
    if n % 10 == 0: #Mitjana cada 10 passos 
        mpress_graf.plot((t,mpress/10)) #Actualiztem el gràfic
        setpress(m_glob_press) #Actualitz el text aquí pq no canvïi tant ràpid
        setPV(m_glob_press*L*L*Lx*2)
        mpress=0
    
    mpress+=press #Per fer-ne la mitjana
    m_glob_press = (n2*m_glob_press+press)/(n2+1) #Computem la mitjana global
    press_graf.plot((t,press)) #Gràfics
    m_glob_press_graf.plot((t,m_glob_press))
    
    #Mesura de la temperatura
    mtemp=0
    mtempx=0
    for i in range(Natoms): 
        mtemp+=((p[i].mag)**2)/3 #Suma de moments de la mitjana de les 3 comps
        mtempx+=(p[i].x)**2 #Suma de moments només component x
    Temp=1/k*mtemp/(Natoms*mass)
    Tempx=1/k*mtempx/(Natoms*mass)
    temp_graf.plot((t,Temp))
    tempx_graf.plot((t,Tempx))
    setnkT(Natoms*k*Temp) #Actualitzem els textos
    settemp(Temp)
    setn(n)
    
    
    n += 1
    n2 += 1
    
    #Activacions amb el botó
    if reset: #Netejar i deixar-ho com al principi (durant les derivades):
        reset = False 
        L=1
        Lx = L/2
        T=T0
        flag=0
        q4contl = 0
        press_vol = []
        temp_vol = []
        press_temp = []
        vol_temp = []
        temp_temp = []
        resetvars(animation,Atoms,p,apos,histo,vdist,press_graf,
              mpress_graf,m_glob_press_graf,temp_graf,tempx_graf)
        inicialitzacio(Atoms,p,apos,histo,nhisto,Ratom)
    if quest4: #Fer l'expansió:
        running = False 
        L=1
        q4contl += 1
        Lx = (L+q4contl*dq4)/2
        
        T=T0
        flag=0
        press_vol = []
        temp_vol = []
        press_temp = []
        vol_temp = []
        temp_temp = []
        resetvars(animation,Atoms,p,apos,histo,vdist,press_graf,
              mpress_graf,m_glob_press_graf,temp_graf,tempx_graf)
        inicialitzacio(Atoms,p,apos,histo,nhisto,Ratom)
        quest4 = False
    if inst: #Capturar les dades de p vs 1/V
        if q4contl==0:
            pV_graf.delete()
        pV_graf.plot((1/(L*L*2*Lx),m_glob_press))
        inv_v_exp.append(1/(L*L*2*Lx))
        p_exp.append(m_glob_press)
        inst = False
    if reset_press: #Reset
        m_glob_press=press
        n2=0
        reset_press = False
    if restart:  #Parar-ho tot i recomençar de nou (botó restart):
        running = False 
        restart = False 
        L=1
        Lx = L/2
        T=T0
        flag=0
        q4contl = 0
        press_vol = []
        temp_vol = []
        press_temp = []
        vol_temp = []
        temp_temp = []
        resetvars(animation,Atoms,p,apos,histo,vdist,press_graf,
              mpress_graf,m_glob_press_graf,temp_graf,tempx_graf)
        inicialitzacio(Atoms,p,apos,histo,nhisto,Ratom)
    
    #Càlcul de les derivades i el coeficient de dilatació
    if running:
        if n==nderiv:
            setflag(flag+1)
            if flag<8: #Dp/DT a V ctant
                press_vol.append(m_glob_press)
                temp_vol.append(Temp)
                resetvars(animation,Atoms,p,apos,histo,vdist,press_graf,
                          mpress_graf,m_glob_press_graf,temp_graf,tempx_graf)
                inicialitzacio(Atoms,p,apos,histo,nhisto,Ratom)
                flag+=1
            elif flag < 28: #Dp/DV a T ctant
                press_temp.append(m_glob_press)
                temp_temp.append(Temp)
                vol_temp.append(L**3)
                if (flag-7)%5==0:
                    L+=deltaL
                    Lx = L/2
                resetvars(animation,Atoms,p,apos,histo,vdist,press_graf,
                          mpress_graf,m_glob_press_graf,temp_graf,tempx_graf)
                inicialitzacio(Atoms,p,apos,histo,nhisto,Ratom)
                flag+=1
            else: 
                #Fem les regressions i mitjanes
                div1=np.polyfit(temp_vol,press_vol,1)[0] #Numerador del coef,
                volum=[]
                suma=0
                for i in range(0,5): suma+=vol_temp[i]
                volum.append(suma/5)
                suma=0
                for i in range(5,10): suma+=vol_temp[i]
                volum.append(suma/5)
                suma=0
                for i in range(10,15): suma+=vol_temp[i]
                volum.append(suma/5)
                suma=0
                for i in range(15,20): suma+=vol_temp[i]
                volum.append(suma/5)
                press2=[]
                suma=0
                for i in range(0,5): suma+=press_temp[i]
                press2.append(suma/5)
                suma=0
                for i in range(5,10): suma+=press_temp[i]
                press2.append(suma/5)
                suma=0
                for i in range(10,15): suma+=press_temp[i]
                press2.append(suma/5)
                suma=0
                for i in range(15,20): suma+=press_temp[i]
                press2.append(suma/5)
                div2=np.polyfit(volum,press2,1)[0] #Denominador del coef
                alp=-div1/div2 #Calculem el coeficient de dilatació
                seta(alp) #Imprimeix el coeficient a la pantalla
                coef_dil.append(alp)
                
                if alp>0 and alp < 0.02225: histo2[bara(alp)][1]+=1
                #Pq l'if anterior? Si nderiv és massa petita llavors dona
                #valors molt distànts que no caben a l'histo (i si intentes 
                #plotejar-los el programa peta). Per T altes cal fer el mateix
                #amb l'histograma de velocitats (o desactivar-lo)
                
                vdist2.data = histo2
                reset = True
        
    
