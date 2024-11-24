import numpy as np
import cmath
import sys
import subprocess

#set the variables
total = 10 # total minimization steps
k1 = 0.01

tip = [0, 166, 151.9596547988024]
# 131 to 152 [0,154, 151.9596547988024]
# 111 to 130 [0, 146, 151.9596547988024]
# 81 to 110 [0, 135, 151.9596547988024]
# 46 to 80 [0, 127, 151.9596547988024]
# 0 to 45 [0, 120.24307698685955, 151.9596547988024]

C11 = 215988.610377874 # in MPa
C22 = 215612.144829691 
C33 = 216578.15561501 

C12 = 163950.700824468
C13 = 163762.511799617
C23 = 163762.511799617

C44 = 25971.0447906336
C55 = 26116.6337858838  
C66 = 25991.3147829455 


inp = sys.argv[1]
j = int(inp)

A = np.array([[C11, C12, C13, 0, 0, 0],
              [C12, C22, C23, 0, 0, 0],
              [C13, C23, C33, 0, 0, 0],
              [0,   0,   0,   C44, 0, 0],
              [0,   0,   0,   0, C55, 0],
              [0,   0,   0,   0, 0, C66]])

#read the data
g=open(f'./dump/boundary_dump_{j-1}.lmp','r')
lineb = g.readlines()
g.close()
n_boundary = int(lineb[3].split()[0])


g = open(f'./dump/inner_dump_{j-1}.lmp','r')
linei = g.readlines()
g.close()
n_inner = int(linei[3].split()[0])

number = linei[5].split()
xhi = float(number[1])
yvalues = linei[6].split()
yl = float(yvalues[0])
yh = float(yvalues[1])
Y_length = yh-yl

zvalues = linei[7].split()
zl = float(zvalues[0])
zh = float(zvalues[1])
Z_length = zh-zl


#get polar coordinate function
def polar(x,y,z):
    b = tip[1]
    c = tip[2]
    r = np.sqrt((y-b)**2+(z-c)**2)
    theta = np.arctan2((z-c),(y-b))

    return x,r,theta


#anisotropic displacement field functions
def get_inverse(A): #this equation gives the inverse of a matrix
    return np.linalg.inv(A)

def solve_eqn(s11,s16,s12,s66,s26,s22): # this one solve the ^4 equation and gives 4 roots
    coeff = [s11, -2*s16, 2*s12+s66, -2*s26, s22]
    a1,a2,a3,a4 = np.roots(coeff)
    #print('a1,a2,a3,a4',a1,a2,a3,a4)
    return a2, a4

def get_pq(s11,s12,s66,s22,a1,a2):
    p1 = s11*a1**2 +s12
    p2 = s11*a2**2 +s12
    q1 = s12*a1 +s22/a1
    q2 = s12*a2 +s22/a2
    return p1, p2, q1, q2

def displace(a1,a2,p1,p2,q1,q2,K,r,th):
    a = 1/(a1-a2)
    d = K*np.sqrt(2*r/np.pi)

    b1 = a1*p2*np.sqrt(np.cos(th)+a2*np.sin(th))
    c1 = a2*p1*np.sqrt(np.cos(th)+a1*np.sin(th))
    x = a*(b1-c1)*d

    b2 = a1*q2*np.sqrt(np.cos(th)+a2*np.sin(th))
    c2 = a2*q1*np.sqrt(np.cos(th)+a1*np.sin(th))
    y = a*(b2-c2)*d

    return x.real*1e10, y.real*1e10   #1e10 to convert the unit



#get all the variables ready
B = get_inverse(A)

sp11 = B[0][0]-B[0][2]*B[2][0]/B[2][2]
sp12 = B[0][1]-B[0][2]*B[2][1]/B[2][2]
sp66 = B[5][5]-B[5][2]*B[2][5]/B[2][2]
sp22 = B[1][1]-B[1][2]*B[2][1]/B[2][2]
sp16 = B[0][5]-B[0][2]*B[2][5]/B[2][2]
sp26 = B[1][5]-B[1][2]*B[2][5]/B[2][2]

a1,a2 = solve_eqn(sp11,sp16,sp12,sp66,sp26,sp22)

p1,p2,q1,q2 = get_pq(sp11, sp12, sp66, sp22,a1,a2)



#loop over every atom in the boundary and inner
#those are for the boundary atoms
ib = []
tb = []
coorb = []
vb = []
for i in range(len(lineb)): 
    if i>=9:
        values = lineb[i].split()
        ib.append(int(values[0]))
        tb.append(values[1])
        coorb.append([float(values[2]), float(values[3]), float(values[4])])
        #vb.append([float(values[5]), float(values[6]), float(values[7])])
        #coorb.append([float(values[2])*xhi, yl+float(values[3])*Y_length, zl+float(values[4])*Z_length])
        #vb.append([float(values[5])*xhi, yl+float(values[6])*Y_length, zl+float(values[7])*Z_length])

new_coorb = []
new_yb = []
new_zb = []
for i in range(len(coorb)):
    x,r,theta = polar(coorb[i][0],coorb[i][1],coorb[i][2])
    d_y, d_z = displace(a1,a2,p1,p2,q1,q2,k1,r*1e-10,theta)
    new_coorb.append([x, coorb[i][1]+d_y, coorb[i][2]+d_z])
    new_yb.append(coorb[i][1]+d_y)
    new_zb.append(coorb[i][2]+d_z)

ylo = min(new_yb)-4
yhi = max(new_yb)+4
zlo = min(new_zb)-4
zhi = max(new_zb)+4


#those are for the inner atoms
ii = []
ti = []
coori = []
vi = []
for i in range(len(linei)):
    if i>=9:
        values = linei[i].split()
        ii.append(int(values[0]))
        ti.append(values[1])
        coori.append([float(values[2]), float(values[3]), float(values[4])])
        #vi.append([float(values[5]), float(values[6]), float(values[7])])
        #coori.append([float(values[2])*xhi, yl+float(values[3])*Y_length, zl+float(values[4])*Z_length])
        #vi.append([float(values[5])*xhi, yl+float(values[6])*Y_length, zl+float(values[7])*Z_length])

new_coori = []
for i in range(len(coori)):
    x,r,theta = polar(coori[i][0],coori[i][1],coori[i][2])
    d_y, d_z = displace(a1,a2,p1,p2,q1,q2,k1,r*1e-10,theta) 
    #new_coori.append([x, coori[i][1], coori[i][2]])
    new_coori.append([x, coori[i][1]+d_y, coori[i][2]+d_z])



#write the displaced atoms
#write displaced boundary atoms:
g = open(f'./boundary_displace/boundary_displace_{j}.lmp','w')
g.close()

g = open(f'./boundary_displace/boundary_displace_{j}.lmp','a')
count = 0
g.write('# Bcc NbTaTiHf oriented X=[10-1] Y=[101] Z=[010].  \n')
g.write('\n')
g.write(f'{n_boundary} atoms \n')
g.write('4 atom types \n')
g.write('\n')
g.write(f'0.000000000000     {xhi}  xlo xhi \n')
g.write(f'{ylo}    {yhi}   ylo yhi \n')
g.write(f'{zlo}    {zhi}   zlo zhi \n')
g.write('\n')
g.write('Masses \n\n')
g.write('1 178.49  # Hf \n')
g.write('2 92.91  # Nb \n')
g.write('3 180.9479  # Ta \n')
g.write('4 47.867  # Ti \n\n')
g.write('Atoms #atomic \n\n')
for i in range(len(lineb)):
    if i>=9:
        #count+=1
        g.write(str(ib[i-9])+' '+tb[i-9]+' '+str(new_coorb[i-9][0])+' '+str(new_coorb[i-9][1])+' '+str(new_coorb[i-9][2])+'\n')

#the following lines are for the dynamics

#g.write('\n')
#g.write('Velocities \n\n')

#for i in range(len(lineb)):
#    if i>=9:
#        xv = vb[i-9][0]
#        yv = vb[i-9][1]
#        zv = vb[i-9][2]
#        g.write(str(ib[i-9])+' '+str(xv)+' '+str(yv)+' '+str(zv)+' \n')

g.close()


#write displaced inner atoms:
g = open(f'./inner_displace/inner_displace_{j}.lmp','w')
g.close()

g = open(f'./inner_displace/inner_displace_{j}.lmp','a')
count = 0
g.write('# Bcc NbTaTiHf oriented X=[10-1] Y=[101] Z=[010].  \n')
g.write('\n')
g.write(f'{n_inner} atoms \n')
g.write('4 atom types \n')
g.write('\n')
g.write(f'0.000000000000      {xhi}  xlo xhi \n')
g.write(f'{ylo}    {yhi}   ylo yhi \n')
g.write(f'{zlo}    {zhi}   zlo zhi \n')
g.write('\n')
g.write('Masses \n\n')
g.write('1 178.49  # Hf \n')
g.write('2 92.91  # Nb \n')
g.write('3 180.9479  # Ta \n')
g.write('4 47.867  # Ti \n\n')
g.write('Atoms #atomic \n\n')
for i in range(len(linei)):
    if i>=9:
        count+=1
        g.write(str(ii[i-9])+' '+ti[i-9]+' '+str(new_coori[i-9][0])+' '+str(new_coori[i-9][1])+' '+str(new_coori[i-9][2])+'\n')
   
# the following lines are for the dynamics

#g.write('\n')
#g.write('Velocities \n\n')

#for i in range(len(linei)):
#    if i>=9:
#        xv = vi[i-9][0]
#        yv = vi[i-9][1]
#        zv = vi[i-9][2]
#        g.write(str(ii[i-9])+' '+str(xv)+' '+str(yv)+' '+str(zv)+' \n')

g.close()
print('finish')
'''
g = open(f'./lr3_{j}.sh','w')
g.close()
g = open(f'./lr3_{j}.sh','a')
g.write('#!/bin/bash \n')
g.write('#SBATCH --partition=lr4 \n')
g.write('#SBATCH --account=pc_damage \n')
g.write('#SBATCH --qos=lr_normal \n')
g.write('#SBATCH --job-name fracture \n')
g.write('#SBATCH --nodes 4 \n')
g.write('#SBATCH --time 00:20:00 \n\n')
g.write('module load gcc/11.3.0 \n')
g.write('module load openmpi/4.1.1-gcc \n')
g.write('module load python \n')
g.write(f'mpirun  ~/compile/lammps-ace-ptm/src/lmp_mpi -in relax.in -v current {j} -v next {j+1}\n')
g.close()

if j<total:
    print('!!!!!')
    argument_to_pass = f"{j}"

# Call the shell script to submit the specified script to the supercluster
    subprocess.run(["./loop_submit.sh", argument_to_pass])
'''
