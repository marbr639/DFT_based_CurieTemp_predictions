import sys
import numpy as np

if len(sys.argv) >= 0:
    input_lst = [float(i) for i in sys.argv[1:]]
    length=len(input_lst)

    nAtoms = input_lst[:length]

c = {'config MAGMOM': []}
for n in nAtoms:
        for i in range(0,int(n)):
            point=np.array((0.0,0.0,0.0));


            x=np.random.rand(4)*2-1;
            x2=np.square(x)

            while (np.sum(x2)>=1.0):
                x=np.random.rand(4)*2-1;
                x2=np.square(x);
            mag = 3
            point[0]=2*(x[1]*x[3]+x[0]*x[2])/np.sum(x2);
            point[1]=2*(x[2]*x[3]-x[0]*x[1])/np.sum(x2);
            point[2]=(x2[0]+x2[3]-x2[1]-x2[2])/np.sum(x2);
            c['config MAGMOM'].append([point[0]*mag, point[1]*mag, point[2]*mag])

MAGMOM_str = ''
MAGMOM = c['config MAGMOM']

for m in MAGMOM:
    MAGMOM_str = MAGMOM_str + '{} {} {} '.format(m[0], m[1], m[2])
    
print(f'MAGMOM = {MAGMOM_str}')
