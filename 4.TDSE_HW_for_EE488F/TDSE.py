import numpy as np
from modl import operator, Input 


#Inputs & Unit conversion

atomic_unit_of_time_to_fs = 2.418884326509E-2
angsL = Input.Lx * 1.88973
l = 1
nstep = Input.nstep
dt = 0.001/atomic_unit_of_time_to_fs


# Make wave fucntion (also construct it's operator)
wave = operator.wave(angsL, Input.ngx, l, dt)
# Save wave function as txt
f = open("wave.txt",'w')
f.write("# t(fs) " + operator.grid(angsL, Input.ngx).plot_grid())
f.write("  0.0000 " + wave.plot_grid())
reflec = np.float64(0)
trans = np.float64(0)
poten = operator.Potential(angsL, Input.ngx)


inten=np.zeros(Input.ngx)
saveprob = []
# Time propagation
for i in range(nstep):
    for j in range(0, 10):
        if i < -1:
            wave.back_step()
        else:
            wave.next_step()
    reflec = 0
    trans = 0 
    #counting
    prob = np.abs(wave.grd)**2
    for k in range(0, Input.ngx):
        if k <= poten.left :
            reflec += prob[k]
        if k > poten.right :
            trans += prob[k]
    t = (i+1)*10 * dt * atomic_unit_of_time_to_fs
    f.write("%.6f" % t + wave.plot_grid())
    saveprob.append('%.6f    %.6f     %.6f\n' %(t,reflec,trans))
f.close()
# Print Potential Shape
f2 = open("Potential.txt",'w')
f2.write(operator.Potential(angsL, Input.ngx).plot_grid())
f2.close()

f3 = open("Probablity.txt",'w')
f3.write("# t(fs)    Reflection  Transmission\n")
f3.writelines(saveprob)
f3.close()

