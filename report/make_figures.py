import numpy as np
import pylab

# simulation parameters
nsteps = 100
N = 4000

# read in thermo information from lammps log
f = open('../data_files/log.lammps','r')

T_lmp = np.zeros(nsteps)
E_lmp = np.zeros(nsteps)
PE_lmp = np.zeros(nsteps)
time_lmp = np.zeros(nsteps)

for i in range(33):
    f.readline()

for i in range(nsteps):
    line = f.readline()
    time_lmp[i] = float(line.split()[0])
    T_lmp[i] = float(line.split()[1])
    PE_lmp[i] = float(line.split()[2])
    E_lmp[i] = float(line.split()[4])

f.close()

# read in thermo information from cbemd output
f = open('../data_files/run_lmp_compare.sh.o736917', 'r')

T_cbemd = np.zeros(nsteps)
E_cbemd = np.zeros(nsteps)
PE_cbemd = np.zeros(nsteps)
time_cbemd = np.zeros(nsteps)

for i in range(nsteps):
    line = f.readline()
    time_cbemd[i] = float(line.split()[0])
    T_cbemd[i] = float(line.split()[3])
    PE_cbemd[i] = float(line.split()[2])
    E_cbemd[i] = float(line.split()[4])
    
f.close()

# read in g(r) infomation for lammps data (this file came from dump.lammps + vmd)
f = open('../data_files/gr_lmp.dat','r')
r_lmp = []
gr_lmp = []

for line in f.readlines():
    r_lmp.append(float(line.split()[0]))
    gr_lmp.append(float(line.split()[1]))

f.close()

# read in g(r) infomation for lammps data (this file came from trajectory + vmd)
f = open('../data_files/gr_cbemd.dat','r')
r_cbemd = []
gr_cbemd = []

for line in f.readlines():
    r_cbemd.append(float(line.split()[0]))
    gr_cbemd.append(float(line.split()[1]))

f.close()

# read in timing results
f = open('../data_files/timing_results.txt')
nprocs_tmp = []
natoms_tmp = []
rs_tmp = []
runtime_tmp = []

for line in f.readlines():
    nprocs_tmp.append(int(line.split()[0]))
    natoms_tmp.append(int(line.split()[1]))
    rs_tmp.append(float(line.split()[2]))
    runtime_tmp.append(float(line.split()[4]))


nprocs = sorted(list(set(nprocs_tmp)))
natoms = sorted(list(set(natoms_tmp)))
rs = sorted(list(set(rs_tmp)))
runtime = np.zeros((len(nprocs), len(natoms), len(rs)))


for (i, j, k, l) in zip(nprocs_tmp, natoms_tmp, rs_tmp, runtime_tmp):
    i1 = nprocs.index(i)
    j1 = natoms.index(j)
    k1 = rs.index(k)
    runtime[i1][j1][k1] = l

f.close()

# idx is the first index over which to start plottting, averaging, etc. (discard first few timesteps for equilibration purposes)
idx = 20
label1 = 'CBEMD'
label2 = 'LAMMPS'

# print results
print 'Average CBEMD temperature: %2.4f' % np.average(T_cbemd[idx:])
print 'Average LAMMPS temperature: %2.4f' % np.average(T_lmp[idx:])

print 'Average CBEMD potential energy/atom: %2.4f' % np.average(PE_cbemd[idx:]/N)
print 'Average LAMMPS potential energy/atom: %2.4f' % np.average(PE_lmp[idx:])

print 'Average CBEMD total energy/atom: %2.4f' % np.average(E_cbemd[idx:]/N)
print 'Average LAMMPS total energy/atom: %2.4f' % np.average(E_lmp[idx:])

# make plots
pylab.plot(time_cbemd[idx:], T_cbemd[idx:], label=label1)
pylab.plot(time_lmp[idx:], T_lmp[idx:], label=label2)
pylab.xlabel('timestep')
pylab.ylabel('temperature')
pylab.legend()
pylab.savefig('T_compare.eps')
pylab.clf()

pylab.plot(time_cbemd[idx:], PE_cbemd[idx:]/N, label=label1)
pylab.plot(time_lmp[idx:], PE_lmp[idx:], label=label2)
pylab.xlabel('timestep')
pylab.ylabel('potential energy/atom (reduced units)')
pylab.legend()
pylab.savefig('PE_compare.eps')
pylab.clf()

pylab.plot(time_cbemd[idx:], E_cbemd[idx:]/N, label=label1)
pylab.plot(time_lmp[idx:], E_lmp[idx:], label=label2)
pylab.xlabel('timestep')
pylab.ylabel('total energy/atom (reduced units)')
pylab.legend()
pylab.savefig('E_compare.eps')
pylab.clf()

pylab.plot(r_cbemd, gr_cbemd, label=label1)
pylab.plot(r_lmp, gr_lmp, label=label2)
pylab.xlabel('r (reduced units)')
pylab.ylabel('g(r)')
pylab.legend()
pylab.savefig('gr_compare.eps')
pylab.clf()

for j in range(len(rs)):
    for i in range(len(natoms)):
        pylab.plot(nprocs, runtime[:,i, j], label = str(natoms[i]))
    pylab.legend()
    pylab.xlabel('number of processors')
    pylab.ylabel('runtime (seconds)')
    pylab.savefig('scaling_rs_%2.2f.eps' % rs[j])
    pylab.clf()
