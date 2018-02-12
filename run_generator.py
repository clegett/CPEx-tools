import mmm
import math

RHOST = 5000
DRRIM = 100
RINCL = 10
NINCL = 500
HOSTOC = [1.64006733801622,1.00010324864145E-05]
RIMOC = [2.117961207,1.1E-07]
INCLOC = [3.42823,3.39901]
wavelength = 700
mypack = mmm.Pack.from_file('example_input/100.dat')

mypack.rescale_pack(RHOST + DRRIM)
mypack.center_pack()
mycluster = mmm.Cluster()
mycluster.grainlist = []

for sphere in mypack.sphere_coords:
    mycluster.grainlist.append(mmm.Grain.new_grain(sphere[0], sphere[1],
    sphere[2], RHOST, DRRIM, RINCL, NINCL))
for grain in mycluster.grainlist:
    grain.set_grain_oc(HOSTOC[0], HOSTOC[1], RIMOC[0], RIMOC[1], INCLOC[0],
                       INCLOC[1])

myrun = mmm.ModelOptions()
myrun.number_spheres = (1 + 1 + NINCL) * len(mycluster.grainlist)
myrun.length_scale_factor = (2*math.pi/wavelength)
myrun.output_file = 'output/' + str(wavelength) + 'nm.dat'
myrun.scattering_coefficient_file = 'output/' + str(wavelength) + 'nm.sc.dat'
runoutput = myrun.get_formatted('required')
print(runoutput)
print('sphere_sizes_and_positions')

for grain in mycluster.grainlist:
    grain.print_rxyznk()
