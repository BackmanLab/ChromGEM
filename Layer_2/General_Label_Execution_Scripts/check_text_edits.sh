#!/bin/bash

size=15

# automatically trade out config-specific filenames in scripts and exchange label size for user-selected size

sizeinru=$(echo "scale=2; $size/10" | bc)

lammpsdump="config-relaxed-${SLURM_ARRAY_TASK_ID}.dump"
editedlammpsdump="edited-config-${SLURM_ARRAY_TASK_ID}.dump"
outputwashfile="BrdU-postwash-config-${SLURM_ARRAY_TASK_ID}.xyz"

sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" getconfig_nooverlaps.in
sed -i "s:^dump 1 ellip custom 100 edited-config.*:dump 1 ellip custom 100 ${editedlammpsdump} id type x y z c_0[*]:" getconfig_nooverlaps.in
sed -i "s:^set type 2 shape .*:set type 2 shape ${sizeinru} ${sizeinru} ${sizeinru}:" getconfig_nooverlaps.in

sed -i "s:^  configfilename = .*:  configfilename = '${lammpsdump}':" InitWalkers.f90
sed -i "s:^  sigmawalkerwalker = .*:  sigmawalkerwalker = ${sizeinru}:" InitWalkers.f90
sed -i "s/\r//g" InitWalkers.f90
gfortran -Ofast InitWalkers.f90 -o InitWalkers
./InitWalkers

sed -i "s:^  configfilename = .*:  configfilename = '${editedlammpsdump}':" ImmobilConversion.f90
sed -i "s:^  protradius = .*:  protradius = ${sizeinru}/2.0d0:" ImmobilConversion.f90
sed -i "s/\r//g" ImmobilConversion.f90

sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" overlapfilter.in
sed -i "s:^set type 2 shape .*:set type 2 shape ${sizeinru} ${sizeinru} ${sizeinru}:" overlapfilter.in
sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" first_sim.in
sed -i "s:^set type 2 shape .*:set type 2 shape ${sizeinru} ${sizeinru} ${sizeinru}:" first_sim.in
sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" continue_sim.in
sed -i "s:^set type 2 shape .*:set type 2 shape ${sizeinru} ${sizeinru} ${sizeinru}:" continue_sim.in
sed -i "s:^set type 3 shape .*:set type 3 shape ${sizeinru} ${sizeinru} ${sizeinru}:" first_sim.in

sed -i "s:^  configfilename = .*:  configfilename = '${editedlammpsdump}':" FinalWash.f90
sed -i "s:^open(unit = 9, file =.*:open(unit = 9, file = '${outputwashfile}'):" FinalWash.f90
sed -i "s:^  protradius = .*:  protradius = ${sizeinru}/2.0d0:" FinalWash.f90
sed -i "s/\r//g" FinalWash.f90

stop

# get version of SREV configs that lack all overlapping nucleosomes

$mpirunexec -np 32 $lammps -in getconfig_nooverlaps.in

# get needed values for formula for temp to input that produces target mobile temp

$mpirunexec -np 32 $lammps -in overlapfilter.in

echo "first ellipsoid particle number"
var=$(sed -n '4,4p' dump.ellipsoid)
echo $var
a=10000
c=$var
d=2.50
e=$(echo "scale=10; $a / ( $a + $c ) * $d" | bc -l)
echo $e
sed -i "s:^fix 1 all nvt temp .*:fix 1 all nvt temp $e $e 1.0:" first_sim.in
rm dump.ellipsoid

# run first iteration of simulation

$mpirunexec -np 32 $lammps -in first_sim.in

for i in $(seq 1 10); do
	# below rm command clears nothing when i==1; kept in anyway since it performs rm when needed for rest of i's in loop range
	rm post-immobil-input.in
	# produce version of label input data that marks labels that should be immobilized at the conclusion of that iteration
	gfortran -Ofast ImmobilConversion.f90
	./a.out
	rm walkers.xyz
	rm dump.ellipsoid*
	# get needed values for formula for temp to input that produces target mobile temp
	sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" overlapfilter_continue.in
	$mpirunexec -np 32 $lammps -in overlapfilter_continue.in
	var=$(sed -n '4,4p' dump.immob)
	var2=$(sed -n '4,4p' dump.ellipsoid)
	echo "immob particle numbers" # print statements monitoring changes in number of mobile particles
	echo $var
	a=$(echo "10000 - ( $var - $var2 )" | bc -l)
	echo $var2
	echo $a
	c=$var
	d=2.50
	e=$(echo "scale=10; $a / ( $a + $c ) * $d" | bc -l)
	echo $e
	sed -i "s:^fix 1 all nvt temp .*:fix 1 all nvt temp $e $e 1.0:" continue_sim.in
	rm dump.immob
	rm dump.ellipsoid
	# run next iteration of simulation
	$mpirunexec -np 32 $lammps -in continue_sim.in
done

gfortran -Ofast FinalWash.f90

./a.out
