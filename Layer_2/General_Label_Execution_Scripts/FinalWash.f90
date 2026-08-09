! TITLE: FinalWash.f90
! AUTHOR: Rivaan Kakkaramadam
! DESCRIPTION: This script takes in the dye label trajectory file outputted by the last iteration of the diffusivity simulation
! to perform a Euclidean distance-based colocalization check and appropriate tagging of all dye labels. The "wash" step of this model
! chromatin sample is simulated by producing a list of dye labels that filters out the dye labels that had not bound to the model chromatin
! by the end of the simulation and thus would get "washed" out of the model chromatin sample.
! Output is the coordinates of only the dye labels that are tagged as colocalized by the end of the full simulation and thus would be 
! the labels that remain on the model chromatin sample post-wash.

program finalwash
  implicit none
  real*8 :: wx, wy, wz, xr, yr, zr, r, srevradius, protradius, x, y, z, kr, xlength, ylength, zlength
  integer :: linecount, saves, i, j, k, n, w, npart, walkercount
  integer :: nt, binding, extractionpoint,srevtype
  integer, allocatable :: bindarray(:)
  real*8, allocatable :: srevx(:), srevy(:), srevz(:), walkx(:), walky(:), walkz(:)
  character(len=100) :: syncommand, wccommand
  character(len=50) :: configfilename
  character(len=15) :: tmp, tmp1
  

  xlength=130 ! define box bounds - ENSURE THIS IS CONSISTENT WITH PARAMS IN SREV CONFIG FILE
  ylength=130 ! xlength, ylength, zlength should be length of simulation box in x,y,z directions
  zlength=130 
  srevradius = 1.00d0/2.0d0
  protradius = 1.00d0/2.0d0
  saves = 21
  extractionpoint = 20000
  configfilename = 'edited-config-8.dump'

  open(unit = 1, file = trim(configfilename))
  do i = 1, 9
     if (i .eq. 4) then
        read(1,*) npart ! read in the number of SREV nucleosomes
     else
        read(1,*)
     endif
  enddo
  walkercount = 10000
  allocate(srevx(npart), srevy(npart), srevz(npart)) ! allocate sufficient memory to arrays storing data on SREV nucleosomes

  do i = 1, (npart) ! read in SREV nucleosome coordinates
     read(1,*) k, srevtype, x, y, z,kr,kr,kr,kr,kr,kr,kr
     srevx(i) = x
     srevy(i) = y
     srevz(i) = z
  enddo

  CALL execute_command_line('rm fort.2')
  CALL execute_command_line('ln -s walkers.xyz fort.2')


  allocate(walkx(walkercount), walky(walkercount), walkz(walkercount))
  allocate(bindarray(walkercount))
  do j = 1, saves ! iterate over flurophore trajectory file
     print *, "j", j
     read(2,*) walkercount
     read(2,*) tmp, tmp1, nt
     print *, "entering walker loop"
     if (nt .eq. extractionpoint) then
     do i = 1, walkercount
        read(2,*) k, wx, wy, wz
        if (k .eq. 3) then ! if the dye label was already tagged as SREV-colocalized there is no need to repeat the distance check to keep it tagged as coloc
           bindarray(i) = k
           walkx(i) = wx
           walky(i) = wy
           walkz(i) = wz
        else
     binding = 2
     do n = 1, npart ! check against SR-EV beads
        xr = wx - srevx(n)
        yr = wy - srevy(n)
        zr = wz - srevz(n)
        xr = xr - xlength*ANINT(xr/xleyngth) ! adjust computed distance for periodic boundary conditions in sim
        yr = yr - length*ANINT(yr/ylength)
        zr = zr - zlength*ANINT(zr/zlength)
        r = sqrt(xr**2 + yr**2 + zr**2)
        if (r < (srevradius+protradius)) then
           binding = 3
        endif
     enddo
     walkx(i) = wx
     walky(i) = wy
     walkz(i) = wz
     bindarray(i) = binding
     endif
     enddo
     else
        do i = 1, walkercount ! this ensures that before program gets to extraction point the file is properly iterated over
           read(2,*) k, wx, wy, wz
        enddo
     endif
enddo

print *, 'finished iterating'

open(unit = 9, file = 'BrdU-postwash-config-8.xyz')
do i = 1, walkercount
	if (bindarray(i) .eq. 3) then
     	write(9,*) bindarray(i), walkx(i), walky(i), walkz(i)
	endif
enddo

print *, 'finished writing'

deallocate(srevx, srevy, srevz, walkx, walky, walkz, bindarray)
  
end program finalwash

  
