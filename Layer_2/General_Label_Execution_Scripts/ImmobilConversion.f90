! TITLE: ImmobilConversion.f90
! AUTHOR: Rivaan Kakkaramadam
! DESCRIPTION: This script takes in the dye label trajectory file outputted by one of the LAMMPS sim files here 
! to perform a Euclidean distance-based colocalization check and appropriate tagging of all dye labels. 
! Output is a new dye label input file where particle type corresponds to the dye label tag (2 = uncolocalized, 3 = colocalized).

program immobilconversion
  implicit none
  real*8 :: wx, wy, wz, xr, yr, zr, r, srevradius, protradius, x, y, z, kr, xlength, ylength, zlength
  integer :: linecount, startline, saves, i, j, k, n, w, npart, walkercount
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
        read(1,*) npart
     else
        read(1,*)
     endif
  enddo
  walkercount = 10000

  allocate(srevx(npart), srevy(npart), srevz(npart))
  
  do i = 1, (npart)
     read(1,*) k, srevtype, x, y, z,kr,kr,kr,kr,kr,kr,kr
        srevx(i) = x
        srevy(i) = y
        srevz(i) = z
  enddo

  CALL execute_command_line('rm fort.2')
  CALL execute_command_line('ln -s walkers.xyz fort.2') ! read in trajectory file of dye labels


  allocate(walkx(walkercount), walky(walkercount), walkz(walkercount))
  allocate(bindarray(walkercount)) ! tracks binding statuses of dye labels
  do j = 1, saves
     print *, "j", j
     read(2,*) walkercount
     read(2,*) tmp, tmp1, nt
     print *, "entering walker loop"
     if (nt .eq. extractionpoint) then
     do i = 1, walkercount
        read(2,*) k, wx, wy, wz
        if (k .eq. 3) then
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
        xr = xr - xlength*ANINT(xr/xlength) ! adjust computed distance for periodic boundary conditions in sim
        yr = yr - ylength*ANINT(yr/ylength)
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
        do i = 1, walkercount
           read(2,*) k, wx, wy, wz
        enddo
     endif
enddo

print *, 'finished iterating'

open(unit = 9, file = "post-immobil-input.in")
  do i = 1, walkercount
     write(9,*) "create_atoms", bindarray(i), "single", walkx(i), walky(i), walkz(i)
  enddo
close(9)

print *, 'finished writing'

 deallocate(srevx, srevy, srevz, walkx, walky, walkz, bindarray)
  
end program immobilconversion

  
