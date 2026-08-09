program initializelabels
  implicit none
  integer :: i
  real*8 :: xr, yr, zr, sigmawalkerwalker, sigmachromatinchromatin, sigmawalkerchromatin, totalfrac
  real*8, allocatable :: allx(:)
  real*8, allocatable :: ally(:)
  real*8, allocatable :: allz(:)
  real*8 :: boxlength,xlength,ylength,zlength,volume, xlowerbound, ylowerbound, zlowerbound
  integer :: npart
  integer :: j
  real*8 :: r
  integer :: startline
  integer :: k
  integer :: linecount
  integer :: walkers
  integer, allocatable :: seed(:)
  integer :: s
  real*8 :: randnum
  real :: pi
  real*8, allocatable :: xi(:)
  real*8, allocatable :: yi(:)
  real*8, allocatable :: zi(:)
  real*8 :: xt
  real*8 :: yt
  real*8 :: zt
  character(len=50) :: configfilename
  character(len=100) :: wccommand
  
  configfilename = 'config-relaxed-8.dump'
  pi = dacos(-1.0d0)
  xlength=130 ! define box bounds - ENSURE THIS IS CONSISTENT WITH PARAMS IN SREV CONFIG FILE
  ylength=130 ! xlength, ylength, zlength should be length of simulation box in x,y,z directions
  zlength=130 !xlowerbound, ylowerbound, and zlowerbound should be same as in SREV config file
  xlowerbound = -65
  ylowerbound = -65
  zlowerbound = -65
  sigmawalkerwalker = 1.0 ! this is diameter of dye label in reduced units (ru) (1 ru = 10 nm)
  sigmachromatinchromatin = 1.0 ! this is diameter of SREV nucleosome face
  sigmawalkerchromatin = (sigmawalkerwalker+sigmachromatinchromatin)/2.0 
  walkers = 10000 ! this is number of dye labels to initialize

  CALL execute_command_line('rm walkerinput.in') ! remove any previously made dye input files
  wccommand = trim('wc -l < ' // trim(configfilename) // ' > wc.txt') ! get number of lines in SREV config file
  print *, wccommand
  CALL execute_command_line(wccommand)
  OPEN(unit=8, file='wc.txt')
  READ(8, *) linecount
  print *, linecount
  rewind(8)
  CALL execute_command_line('rm wc.txt')
  close(8)
  open(unit = 4, file = trim(configfilename))
  startline = 9
  do i = 1, startline
     if (i == 4) then
        read(4,*) npart ! read in number of SREV nucleosomes
     else
        read(4,*)
     endif
  enddo

  volume=xlength*ylength*zlength

  allocate (allx(npart), ally(npart), allz(npart)) ! allocate enough space to arrays meant for storing SREV nucleosome coords

  do i = 1, npart
     read(4,*) k,k,allx(i),ally(i),allz(i) ! read in SREV nucleosome coords from input SREV config file
  enddo

  close(4)
  
  print *, "finished reading in SREV config data"

  allocate(xi(walkers)) ! allocate enough space to arrays meant for storing initialized dye label coords
  allocate(yi(walkers))
  allocate(zi(walkers))

  print *, 'Starting random seed call'

  CALL RANDOM_SEED (size=s) ! using RANDOM_SEED & RANDOM_NUMBER() to generate the random number since some sites recommend that over RAND: https://gcc.gnu.org/onlinedocs/gfortran/RAND.html
  allocate(seed(s)) ! RANDOM_SEED has thing where seed must be array of some minimum size I won't know until I call RANDOM_SEED
  seed = 1

  open(5, file = 'walkerinput.in', status = 'new')
  CALL RANDOM_SEED(PUT=seed)
  do i = 1, walkers
     if (MOD(i,100) .eq. 0) then 
        print *, i ! meant to report progress of initialization script to user (ex: "100" means it is at loop i=100)
     endif
71   CALL RANDOM_NUMBER(randnum)  ! 71 is label for GOTO statements
     xt = randnum*xlength + xlowerbound
     CALL RANDOM_NUMBER(randnum)
     yt = randnum*ylength + ylowerbound
     CALL RANDOM_NUMBER(randnum)
     zt = randnum*zlength + zlowerbound
     do j = 1, npart ! check generated dye label coord against SR-EV beads for overlaps
        xr = xt - allx(j)
        yr = yt - ally(j)
        zr = zt - allz(j)
        xr = xr - xlength*ANINT(xr/xlength) ! adjust computed distance for periodic boundary conditions in simulation
        yr = yr - ylength*ANINT(yr/ylength)
        zr = zr - zlength*ANINT(zr/zlength)
        r = sqrt(xr**2 + yr**2 + zr**2)
        if (r < (0.9*sigmawalkerchromatin)) then
           go to 71
        endif
     enddo
     do k = 1, i-1 ! check generated dye label coord against already initialized dye labels
        xr = xt-xi(k)
        yr = yt-yi(k)
        zr = zt-zi(k)
        xr = xr - xlength*ANINT(xr/xlength) ! adjust computed distance for periodic boundary conditions in simulation
        yr = yr - ylength*ANINT(yr/ylength)
        zr = zr - zlength*ANINT(zr/zlength)
        r = sqrt(xr**2 + yr**2 + zr**2)
        if (r < (0.9*sigmawalkerwalker)) then
           go to 71
        endif
     enddo
     xi(i) = xt
     yi(i) = yt
     zi(i) = zt
     write(5,fmt="(a,1x,i7.1,1x,a,1x,f16.8,1x,f16.8,1x,f16.8,1x)") "create_atoms", 2, "single", xt,yt,zt ! write label coords to file for input into sim
  enddo

  close(5)

end program initializelabels
