c **************************************************************
c This file contains the:  main (SINGLE PROCESSOR JOBS ONLY,
C                                FOR PARALLEL JOBS USE pmain)
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
C CALLS: init_energy,init_molecule
C CALLS TASK SUBROUTINE: anneal,canon,elp,minim,mulcan_par,
c                        mulcan_sim,partem_s, or regul
C CAN ALSO CALL MEASUREMENT ROUTINES: cnteny,contacts,helix,hbond,
C                                    outpdb,outvar,rgyr,
C                                    rmsinit and rsmdfun,zimmer
c $Id: main.f 334 2007-08-07 09:23:59Z meinke $
c **************************************************************
      
      program main

      include 'INCL.H'
      include 'INCP.H'
      include 'omp_lib.h'

      common/updstats/ncalls(5),nacalls(5)
      character*80 libdir, seqfile, varfile
      character grpn*4,grpc*4
      character(8) x1
      logical lrand,bgsposs
	  real inicio,final,namefile1
 

c      character*10 b(3)
c      integer date_time(8)
c      character*50 filename
c      real tiempo
c	  integer time

c =================================================== Energy setup

c            Directory for SMMP libraries
c     Change the following directory path to where you want to put SMMP
c     libraries of residues. 
      libdir='./SMMP/'

c      The switch in the following line is now not used.
      flex=.false.        ! .true. for Flex  / .false. for ECEPP

c     Choose energy type with the following switch instead ...
      ientyp = 0
c        0  => ECEPP2 or ECEPP3 depending on the value of sh2
c        1  => FLEX 
c        2  => Lund force field
c        3  => ECEPP with Abagyan corrections
c

      sh2=.true.         ! .true. for ECEPP/2; .false. for ECEPP3
      epsd=.false.        ! .true. for  distance-dependent  dielectric
                          !  permittivity

      itysol= 0    !  0: vacuum
                   ! >0: numerical solvent energy
                   ! <0: analytical solvent energy & gradients

      call init_energy(libdir)

c ================================================= Structure setup

      grpn = 'nh2' ! N-terminal group
      grpc = 'cooh'! C-terminal group

      iabin = 1  ! =0: read from PDB-file
                 ! =1: ab Initio from sequence (& variables)
      seqfile='EXAMPLES/1a13.seq'
      varfile='EXAMPLES/1a13F.var'
!       varfile = ' '
      
      ntlml = 0
      write (*,*) 'Solvent: ', itysol

c     Initialize random number generator.
	CALL SYSTEM_CLOCK(seed)
	call sgrnd(seed)
	call srand(seed)

c     Initialize random number generator.
c      call sgrnd(31433)
      
      if (itysol.eq.0.and.ientyp.eq.3) then
         print *,'Can not use Abagyan entropic corrections without '
         print *,'solvent term. '
         stop
      endif

      call init_molecule(iabin,grpn,grpc,seqfile,varfile)

c Decide if and when to use BGS, and initialize Lund data structures 
      bgsprob=0.90   ! Prob for BGS, given that it is possible
c upchswitch= 0 => No BGS 1 => BGS with probability bgsprob 
c 2 => temperature dependent choice 
      upchswitch=1
      rndord=.true.
      call init_lund
      if (ientyp.eq.2) call init_lundff
      if (ientyp.eq.3) call init_abgn
      

c ========================================  Add your task down here

	  
		call cpu_time(inicio)
			call anneal(energia)
		call cpu_time(final)
      			total = final - inicio
	
	namefile1=energia*10000
	write(x1,'(I8)')  int(namefile1)
    
      	write(*,*) 'ENERGY FOUND',energia,total
	open(16, file='./RESULTS/Energy'//x1//'.txt',status='new')
	write(16,*) 'ENERGY FOUND',energia,total
	close (16)

c ========================================  End of main      
       end
