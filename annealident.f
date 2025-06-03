c *************axvr*************************************************
c
c This file contains the subroutines:  anneal
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku Hu
c
c **************************************************************

      subroutine  annealident(countvect,jv2a,posvect02,resvect)
C --------------------------------------------------------------
C PURPOSE: IDENTIFY
C
C CALLS: none
C------------------------------------------

      include 'INCL.H'
 
	integer jv1, jv2a, resinit
	integer, dimension(1:nvr) :: resvect
	integer, dimension(1:nvr) :: posvect
	integer, dimension(1:nvr) :: posvect02  !!!!!lateral chain
	integer, dimension(1:nvr) :: countvect
	integer, dimension(1:nvr) :: identvect !!!!!logical identifier
c=================================
c Values for angles
c=================================
ccccccccccccccccccccccc
	jv1 = 1
	do while(jv1.le.nvr)	
	identvect(jv1) = 0
	posvect02(jv1) = 0
	jv1 = jv1 + 1
	enddo
c	write(*,*) countvect
c	stop
ccccccccccccccccccccccccccccc
		jv1= 0
		jv2a= 0
		jvcountaa= 1  !!!count residues

          do while(jv1.le.nvr)
		jv1=jv1 + 1
                vrol = vlvr(jv1)  ! save old in vrol
	  if (nmvr(jv1).eq.'psi') then 	!!!initial each residue	with psi
		jvcountaa = jvcountaa + 1
		dvpsi=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
           write(*,*) 'ANGLE PSI', dvpsi, jv1, jvcountaa
ccccc           vlvr(jv1) = dvpsi
	   else if (nmvr(jv1).eq.'omg') then   		
		dvomg=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
c           write(20,*) 'ANGLE OMG', dvomg, jv1, jvcountaa
           write(*,*) 'ANGLE OMG', dvomg, jv1, jvcountaa
ccccc           vlvr(jv1) = dvomg
	   else if (nmvr(jv1).eq.'x1') then
 		dvx1=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x1', dvx1, jv1, jvcountaa
           write(*,*) 'ANGLE x1', dvx1, jv1, jvcountaa
cccccc           vlvr(jv1) = dvx1
	   else if (nmvr(jv1).eq.'x2') then 		
		dvx2=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x2', dvx2, jv1, jvcountaa
           write(*,*) 'ANGLE x2', dvx2, jv1, jvcountaa
cccccc           vlvr(jv1) = dvx2
	   else if (nmvr(jv1).eq.'x21') then 		
		dvx21=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1	
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x21', dvx21, jv1, jvcountaa
           write(*,*) 'ANGLE x21', dvx21, jv1, jvcountaa
cccccc          vlvr(jv1) = dvx21
	   else if (nmvr(jv1).eq.'x22') then 		
		dvx22=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x22', dvx22, jv1, jvcountaa
           write(*,*) 'ANGLE x22', dvx22, jv1, jvcountaa
cccccc           vlvr(jv1) = dvx22
	   else if (nmvr(jv1).eq.'x3') then 		
		dvx3=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x3', dvx3, jv1, jvcountaa
           write(*,*) 'ANGLE x3', dvx3, jv1, jvcountaa
cccccc           vlvr(jv1) = dvx3
	   else if (nmvr(jv1).eq.'x31') then 		
		dvx31=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x31', dvx31, jv1, jvcountaa
           write(*,*) 'ANGLE x31', dvx31, jv1, jvcountaa
ccccc           vlvr(jv1) = dvx31
	   else if (nmvr(jv1).eq.'x32') then 		
		dvx32=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x32', dvx32, jv1, jvcountaa
           write(*,*) 'ANGLE x32', dvx32, jv1, jvcountaa
ccccc           vlvr(jv1) = dvx32
	   else if (nmvr(jv1).eq.'x4') then 		
		dvx4=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x4', dvx4, jv1, jvcountaa
           write(*,*) 'ANGLE x4', dvx4, jv1, jvcountaa
ccccc           vlvr(jv1) = dvx4
	   else if (nmvr(jv1).eq.'x5') then 		
		dvx5=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x5', dvx5, jv1, jvcountaa
           write(*,*) 'ANGLE x5', dvx5, jv1, jvcountaa
ccccc           vlvr(jv1) = dvx5
	   else if (nmvr(jv1).eq.'x6') then 		
		dvx6=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x6', dvx6, jv1, jvcountaa
           write(*,*) 'ANGLE x6', dvx6, jv1, jvcountaa
ccccc           vlvr(jv1) = dvx6
	   else if (nmvr(jv1).eq.'x61') then 		
		dvx61=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x61', dvx61, jv1, jvcountaa
           write(*,*) 'ANGLE x61', dvx61, jv1, jvcountaa
ccccc           vlvr(jv1) = dvx61
	   else if (nmvr(jv1).eq.'x62') then 		
		dvx62=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
		jv2a= jv2a + 1
		posvect02(jv2a) = jv1
		identvect(jv1) = identvect(jv1) + 1
c           write(20,*) 'ANGLE x62', dvx62, jv1, jvcountaa
           write(*,*) 'ANGLE x62', dvx62, jv1, jvcountaa
ccccc           vlvr(jv1) = dvx62
	   else if (nmvr(jv1).eq.'phi') then 		
		dvphi=vrol
		resvect(jv1) = jvcountaa
		posvect(jv1) = jv1
c           write(20,*) 'ANGLE PHI', dvphi, jv1, jvcountaa
           write(*,*) 'ANGLE PHI', dvphi, jv1, jvcountaa
ccccc           vlvr(jv1) = dvphi
	   else if (nmvr(jv1).eq.'pst') then 		
		dvpst=vrol
		resvect(jv1) = 0
		posvect(jv1) = jv1
c           write(20,*) 'ANGLE PST', dvpst, jv1
           write(*,*) 'ANGLE PST', dvpst, jv1
ccccc           vlvr(jv1) = dvpst
	   else if (nmvr(jv1).eq.'omt') then 		
		dvomt=vrol
		resvect(jv1) = 0
		posvect(jv1) = jv1
c           write(20,*) 'ANGLE OMT', dvomt, jv1
           write(*,*) 'ANGLE OMT', dvomt, jv1
c           vlvr(jv1) = dvomt
	   endif	
	   enddo
		
ccccccccccccc

	countvect = identvect

      end
