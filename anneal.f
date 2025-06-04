c *************axvr*************************************************
c
c This file contains the subroutines:  anneal
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku Hu
c
c **************************************************************

      subroutine  anneal(energia)

C --------------------------------------------------------------
C PURPOSE: SIMULATED ANNEALING SEARCH OF LOWEST-POTENTIAL-ENERGY
C          CONFORMATIONS OF PROTEINS
C
C CALLS: addang,energy,metropolis,outvar,outpdb,rgyr,setvar,zimmer
C
C ---------------------------------------------------------------
C  SIMULATION WITH THE CHAOTIC FUNCTION
C ---------------------------------------------------------------

C------------------------------------------
C---------------------GRSA-----------------
C------------------------------------------

      include 'INCL.H'

      character(8) x2
      logical lrand
      parameter(lrand=.true.)
	logical bandera
	logical bandera1
	logical bandera2
	logical bandera3
	logical bandera4

        dimension vlvrm(mxvr)
	dimension idvrm(mxvr)
	dimension axvrm(mxvr)

        dimension vlvr_aux(mxvr)
	dimension idvr_aux(mxvr)
	dimension axvr_aux(mxvr)

	real*8 energia,e,bangl, energiapr
	real inicio,final,total
	real inicioT,namefile
	integer seed1,p,THREADS,ID,MaxThreads
        integer countbest
	integer, dimension(1:nvr) :: posvect02,resvect
	integer, dimension(1:nvr) :: countvect
 	real, dimension(1:nvr) :: normvect
      
      
        bandera = .true.
 	bandera1 = .true.
 	bandera2 = .true.
 	bandera3 = .true.
 	bandera4 = .true.
	val_pen = 100.0d0

c	seed=86456
	paro_recta = 0.000000000001d0
	paro_recta1 = 0.0000001d0
	recta = 0
	C = 0
	D = 0
	cnt = 1
	cnt2 = 0


	inicioT = secnds(0.0)
       nresi=irsml2(1)-irsml1(1) + 1     
c ==================
c == random start ==
c ==================
	
	!!!!Vector vlvr provide by file var (fragments)
        eol = energy()
	ymin = eol
 	epsilon = eol

	call annealident(countvect,jv2a,posvect02,resvect)
c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
c      do i=1,ntlml
c        call outpdb(i,11)
c      end do
c==========================
c==SIMULATING ANNEALING====
c==========================
      write(*,*) 'RECOCIDO SIMULADO'
	call cpu_time(start)
	
      write(*,*) 'Initial Energy', eol
c ====================================================
c == Simulating Annealing-Initial Parameters==========
c ====================================================
	
c ====================================================
c ========== Peptide Parameters ==============
c ====================================================	
          
	deltamin = 10.0d0
    
        currtem = 5990.5068258550d0 !Initial Temperature 1a13
        paro =  0.0000000126308167407501!Final Temperature (stop) 1a13	

	temp_aureo = currtem*0.618
	temp_aureo1 = temp_aureo*0.618
	temp_aureo2 = temp_aureo1*0.618
	temp_aureo3 = temp_aureo2*0.618
	temp_aureo4 = temp_aureo3*0.618

        alfa = 0.70
	blmax = 90.0d0

	bbeta = 1.00073521131733 !1a13 parameter

	temp1=0.1  !!!Reheat(First)
	temp2=1.9E-4 !!!Reheat(Second)
	temp3=1.9E-8 !!!Reheat(Third)
  	ban1=1
	ban2=1
	ban3=1
	ban4=1
	band5=1
	band6=1	
	band7=1

c =================================
c == Simulating Annealing Start ==
c =================================
      do while(currtem.GE.paro)
        propon = 0.0d0
        acepta = 0.0d0
        ycurr = 0.0d0
c =====================================
c == Metropolis Start  		 ==
c == Markov Chain (Increasing)	 ==
c =====================================
		numhits=1
        do while(propon.LE.NINT(blmax))
          propon = propon + 1.0d0
          ycurr = eol          
			
	if (ban2.EQ.1) then
          ban2 = 0
          write(*,*) "Bioinspired-JSOA"
        end if

	  call metropolisbio(eol,currtem,acepta)
c	  write(*,*) "Bioinspired-energy", eol


c =================================================
c == Store the local lowest-energy conformation  ==
c =================================================
          if (eol.LT.ycurr) then
            ycurr = eol
          end if
c =================================================
c == Store the global lowest-energy conformation ==
c =================================================
          if (eol.LT.ymin) then
            ymin = eol
            write(*,*) 'MINIMA:', currtem, ymin
c==================================================
c   Stochastical Equilibrium
c==================================================
			if ((currtem.LT.val_pen).AND. (cnt.LT.3)) then
				D = D + ymin
				C = C + cnt * ymin
				cnt = cnt + 1
c				write(*,*) 'leer las variables'
			end if

			if (cnt.EQ.3) then
		        recta = (((12 * C)-(6*(cnt-1)*D))/(cnt**3 - cnt))	
			cnt = 1	

			 if (recta.LT.paro_recta) then
c				write(*,*) 'bye!!'			
				currtem = currtem*0.618
				bangl= bangl+1
		         endif
			
			 if (bangl.EQ.2) then
			        currtem = paro
			 endif			  

			endif
c ==================================
c == Respaldo del punto minimo a  ==
c == a idvrm, axvrm, vlvrm	    ==
c ==================================
       	do i=1,nvr
		  idvrm(i)  = idvr(i)
		  iv        = idvr(i)
		  axvrm(iv) = axvr(iv)
		  vlvrm(iv) = vlvr(iv)
       	enddo
          end if
        enddo
c ========================
c == Termina Metropolis ==
c ========================

c ==============================
c == Disminuye la Temperatura ==
c ==============================
        currtem = alfa * currtem
	write(*,*) currtem,ymin
c ==============================
c Golden Ratio 
c ==============================

	if ((currtem.LT.temp_aureo).AND. (bandera)) then
		alfa = 0.75d0
		bandera = .false.
	else if ((currtem.LT.temp_aureo1).AND. (bandera1)) then
		alfa = 0.80d0
		bandera1 = .false.
	else if ((currtem.LT.temp_aureo2).AND. (bandera2)) then
		alfa = 0.85d0
		bandera2 = .false.
	else if ((currtem.LT.temp_aureo3).AND. (bandera3)) then
		alfa = 0.90d0
		bandera3 = .false.
	else if ((currtem.LT.temp_aureo4).AND. (bandera4)) then
		alfa = 0.95d0
		bandera4 = .false.
c		Reannealing con Golden Ratio
		currtem = currtem+currtem
	endif

      	if ((currtem.LE.temp2) .AND. (ban1.EQ.1) ) then
	    ban1 = 0
    	    write(*,*)bandera
	    write(*,*)'First REHEAT', currtem,ymin
	    currtem = currtem+currtem+currtem       	
        end if
	
	if ((currtem.LT.paro) .AND. (ban3.EQ.1) ) then
	    ban3 = 0
    	    write(*,*)bandera
	    write(*,*)'Residue REHEAT', currtem,ymin
	    currtem = temp1
c	    energiapr = ymin
        end if


c ==============================
c == CADA DE MARKOV CRECIENTE ==
c ==============================
        blmax = bbeta * blmax

      enddo
	energia = ymin
c ====================================
c == Termina el Simulated Annealing ==
c ====================================
	  namefile=energia*10000
	  write(x2,'(I8)') int(namefile) 	
          open(16, file='./RESULTS/Estruc'//x2//'.pdb',status='new')	  

c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
       do i=1,ntlml
         call outpdb(i,16)
       end do

       close(16)

      end
