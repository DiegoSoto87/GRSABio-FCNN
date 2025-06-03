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
C  SIMULACION CON LA FUNCION CAOTICA
C ---------------------------------------------------------------

      include 'INCL.H'
     


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

	real*8 energia,e,bangl
	real inicio,final,total,real,inicioT
	integer seed1,p,THREADS,ID,MaxThreads
 	bandera = .true.
 	bandera1 = .true.
 	bandera2 = .true.
 	bandera3 = .true.
 	bandera4 = .true.
	val_pen = 100.0d0

	seed=86456
	paro_recta = 0.000000000001d0
	recta = 0
	C = 0
	D = 0
	cnt = 1
	cnt2 = 0

c	call clock(inicio)


	inicioT = secnds(0.0)
      nresi=irsml2(1)-irsml1(1) + 1     
c ==================
c == random start ==
c ==================

c   generar población inicial
	
        MaxThreads = OMP_GET_MAX_THREADS()	     ! runtime environment returns 4
	write(*,*) 'Numero de hilos',MaxThreads	
	call    OMP_SET_NUM_THREADS(MaxThreads) ! set to 4 threads

	!$OMP PARALLEL
      if(lrand) then	
	!$OMP PARALLEL DO 
	do p=1,9
        	do i=1,nvr
         	iv=idvr(i)  ! provides index of non-fixed variable
        	 e = rand(int(inicioT))
        	 dv=axvr(iv)*e
	         vr=addang(pi,dv)

         	vlvr(iv)=vr
        	enddo

	write(*,*) 'iteración',p
	enddo
	!$OMP END PARALLEL DO
	write(*,*) 'iteración',p
      	end if		
	!$OMP END PARALLEL
	stop

      eol = energy()
	ymin = eol
 	epsilon = eol

c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
c      do i=1,ntlml
c        call outpdb(i,11)
c      end do


c==========================
c==	RECOCIDO SIMULADO  ==
c==========================
      write(*,*) 'RECOCIDO SIMULADO'
	call cpu_time(start)

c ====================================================
c == Parametros iniciales para el Recocido Simulado ==
c ====================================================

	currtem = 226624331780759000000000000000000000.0d0
	temp_aureo = 140053837040509000000000000000000000.0d0
	temp_aureo1 = 86553271291034600000000000000000000.0d0
	temp_aureo2 = 53489921657859400000000000000000000.0d0
	temp_aureo3 = 33056771584557100000000000000000000.0d0
	temp_aureo4 = 20429084839256300000000000000000000.0d0

      paro = 0.00000014573301692707600
      alfa = 0.70d0
      blmax = 360.0d0
	bbeta  =  1.000806052

	temp1=1.9E-4
  	ban1=1
	ban2=1
	ban3=1
	ban4=1
	band5=1
	band6=1	
	band7=1

c===========================


c =================================
c == Inicia el Recocido Simulado ==
c =================================
      do while(currtem.GE.paro)
        propon = 0.0d0
        acepta = 0.0d0
        ycurr = 0.0d0
c =====================================
c == Inicia Metropolis  		 ==
c == Cadena de Markov Creciente	 ==
c =====================================
        do while(propon.LE.NINT(blmax))
          propon = propon + 1.0d0
          ycurr = eol

c      	if ((currtem.LE.TempEA) .AND. (ban1.EQ.1) ) then
c		cuenta = cuenta+1 
c		ban1 = 0
c    	write(*,*)bandera
c		write(*,*)'REHEAT', currtem,ymin

c		currtem = currtem+currtem
c           end if

          call metropolis(eol,currtem,acepta)
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
				write(*,*) 'leer las variables'
			end if

			if (cnt.EQ.3) then
		       recta = (((12 * C)-(6*(cnt-1)*D))/(cnt**3 - cnt))	
			cnt = 1	

			 if (recta.LT.paro_recta) then
				write(*,*) 'bye!!'			
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
c =======================
c == Termina Metropols ==
c =======================

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
c		Reanealing con Golden Ratio
		currtem = currtem+currtem
	endif


c ==============================
c == CADA DE MARKOV CRECIENTE ==
c ==============================
        blmax = bbeta * blmax

      enddo
	energia = ymin
c ====================================
c == Termina el Simulated Annealing ==
c ====================================


c      open(11, file='./RESULTADOS/ESTRUCTURA.pdb')
c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
c      do i=1,ntlml
c        call outpdb(i,11)
c      end do

c      close(11)

      end
