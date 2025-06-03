c **************************************************************
c
c This file contains the subroutines:  metropolis
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku Hu
c
c **************************************************************

      subroutine metropolis(eol,currtem,acepta)

c===============================================================
c== SUBROUTINE FOR METROPOLIS UPDATE OF CONFIGURATIONS 	  ==
c==									 	  ==
c== CALLS: energy,addang,(rand),					  ==
c== dummy (function provided as argument)				  ==
c===============================================================


c------------GRSA-----------
      include 'INCL.H'

	real*8 e,t_col
	real KE,delta_i,delta_l,KE1,delta1,KE1_inter,KE2_inter,KE_ch1    ! parametros para el algoritmo CRO
	real minhits,theta,delta,difhits
	dimension vlvr1(mxrs),vlvr2(mxrs),c_vlvr(mxrs),c2_vlvr(mxrs)
	dimension vlvrv(mxrs)
	real chil1_vlvr(mxrs),chil2_vlvr(mxrs)
	integer desp1,desp2,Temp,cruza
	real time_i,eps,c1_eol,c2_eol,nval,betta
	real delta_ch1,delta_ch2,Einter,Einter1,Einter2


    	theta=15	
		minhits=20	
		eps=0.5	
        
		nval=0.8
        
		betta=0.8	
	
	delta_l=rand()	
	time_i = secnds(0.0)
c================================
c== Get Proposal configuration ==
c================================
c	collision type

	t_col  = rand()
    	
	if (t_col.LT.1.0) then	! collision on walls 
cc	write(*,*) ' collision on walls'
c================================
c== Decomposition ===
c================================
c		vlvrv= vlvr
		jv = idvr(1+int(nvr*rand()))  ! select var.
		vrol = vlvr(jv)  ! save old
        	e = rand()
        	dv = axvr(jv)*e
        	vlvr(jv) = addang(vrol,dv)
		
		enw = energy()
	        delta = eol
	        delta_i = enw 

c	si la molecula permanece estable durante un numero de hits
		if(abs(enw-eol).LE.eps) then
			numhits=numhits+1
			difhits=minhits-numhits
c   	realizamos perturbación circular		
			if (difhits.lt.theta) then
c				write(*,*) 'entro también aquí!!!'
c				write(*,*) 'perturbacion circular'				
				desp1=-86+int(mxrs*rand()*2)
				vlvr1=cshift(vlvr,desp1)
				desp2=-86+int(mxrs*rand()*2)
				vlvr2=cshift(vlvr,desp2)

				c_vlvr=vlvr	! guardamos la copia del original

				vlvr = vlvr1    ! asignamos la configuración propuesta y calculamos su energía
				c1_eol = energy()
				vlvr = vlvr2
				c2_eol = energy()

c	se verifica cual configuración es mejor					
				if (c1_eol.le.enw.and.c1_eol.le.c2_eol) then
					vlvr=vlvr1
					delta_i = energy()
c					write(*,*)'ganador',delta_i														
				else if (c2_eol.le.enw.and.c2_eol.le.c1_eol) then
					vlvr=vlvr2
					delta_i = energy()	
c					write(*,*)'ganador',delta_i			
				else
					vlvr=c_vlvr
					delta_i = energy()								
c					write(*,*)'ganador',delta_i
				endif
												
			endif
c			numhits=numhits+1
		else
			numhits=1
		endif
	
c================================
c== check acceptance criteria ===
c================================
	        if (delta+KE.GE.delta_i) then
	          eol = enw
	          acepta = acepta + 1.0d0
		KE=((delta+KE)-delta_i)*delta_l        
	        else
	          vlvr(jv) = vrol
	        endif

	else		! inter-melecular collisions 
		write(*,*) ' inter molecular collision'
c		
		c2_vlvr=vlvr
		jv = idvr(1+int(nvr*rand()))  ! select var.
		vrol = vlvr(jv)  ! save old
       		e = rand()
	      	dv = axvr(jv)*e
	     	vlvr(jv) = addang(vrol,dv)
		
c		se tienen que generar dos hijos y calcular sus energías

		cruza =1+int(mxrs*rand())

		do i=1,mxrs
        		if (i.le.cruza) then 					
				chil1_vlvr(i)=c2_vlvr(i)
				chil2_vlvr(i)=vlvr(i)
			else
				chil1_vlvr(i)=vlvr(i)
				chil2_vlvr(i)=c2_vlvr(i)
			endif	          	
	       	enddo

c		criterio de sistesis
		if (KE1_inter.ge.betta) then

        	write(*,*)'entra a sintesis...'

			delta=energy()											
			c_vlvr=vlvr
			vlvr=c2_vlvr
			delta_i=energy()
			vlvr=chil1_vlvr	
			delta_ch1=energy()
			
			if (delta+delta_i+KE1_inter+KE2_inter.ge.delta_ch1) then
				KE_ch1= (delta+delta_i+KE1_inter+KE2_inter)-delta_ch1
				vlvr=chil1_vlvr
				delta=delta_ch1			
			else	

			endif				

		else

c		IntermolecularIneffectiveCollision

			write(*,*)'entra a intermolecular...'
 
			delta=energy()											
			c_vlvr=vlvr
			vlvr=c2_vlvr
			delta_i=energy()
			vlvr=chil1_vlvr	
			delta_ch1=energy()
			vlvr=chil2_vlvr
			delta_ch2=energy()

			Einter1=delta+delta_i+KE1_inter+KE2_inter
			Einter2=delta_ch1+delta_ch2		
			Einter=Einter1-Einter2

				if (Einter1.ge.Einter2) then
					KE1_inter = Einter*delta_l
					KE2_inter = Einter*(1-delta_l)					

c					se verifica cual configuración es mejor					
					if (delta_ch1.le.delta_ch2) then
						vlvr=chil1_vlvr
						eol = energy()
c						write(*,*)'ganador',delta_i														
					else
						vlvr=chil2_vlvr
						eol = energy()								
c						write(*,*)'ganador',delta_i
					endif
															
				else
					vlvr=c_vlvr
					eol=energy()	
				endif

		end if

	  endif			

      return
      end
