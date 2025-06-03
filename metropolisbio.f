c **************************************************************
c
c This file contains the subroutines:  metropolis
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku Hu
c
c **************************************************************

      subroutine metropolisbio(eol,currtem,acepta)

c===============================================================
c== SUBROUTINE FOR METROPOLIS UPDATE OF CONFIGURATIONS 	  ==
c==									 	  ==
c== CALLS: energy,addang,(rand),					  ==
c== dummy (function provided as argument)				  ==
c===============================================================


c------------GRSA-----------
      include 'INCL.H'

	     real*8 e
       dimension aux_vlvr(mxvr)
       dimension vlvrBio(mxvr) 
        
c================================
c== Get Proposal configuration ==
c================================
	aux_vlvr = vlvr  ! save old
        call bioinspired(vlvrBio,enerBio)
c        write(*,*) "energyBio",enerBio
        
        !!Evaluate the Bioinspired Solution
        vlvr = vlvrBio        
        enw = energy()
        delta = enw - eol
	deltar = abs(delta)
	if (deltar.EQ.0) then
	 deltar = 1
	end if
	if (deltar.LT.deltamin) then
         deltamin = deltar
        else if (deltar.GT.deltamax) then
         deltamax = deltar
        end if
c================================
c== check acceptance criteria ===
c================================
        if (delta.LE.0.0d0) then
          eol = enw
          acepta = acepta + 1.0d0
        elseif (exp(-delta/currtem).GT.rand()) then
          eol = enw
          acepta = acepta + 1.0d0
        else
          vlvr = aux_vlvr !!!Return the old energy
        endif
      return
      end


