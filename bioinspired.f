c *************axvr*************************************************
c
c This file contains the subroutines:  bioinspired
c Strategies of Jumping Spider Algorithm Optimization 
c Copyright XXXX-XXXX  Hernán Peraza Vazquez
c                      
c
c **************************************************************
      subroutine bioinspired(vlvrBio,enerBio)

       include 'INCL.H'

       dimension aux_vlvr(mxvr)
       dimension o(10)
       dimension fit(10)
       dimension agent_energy(10)
       dimension agent_vlvr(10,mxvr)
       dimension positions(10,mxvr)
       dimension theBestVct(mxvr)
       dimension theWorstVct(mxvr)
       dimension vlvrBio(mxvr) 
       dimension Vm(10,mxvr)     
       integer i,j,r,r1,r2,c,m,e_max,e_min
       integer Agents,Pc,Pm !!Agents
       integer getBinary,Maxiteration
       integer startIteration
       real*8 e,cauchy,walk,maxWorst,minBest,grsa2ener

       Agents = 10
       Maxiteration = 50
       !!write(*,*) "Bioinspired"
c      Initialization
       aux_vlvr = vlvr
       grsa2ener = energy()

       do i=1, Agents		 	
            aux_vlvr = vlvr   
            jv = idvr(1+int(nvr*rand()))
            e = rand()
            dv = axvr(jv)*e 
            vlvr(jv) = addang(vrol,dv)
            agent_energy(i) = energy()
            agent_vlvr(i,:) = vlvr
            vlvr = aux_vlvr
       enddo
        !!!The best energy (minimum)
        e_min= minloc(agent_energy,DIM=1)
        theBestVct = agent_vlvr(e_min,:)  !!!The best vector
        !!!!The worst energy(maximum)
	e_max= maxloc(agent_energy,DIM=1)
        theWorstVct = agent_vlvr(e_max,:) !!!The worst vector
	
	!!!!Pheromone "o"
	do i=1, Agents
		minBest = agent_energy(e_min)
		maxWorst = agent_energy(e_max)
		fit(i) = agent_energy(i)
		o(i) = (maxWorst-fit(i))/(maxWorst-minBest)
	enddo
	
c     Bioinspired Algorithm
      gravity= 9.80665 !! %m/seg^2 
      c1=0.5
      c2=0.3
      vo=100
      positions = agent_vlvr !!!Solutions
      Vm=1
      do startIteration=1, MaxIteration  !!!Start iterations
      do r=1, Agents   !!!SearchAgents
      if (rand().lt.c1) then !!! Attack, persecution and jumping on the prey
      ale1=rand()
      	if (rand().lt.c1) then
      	ale2=rand()
      	 !! Jumping on the prey represented by equation of projectile motion.
      	 radians_value= (90*pi*rand())/180
      	 Vm(r,:) = (positions(r,:)*tan(radians_value)) - 
     # 	 ((gravity*positions(r,:)**2)/(2*vo**2)*(cos(radians_value)**2)) 
      	else
      	 !! Persecution represented by the uniformly accelerated rectilinear motion.
      	 bander=1      	 
      	 do while(bander.eq.1) 
		r1 = idvr(int(1+(int(Agents-1)*rand()))) !!!!!!
		if (r.ne.r1) then
			bander=0	
		endif
	 end do
	Vm(r,:) = 0.5*(positions(r,:) - positions(r1,:))
	endif
	
      else !!!Searching for prey
	if (rand().lt.c1) then  !!! Global search
		c=1; m=0
        	cauchy = c*tan(pi*(rand()-0.5)) + m
        	Vm(r,:) = theBestVct + (theBestVct-theWorstVct)*cauchy 
        else                    !!! Local search
        	walk= -2 + 4*rand() !!!! -2 < d < 2 Uniformly distributed pseudorandom numbers
                e= rand(); !!! Normally distributed pseudorandom numbers.
                Vm(r,:)= theBestVct + walk*(0.5-e);
        endif
      endif	!!END !! Attack, persecution and jumping on the prey
     
      if (o(r).le.c2) then  !!Jumping spider pheromone rates(o)
      	bander=1
      	do while(bander.eq.1)
      	       r1 = idvr(int(1+(int(Agents-1)*rand()))) !!!!!!
               r2 = idvr(int(1+(int(Agents-1)*rand()))) !!!!!!
      		if (r1.ne.r2) then
      		bander=0
      		endif
      	end do
      	
      	if (rand().lt.c1) then
      		getBinary=0
      	else
      		getBinary=1
      	endif     	
      	Vm(r,:)=theBestVct+(positions(r1,:)-
     # 	((-1)**getBinary)*positions(r2,:))/2
      endif  !!!END !!!!Jumping spider pheromone rates (o)
      !!!!!!!!!SearchAgents concluded!!!!!!!!!!!!
      
      !!!! Evaluate new solutions
      aux_vlvr=vlvr !!!Save th global Solution
      vlvr=Vm(r,:)
      Fnew=energy() !!!Vm(r,:) is evaluated
      vlvr=aux_vlvr !!!Return the global Solution
      
      if (Fnew.le.fit(r)) then   !!!!!!if Fnew <= Fitness(r)
        positions(r,:)= Vm(r,:)
        fit(r)= Fnew
      endif
      if (Fnew.le.minBest) then
         theBestVct= Vm(r,:)   
         minBest= Fnew;           
      endif
      
      enddo  !!!SearchAgents (END)
                     
      !!!!Update the worst vector and its energy (Vector)
      e_max= maxloc(fit,DIM=1)
      theWorstVct = positions(e_max,:) !!!The new worst vector
      !!!!!Obtain the worst vector's energy 
      aux_vlvr=vlvr
      vlvr=theWorstVct
      maxWorst=energy()
      vlvr=aux_vlvr
      !!!Update the Pheromons "o" 
      do i=1, Agents
      	o(i) = (maxWorst-fit(i))/(maxWorst-minBest)
      enddo !!!!!!!!END Update the Pheromons "o"
      
      enddo !!!!!!END iterations
      
      enerBio = minBest
      vlvrBio = theBestVct
      
      end
