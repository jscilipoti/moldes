      SUBROUTINE GENESTINTRAM (IFAM,MDV,NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,
     *NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,M,J) ! GEneración de ESTructuras 
C											   INTermedias RAMificadas
      PARAMETER (NSVA=20,NMAX=60000,NESTM=7100)
      IMPLICIT real*8 (A-H,O-Z)
      INTEGER MUD(NSVA,0:NMAX),ncomb(0:10,1001,5),
     *		NFAv(0:NESTM),K1maxv(0:NESTM)
	integer X3,X4,X2,X3min,X3max,X2max,Mv(0:NMAX),Kv(0:NESTM,3)
      common/estintram/MUD,ncomb,NFAv,K1maxv,Mv,Kv
c
c		MUD es la matriz donde se va almacenando la composición 
c		por subgrupos de cada nueva estructura intermedia generada.
c		Se generan J IMS's a partir de M meta-IMS's
C
	do N=1,10
		i=0
		do i1=0,N
			do i2=0,N-i1
				do i3=0,N-i1-i2
					do i4=0,N-i1-i2-i3
						i5=N-i1-i2-i3-i4
						i=i+1
						ncomb(N,i,5)=i1
						ncomb(N,i,4)=i2
						ncomb(N,i,3)=i3
						ncomb(N,i,2)=i4
						ncomb(N,i,1)=i5
					end do
				end do 
			end do
		end do
	end do
c
c	La matriz ncomb(10,1001,5) contiene todas las combinaciones de hasta 10 grupos a partir
c	de 5 grupos distintos. Por ejemplo, ncomb(7,263,2)=3 indicaría que la combi- 
c	nación No. 263 de 7 grupos a partir de 5 distintos contiene 3 veces el segundo.
c
      M=-1	! Para que no sea tenida en cuenta la IMS "vacia"
C			  que luego no podría evaluarse. Esta se genera después
c			  en molde2 y se agrega a las IMS ya aceptadas.	
	j=-1

	lugJ4=NGK4+1
	
	lugK3=lugJ4+NGJ4
	
	lugK3J2=lugK3+NGK3
	
	lugK2J3=lugK3J2+NGK3J2
	
	lugK1J4=lugK2J3+NGK2J3
	
	lugJ3=lugK1J4+NGK1J4
	
	lugK2=lugJ3+NGJ3
	
	lugK2J2=lugK2+NGK2
	
	lugK1J3=lugK2J2+NGK2J2
	
	lugJ2=lugK1J3+NGK1J3
	
	NGX4=NGJ4+NGK4
      
	NGX3=NGJ3+NGK3+NGK3J2+NGK2J3+NGK1J4 !número de grupos trivalentes
	
	NGX2=NGJ2+NGK2+NGK2J2+NGK1J3 ! número de grupos divalentes

c	NG* es el número de grupos corresp. al meta-grupo * selecc. por el usuario
	if (ifam.eq.5) then
		NFAmin = 0	
	else
		NFAmin = 2
	end if
	if (NGX3.eq.0.and.NGX4.eq.0)then !cant de grupos tri y tetravalentes = 0
		NFAmax=NFAmin !Si no se eligen grupos de ramificación sólo se pasa una
		              !vez por el bucle 1 (no hay ramificaciones)
	else
		if (ifam.eq.3) NFAmax=8 !El máx núm de grupos terminales permitidos
		if (ifam.eq.4) NFAmax=4
		if (ifam.eq.5) NFAmax=6
	end if	
	do NFA=NFAmin,NFAmax !bucle 1----
		if (NGX4.eq.0) then ! NGJ4+NGK4 = NGX4
			if (ifam.eq.5) then  
				X3min=NFA        
			else
				X3min=NFA-2 !todas las ramificaciones son debidas a grupos trivalentes
			end if
		else
			X3min=0
		end if
		if (NGX3.eq.0)then !si no se eligieron grupos trivalentes, sólo se
			X3max=0        !pasa una sola vez por el bucle 2 
		else
			if (ifam.eq.5) then
				X3max=NFA
			else
				X3max=NFA-2
			end if
		end if	
		do X3=X3max,X3min,-2 !bucle 2----
		  if (X3.eq.-1) then
			goto 61 !finaliza el ciclo actual del bucle 1
		  end if
		  if (ifam.eq.5) then
			X4=(NFA-X3)/2
		  else
			X4=(NFA-2-X3)/2
		  end if
!---------Esto se agregó para evitar que se generen SMS con un grupo terminal de más		  
		  imp= modulo(nfa,2)
		  if(NGX4.NE.0.and.ifam.NE.5.and.NGX3.EQ.0.and.imp.NE.0) goto 61
!-------------
		  if (ifam.eq.4) then
			J2max=6-X3-X4-NFA
		  else
			J2max=12-X3-X4-NFA
		  end if
c Se eliminan los casos (posibles para NFA=8) que pasarían el Xtotal
c máximo de 12 al completarse la SMS.
		  if (J2max.lt.0) then
			goto 71 !finaliza el ciclo actual del bucle 2
		  end if
c La meta-IMS se completará con distintos X2 (cantidades de grupos duales)
		  if (NGX2.eq.0)then !límite para bucle 8
			X2max=0
		  else
			X2max=J2max !determina el máx n° de grup divalentes
		  end if	
C Nivel K ascendente para la parte ramificada de las meta-IMS generadas
		  if (NGK4.eq.0) then !límite para bucle 3
			K4max=0
		  else
			K4max=X4
		  end if
		  if (NGJ4.eq.0) then
			K4min=X4
		  else
			K4min=0
		  end if
		  do K4=K4min,K4max !bucle 3 -----
			    iK4=NCOMBIN(K4,NGK4) !función
		        J4=X4-K4
			    if (NGK3.eq.0) then !Límites para bucle 4
				    K3max=0
			    else
				    K3max=X3
			    end if
			    if (NGK3.eq.NGX3) then
				    K3min=X3
			    else
				    K3min=0
			    end if
			    do K3=K3min,K3max !bucle 4
				    iK3=NCOMBIN(K3,NGK3)
				    if (NGK3J2.eq.0) then
		    		    K3J2max=0
				    else
					    K3J2max=X3-K3
				    end if
				    if (NGK3J2.eq.NGX3-NGK3) then
					    K3J2min=X3-K3
				    else
					    K3J2min=0
				    end if
				    do K3J2=K3J2min,K3J2max !bucle 5 (trivalence)
					    iK3J2=NCOMBIN(K3J2,NGK3J2)
					    if (NGK2J3.eq.0) then
						    K2J3max=0
					    else
						    K2J3max=X3-K3-K3J2
					    end if
					    if (NGK1J4+NGJ3.eq.0) then
						    K2J3min=X3-K3-K3J2
					    else
						    K2J3min=0
					    end if
					    do K2J3=K2J3min,K2J3max !bucle 6 (trivalence)
						    iK2J3=NCOMBIN(K2J3,NGK2J3)
						    if (NGK1J4.eq.0) then
							    K1J4max=0
						    else
							    K1J4max=X3-K3-K3J2-K2J3
						    end if
						    if (NGJ3.eq.0) then
							    K1J4min=X3-K3-K3J2-K2J3
						    else
							    K1J4min=0
						    end if
					        do K1J4=K1J4min,K1J4max !bucle 7 (trivalence)
						        iK1J4=NCOMBIN(K1J4,NGK1J4)
						        J3=X3-K3-K3J2-K2J3-K1J4
						        Kram=K4+K3+K3J2+K2J3+K1J4 
						        NJF=2*(J4+K1J4)+J3+K2J3   !¿por qué no usa todos los tri y tetravalntes?
						        Jram=J4+J3+K1J4+K2J3+K3J2
						        if (ifam.eq.3) then
							        NJF = NJF + 2
						        else if (ifam.eq.4) THEN
							        NJF = NJF + 2
							        Kram=Kram+1 !por el K1J2 del ACCH2	
							        Jram=Jram+1 !por el K1J2 del ACCH2	
						        end if		
c Los siguientes if eliminan los casos en los cuales aún completándose 
c con todos J2 hasta X total = 12, la meta-IMS no sería factible.
						        if (Kram.gt.NJF) then
							        if (2*Kram.gt.J2max+Jram+NJF) goto 81
						        else if (Kram.gt.Jram+J2max)then
							        goto 81
						        end if
c Para que el anillo de las estructuras cíclicas tenga un mínimo de 3 grupos:
						        if (ifam.eq.5) then
							        if (X3+X4.lt.3) then
								        X2min = 3 - X3 - X4
							        else
								        X2min = 0
							        end if				
						        else
							        X2min = 0
						        end if
! Hasta aquí ya están definidos todos los núcleos de ramificación de las estructuras.
							    do X2=X2min,X2max !bucle 8 
c Dado un X2, nivel K ascendente hasta el límite de factibilidad
			                        if (NGK2.eq.0) then
				                        K2max=0
				                    else
					                    K2max=X2
				                    end if
				                    if (NGK2.eq.NGX2) then
					                    K2min=X2
				                    else
					                    K2min=0
				                    end if
					                do K2=K2min,K2max !bucle 9
							            iK2=NCOMBIN(K2,NGK2)
								        K=Kram+K2
								        J2=X2-K2 !J2 es el J2 máximo p/ K2
								        if (K.gt.NJF) then
								            if (2*K.gt.J2+Jram+NJF) goto 91
								        else if (K.gt.Jram+J2)then
									        goto 91
								        end if
								        if (NGK2J2.eq.0) then !límites bucle 10
								            K2J2max=0
								        else
									        K2J2max=X2-K2
								        end if
								        if (NGK1J3+NGJ2.eq.0) then
								            K2J2min=X2-K2
								        else
									        K2J2min=0
								        end if
					  			        do K2J2=K2J2min,K2J2max !bucle 10 (dual valence)
					        			    iK2J2=NCOMBIN(K2J2,NGK2J2)
					        			    K=Kram+K2+K2J2 
					        			    J2=X2-K2 !J2 p/ eval.= J2+K2J2
				        			        if (K.gt.NJF) then
				        			    	    if (2*K.gt.J2+Jram+NJF) goto 101
				        			        else if (K.gt.Jram+J2)then
				        				        goto 101
						        	        end if
				                            if (NGK1J3.eq.0) then
        					                    K1J3max=0
		        		                    else
				        	                    K1J3max=X2-K2-K2J2
				                            end if
				                            if (NGJ2.eq.0) then
					                            K1J3min=X2-K2-K2J2
			                                else
					                            K1J3min=0
				                            end if
								            do K1J3=K1J3min,K1J3max !bucle 11 (dual valence)
									            iK1J3=NCOMBIN(K1J3,NGK1J3)
								                K=Kram+K2+K2J2+K1J3
									            J2P=X2-K2-K2J2-K1J3 !Number of pure J2
							                    if (K.gt.NJF+K1J3) then !¿por qué "+K1J3"?
								                    iJ=J2P+K2J2+Jram+NJF+2*K1J3
								                    if (2*K.gt.iJ) then 
									                    goto 111 !al final del bucle 11
								                    else
									                    K1max=(iJ-2*K)/2
								                    end if
							                    else 
								                    iJ=Jram+K1J3+K2J2+J2P
								                    if (K.gt.iJ)then
									                    goto 111
								                    else
									                    K1max=iJ-K
								                    end if
							                    END IF
							                    if (K1max.gt.NFA-0) then
								                    K1max=NFA-0	! para formar pre-SMSs
							                    end if
	M=M+1
	NFAv(M)=NFA-2 ! para formar pre-SMSs
	K1maxv(M)=K1max
	Kv(M,1)=K-NJF-K1J3
	Kv(M,2)=2*K-iJ
	Kv(M,3)=K-iJ
C	Ahora se generan las IMS concretas a partir de la meta-IMS número M
		do 08 i1=1,iK4 !bucle 12
			do 07 i2=1,iK3 !bucle 13
				do 06 i3=1,iK3J2 !bucle 14
					do 05 i4=1,iK2J3 !bucle 15
						do 04 i5=1,iK1J4 !bucle 16
							do 03 i6=1,iK2 !bucle 17
								do 02 i7=1,iK2J2 !bucle 18
									do 01 i8=1,iK1J3 !bucle 19
									J=J+1
									Mv(J)=M
									if (NGK4.ne.0) then
										if (K4.eq.0) then
											do i=1,NGK4
												MUD(i,j)=0
											end do
										else
											do i=1,NGK4
											MUD(i,j)= ncomb(K4,i1,i)
											end do
										end if
									end if
									if (NGJ4.ne.0) then
											MUD(lugJ4,j)=J4
									end if
									if (NGK3.ne.0) then
										if (K3.eq.0) then
											do i=1,NGK3
												MUD(lugK3+i-1,j)=0
											end do
										else
											do i=1,NGK3
									MUD(lugK3+i-1,j)= ncomb(K3,i2,i)
											end do
										end if
									end if
									if (NGK3J2.ne.0) then
										if (K3J2.eq.0) then
											do i=1,NGK3J2
												MUD(lugK3J2+i-1,j)=0
											end do
										else
											do i=1,NGK3J2
								MUD(lugK3J2+i-1,j)= ncomb(K3J2,i3,i)
											end do
										end if
									end if
									if (NGK2J3.ne.0) then
										if (K2J3.eq.0) then
											do i=1,NGK2J3
												MUD(lugK2J3+i-1,j)=0
											end do
										else
											do i=1,NGK2J3
								MUD(lugK2J3+i-1,j)= ncomb(K2J3,i4,i)
											end do
										end if
									end if
									if (NGK1J4.ne.0) then
										if (K1J4.eq.0) then
											do i=1,NGK1J4
												MUD(lugK1J4+i-1,j)=0
											end do
										else
											do i=1,NGK1J4
								MUD(lugK1J4+i-1,j)= ncomb(K1J4,i5,i)
											end do
										end if
									end if
									if (NGJ3.ne.0) then
											MUD(lugJ3,j)=J3
									end if
									if (NGK2.ne.0) then
										if (K2.eq.0) then
											do i=1,NGK2
												MUD(lugK2+i-1,j)=0
											end do
										else
											do i=1,NGK2
									MUD(lugK2+i-1,j)= ncomb(K2,i6,i)
											end do
										end if
									end if
									if (NGK2J2.ne.0) then
										if (K2J2.eq.0) then
											do i=1,NGK2J2
												MUD(lugK2J2+i-1,j)=0
											end do
										else
											do i=1,NGK2J2
								MUD(lugK2J2+i-1,j)= ncomb(K2J2,i7,i)
											end do
										end if
									end if
									if (NGK1J3.ne.0) then
										if (K1J3.eq.0) then
											do i=1,NGK1J3
												MUD(lugK1J3+i-1,j)=0
											end do
										else
											do i=1,NGK1J3
								MUD(lugK1J3+i-1,j)= ncomb(K1J3,i8,i)
											end do
										end if
									end if
									if (NGJ2.ne.0) then
											MUD(lugJ2,j)=J2P
									end if
 01									continue
 02								continue
 03							continue
 04						continue
 05					continue
 06				continue
 07			continue
 08		continue
c			 
								            end do !bucle 11
111										end do !bucle 10
101									end do !bucle 9
 91								end do !bucle 8
 							end do !bucle 7
 81						end do !bucle 6
					end do !bucle 5
				end do !bucle 4
		  end do !bucle 3
 71		end do !bucle 2
 61	end do !bucle 1
	if (NGK2.ne.0) then
		M=M+1
		NFAv(M)=0 ! para formar pre-SMSs
		K1maxv(M)=0
		do i=lugK2,lugK2+NGK2-1
			J=J+1
			do L=1,lugJ2
				MUD(L,j)=0
			end do
			MUD(i,j)=1
		end do
	end if
	if (NGK3.ne.0) then
		M=M+1
		NFAv(M)=1 ! para formar pre-SMSs
		K1maxv(M)=0
		do i=lugK3,lugK3+NGK3-1
			J=J+1
			do L=1,lugJ2
				MUD(L,j)=0
			end do
			MUD(i,j)=1
		end do
	end if
	if (NGK4.ne.0) then
		M=M+1
		NFAv(M)=2 ! para formar pre-SMSs
		K1maxv(M)=0
		do i=1,NGK4
			J=J+1
			do L=1,lugJ2
				MUD(L,j)=0
			end do
			MUD(i,j)=1
		end do
	end if
      return
	end

	integer function NCOMBIN(n,NG)
C
C	Esta función devuelve el número de combinaciones posibles para obtener n grupos
C	a partir de NG grupos diferentes. Por ej.: NCOMBIN(3,2)=4 --> AAA,AAB,ABB,BBB
C	NCOMBIN = (NG+n-1)!/(n!*(NG-1)!) = (NG-1+n)*(NG-2+n)*(NG-3+n)*...*n/(NG-1)!
	IF (NG.LT.2) THEN
		NCOMBIN=1
	ELSE
	ifact=1
	NUM=1
	do i=1,(NG-1)
		NUM=NUM*(NG+n-i)
		ifact=ifact*i
	end do
	NCOMBIN=NUM/ifact		
	END IF
      return
	end