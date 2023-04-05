module Design
contains
 subroutine STRUCTURE_GENERATOR (JIST,MS,NC,NGSSV,NGSDV,SALIR,PROP)
!-------------------------------------------------
!   Descripci�n
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
    USE GRPS 
    use StructuresDesign
    USE PropertiesData
    use constantes
    use Input
    use Evaluation
    !USE PropertiesData
    !USE PROPI
    parameter (NESTM=7100,NA=150,NSVA=20,NMSV=20,NSCM=10,NAA=4)
    implicit real*8(A-H,O-Z)
    integer,intent(inout)::MS(NCOM,DiffStructGroups,2)
!...VARIABLES DE ENTRADA
    integer,intent(in)::NC,NGSSV(NSVA,5),NGSDV(NSVA,5)
     ! real*8,intent(in)::TEMP1
    logical,intent(in)::PROP
!...VARIABLES DE SALIDA
    integer,intent(out)::JIST
    logical,intent(out)::SALIR   
!...VARIABLES INTERNAS Y COMMONS
    type(Compound),pointer::recorreSolutes
    type(CalculatedProperties),pointer::ListProperties
    type(CalculatedPerformance),pointer::ListPerformance      
      !type(FinalStructure),pointer::IMSs !queda pendiente mdificar el c�digo para agregar esto
!   Integers
    integer::NGSDV1(NSVA,5),NUD(NSVA,NMAX),NUDC(NMSV),IIMS(10,2),&
        NUR(NA,NMAX),NUS(NSVA,NMAX),NSV1(3,NMAX),NSV2(3,NMAX),&
        NGSSV1(NSVA,5),nIMS(0:NMAX),Kpv(NMAX,3),nsmsv(0:nmax),&
        nsm1v(0:nmax),iterm(nmax,10,2),niterm(nmax),Mprv(NMAX),&
        Mprev(0:NMAX),H,HSR,ims_test(2,NMAX),NUD1(NSVA,NMAX),MSIN(NSCM,2),&
        NUD2(NSVA,NMAX),itermp(nmax,nscm),K1maxpv(NMAX),TERM(NMAX,NSCM,2)
    integer::bulky(NESTM,DiffStructGroups,2)

!   Characters
    character*1 RESP
    character*21 FAMILY(6)
!   Logical
    logical SOLDES,CONTAR,SOLACE,SMS,sm1,rep,NOACE,INV,FinalStruct
!   Commons
    COMMON/IUS/IMS(NSCM,2)
    COMMON/PUNSUB/NPUNT(NA),NGRUP(NA),NCANT
    COMMON/CONT/ PMG(NMG),DELTC(NMG),DELV(NMG),DELPC(NMG),DELVI(NMG),DELN(NMG),DELTCA(NMG)
    common/estintram/MUD(NSVA,0:NMAX),ncomb(0:10,1001,5),NFAv(0:NESTM),K1maxv(0:NESTM),Mv(0:NMAX),Kv(0:NESTM,3)
    common/NGX/NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1
    COMMON/EVAL/INV
    COMMON/EXTDIS/PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
    DATA family /'      aromatic groups','     molecular groups',&
                    3*'  dual valence groups','single valence groups'/

!...Sentencias   
    nullify(FMSs)
    nullify(recorreSolutes)
    allocate(recorreSolutes)
    nullify(ListPerformance)
    allocate(ListProperties)
    SALIR=.FALSE.
    NOACE=.FALSE.
    NUS(:,:)=0
    temp1 = InputProblem01%T
    mop = InputProblem01%mop
    contador=0

!
!	 Aqui empeze a corregir	(Liliana)
!        GENERATION OF INTERMEDIATE STRUCTURES
!
	if (IFAM/=2) then
		write (6,"(' ','***** GENERATING INTERMEDIATE STRUCTURES...')")
	else
		write (6,"(' ',/,' ***** EVALUATING SINGLE GROUP SOLVENTS...')") 
	endif
!
	if(ifam<=2) then	! Para arom�ticos o single substances
	    FinalStruct=.False.
	    call GENESTINT (0,MDV2,NUD2,JISL) !estructuras alif�ticas de las ramificaciones
		call alif_fact (ngdv,ngsdv,mdv1,mdv2,nud2,jisl,ims_test) !elimina las estructuras intermedias alif�ticas no factibles
		call GENESTINT (IFAM,MDV1,NUD1,J) !estructuras arom�ticas
		JEI = J !IMPRIME EL NUMERO DE INTERMIEDOS GENERADOS
		JIST=J
		JIS=0
		IF (IFAM.EQ.2) THEN !single substances
			JIS = JIST
		ELSE !arom�ticos
			DO 1205 J=1,JIST !ac� JIST es el n�mero de IMS's
				SOLDES = .FALSE.
				DO 131 M=1,NSCM !NSCM es un par�metro (=10)
					DO 132 K=1,2
						MS(NC,M,K) = 0
  132				CONTINUE
  131			CONTINUE
				ISF=0
				DO 133 I=1,MDV1 !bucle que guarda en MS la estructura intermedia
					IF (NUD1(I,J).EQ.0) cycle
				    ISF=ISF+1 ! cant de subgrup q participan en el dise�o
					MS(NC,ISF,1) = NGDV(I)  !n�mero identificatorio 
					MS(NC,ISF,2) = NUD1(I,J) !cantidad
					DO 134 K=1,5
						NGSDV1(ISF,K)=NGSDV(I,K)
  134				CONTINUE
  133			CONTINUE
				CALL AROMI (IPAREQ,MOP,ISF,NGSDV1,SOLDES,natot)
 				IF (.NOT.(SOLDES)) THEN 
  200				PMEI=0
					IF(.NOT.PROP)THEN
					    call gc_special_groups(ms(3,:,:),bulky,nbulky)
                        contar=.true.
                        recorreSolutes => InputProblem01%MixtureInput%Solutes
                        it=1					    
					    do while (contar)
					        ms(1,:,:) = recorreSolutes%Formula
					        BPSolute = recorreSolutes%BoilingPoint
                            call Evaluate_Mixture(.True.,FinalStruct,MS,NC,temp1,&
							BPSolute,ifam,NPSWIP,NPSBI,contar,ListPerformance,it)
                            recorreSolutes => recorreSolutes%next
                            if(.not.associated(recorreSolutes)) exit
        	                    it=it+1
                        enddo    
            	        if(contar)then
                            jis = jis + 1
                                !call Incorporate_Solvent (mop,ms(3,:,:),ListProperties,ListPerformance,FMSs)
                        endif
					else
					    CONTAR=.True.
				    ENDIF
					IF (CONTAR) THEN
						DO 143 I=1,MDV1
						    NUD(I,JIS)=NUD1(I,J)
 143					CONTINUE
				    END IF
					IF (JIS.EQ.NMAX) THEN
 1900 				    continue
						WRITE (6,880) NMAX,FAMILY(IFAM)
						READ (5,890) RESP
						IND1 = INDEX('CEce',RESP)
						IF (IND1.EQ.0) THEN
						    GO TO 1900
						ELSE IF ((RESP.EQ.'E').OR.(RESP.EQ.'e')) THEN
						    SALIR=.TRUE.
						    GO TO 1111
						ELSE
						    GO TO 1910
						END IF
				    END IF
		        END IF
1205		CONTINUE !se finalizan de generar la IMS

1910		WRITE (6,2050)

		END IF !arom�ticos
        
        CALL COMB_TERM (MSV,NUS,KST,JST,IST)
        CALL COMB_INTER (JISL,NUR,KSR,JSR,ISR,HSR)
		IST = K !Revisar esto!!!
		IZ = 0
		JIST = 0
		NARES = 1
		IFAM1 = IALAR(1)
		FinalStruct=.True. 

1601	DO 1600 J=1,JIS !generador de SMS para ifam <= 2
! genera todas las sms posibles a partir de la ims n� J		
			DO 1548 M=1,NSCM
				DO 1548 K=1,2
					MS(NC,M,K) = 0
1548		CONTINUE
            ims(:,:)=0
			ISI = 0
			DO 1605 I=1,MDV !copia al vector IMS la IMS guardada en NUD (:,J)
				IF (NUD(I,J).EQ.0) cycle
				ISI=ISI+1 !ISI=n�mero de subgrupos distintos en la sms
				IMS(ISI,1) = NGDV(I) !n�mero identificatorio
				IMS(ISI,2) = NUD(I,J)!cantidad
				DO 144 K=1,5
					NGSDV1(ISI,K) = NGSDV(I,K)
 144				CONTINUE
1605			CONTINUE
			KTEST=0
			MTEST=0
			JTEST=0
			DO 1145 I=1,ISI
				MTEST = IMS(I,2)*NGSDV1(I,1)+MTEST
				JTEST = IMS(I,2)*NGSDV1(I,2)+JTEST
				KTEST = IMS(I,2)*NGSDV1(I,3)+KTEST
1145		CONTINUE
			IF (IFAM.EQ.1) THEN !aromatic solvents
			    CALL AROMN (IPAREQ,ISI,NAC,NACCH,NACCH2,NATOT) !Averigua la cant de grupos terminales que deben agregarse
			    IF (NATOT.EQ.0) THEN !El sig bloque de sentencias condicionales asigna los l�mites para los
				    IIKST = 1        !bucles que agregan ramificaciones alif�ticas (IIKSR y IKSR) y grupos 
			        IKST = 1         !terminales (IIKST y IKST).
			        IIKSR=2
			        IKSR=1
			    ELSEIF (NATOT.EQ.1) THEN !por qu� IKSR toma estos valores?
				    IIKST = 1
			        IKST = MSV
			        IIKSR = 1
			        IKSR = KSR
			    ELSEIF (NATOT.EQ.2) THEN
			        IIKST = 1
			        IKST = KST					    
			        IIKSR = KSR + 1 !desde la primer combinaci�n de dos estructuras
			        IKSR = JSR !hasta la �ltima combinaci�n de dos estructuras
			    ELSEIF (NATOT.EQ.3) THEN
			        IIKST = KST + 1
      			    IKST = JST
      			    IIKSR = JSR + 1
      			    IKSR = ISR
      			ELSEIF (NATOT.EQ.4) THEN 
      			    IIKST = JST + 1
      			    IKST = IST
      			    IIKSR = JSR + 1
      			    IKSR = HSR
      			else
      			    cycle
                ENDIF
			ELSE IF (IFAM.EQ.2) THEN !single substance groups
				IKST = 1
				NATOT = 0
			END IF
!Bucle que agrega cadenas alif�ticas seg�n la cantidad de ramificaciones que
!el anillo permita. La primer iteraci�n (IIKSR-1) no agrega ninguna cadena alif�tica
      		DO 2625 H=IIKSR-1,IKSR
      			ISF=ISI
!El siguiente bucle elimina la ramificaci�n alif�tica 
!que pudo cargarse en la iteraci�n anterior
        	  	DO L=1,10
			    	IIMS(L,1) = IMS(L,1)
					IIMS(L,2) = IMS(L,2)
			  	ENDDO
!se agrega(n) la(s) ramificaci�n(es) alif�tica(s)
      		  	IF (H.GE.IIKSR) THEN
      		    	nCH2=0
           	    	DO M=1,JISL
      			    	IF (NUR(M,H).EQ.0) cycle
      			    	NUDC(:)=NUD2(:,M)
     		            CALL Ingresar_EstInter (MDV1,IIMS,NUDC,ISF,NGDV,NUR(M,H),iCH2)
                     	ISF=ISF+1
                     	nCH2=nCH2+iCH2
      			    ENDDO
!La siguiente sentencia impide que se contin�e con una IMS cuando el enlace entre el grupo (AC)
!y el grupo (CH2) es inevitable
      			  	if(nCH2.gt.(2*NACCH+NACCH2))cycle
                ENDIF
!Bucle para agregar todas las combinaciones posibles de grupos terminales seleccionados
				DO 1625 K=IIKST,IKST 
					SOLDES = .FALSE.
					DO 1777 L=1,10
						MS(NC,L,1) = IIMS(L,1)
						MS(NC,L,2) = IIMS(L,2)
1777				CONTINUE
					IF (NATOT.EQ.0) THEN !si no debe agregarse ning�n GrupTerm
						IST = ISF !cantidad de grupos distintos en la estructura
					ELSE IF (NATOT.EQ.1) THEN !si debe agregarse un GrupTerm
						CALL INGRESAR_TERMINAL (MS,ISF,NGSV(K),1) !!!
						IST = ISF + 1 !Cantidad de grupos distintos en la estructura
						IF (IFAM.EQ.1.and.H.lt.IIKSR)CALL AROM1 (IPAREQ,IST,NACCH,NACCH2,SOLDES)
					ELSE IF (NATOT.GE.2) THEN
						NSV = 0
						DO 150 I=1,MSV
							IF (NUS(I,K).EQ.0) cycle
					    	NSV=NSV+1
							NSV1(NSV,K)=NGSV(I)	! n�mero identificatorio de subgrupo
							NSV2(NSV,K)=NUS(I,K)  ! cantidad
							NGSSV1(NSV,1)=NGSSV(I,1) !propiedades combinat.
							NGSSV1(NSV,2)=NGSSV(I,2)
							NGSSV1(NSV,3)=NGSSV(I,3)
 150					CONTINUE
                    	DO I=1,NSV
                        	CALL INGRESAR_TERMINAL (MS,ISF+I-1,NSV1(I,K),NSV2(I,K)) !!!!
                    	ENDDO
                    	IST=ISF+NSV
						IF(IFAM.EQ.1.and.H.lt.IIKSR) THEN
							CALL AROM1 (IPAREQ,IST,NACCH,NACCH2,SOLDES) !controla casos especiales
							CALL AROM2 (K,NACCH,NSV,NSV2,NGSSV1,SOLDES)
						END IF
					END IF
				! CALL MOSTRAR_COMPONENTE (IST) !Subrutina para checkear por consola las estructuras generadas
					IF(.NOT.(SOLDES)) THEN

                    	call gc_special_groups(ms(3,:,:),bulky,nbulky)
                    	call Evaluate_Pure (FinalStruct,ms(3,:,:),temp1,ifam,NSEW,solace,ListProperties)
                    	recorreSolutes => InputProblem01%MixtureInput%Solutes
                    	it=1	
                    	nullify(ListPerformance)
						do while (solace)
					    	ms(1,:,:) = recorreSolutes%Formula
					    	BPSolute = recorreSolutes%BoilingPoint
                        	call Evaluate_Mixture(.True.,FinalStruct,MS,NC,temp1,BPSolute,ifam,NSWIP,NSBI,solace,ListPerformance,it)
                        	recorreSolutes => recorreSolutes%next
                        	if(.not.associated(recorreSolutes)) exit					
                        	it=it+1
                    	enddo    
                    
                    	if(solace)then
                        	jist = jist + 1
                        	call Incorporate_Solvent (mop,ms(3,:,:),ifam,ListProperties,ListPerformance)
                    	endif 

						IF (JIST.GE.NMAX) THEN
1940					    WRITE (6,930) NMAX,FAMILY(6)
							READ (5,890) RESP
							IND1 = INDEX('CEce',RESP)
							IF (IND1.EQ.0) THEN
								GO TO 1940
							ELSE IF ((RESP.EQ.'E').OR.(RESP.EQ.'e')) THEN
								SALIR=.TRUE.
								GO TO 1111
							ELSE
								GO TO 950
							END IF
				    	END IF
					END IF
1625			CONTINUE
2625        CONTINUE
1600	CONTINUE
		nfsg=nsew+nswip+nsbi
	else ! para alif�ticos, mixtos o c�clicos (ifam 3,4 y 5)
!			que ahora usan Genestintram en lugar de genestint
!
!      		ntcp: N�mero de grupos terminales para completar una pre-SMS			
		IF (IFAM.eq.4) THEN
			ntcp=1
		ELSE
			ntcp=2 ! Para ifam 3 o 5
		END IF

		CALL GENESTINTRAM (IFAM,MDV,NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,&
     						NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,M,J)
		JEI = J
		MEI = M
!
!         GENERATION OF THE MS SET FOR THE INTERMEDIATE STRUCTURES
!
		JIST=J
		JIS=0
		NIMS(0)=0
		Mpre=0
		Npre=0
		FinalStruct=.False.
		NPSEW=0
		NPSBI=0
		NPSWIP=0
		ns=0
		n1=0
		DO 205 J=1,JIST !ac� JIST es el n�mero de IMS's
			SOLDES = .FALSE.
			sms = .FALSE.
			sm1 = .FALSE.
			DO 31 M=1,NSCM !NSCM = 10 (par�metro)
				DO 32 K=1,2
					MS(NC,M,K) = 0
 32					CONTINUE
 31				CONTINUE
			if (NFAv(Mv(j)).eq.-2) then ! ya es SMS (c�cl. no ram.)
				sms = .TRUE.
				ikst=1
				do L=1,msv !cantidad de subgrupos seleccionados por el usuario
					itermp(1,L)=0
				end do
				Mpre= Mpre + 1
				Mprv(ikst)=Mpre			
			else if (NFAv(Mv(j)).le.0) then	! (NFAv = 0 o -1)
				ikst=1
				do L=1,msv
					itermp(1,L)=0
				end do
				Mpre= Mpre + 1
				Mprv(ikst)=Mpre			
				K1maxpv(Mpre)=K1maxv(Mv(j))
				Kpv(Mpre,1) = Kv(Mv(j),1)
				Kpv(Mpre,2) = Kv(Mv(j),2)
				Kpv(Mpre,3) = Kv(Mv(j),3)
				if (NFAv(Mv(j)).eq.-1) then
					sm1 = .TRUE.
				end if
			else ! NFA>0 p/formar pre SMS
				if(Mv(j).ne.Mv(j-1)) then
					ikst=0
					if (NGK1.eq.0) then
						K1max=0
					else
						K1max=K1maxv(Mv(j))
					end if
					if (NGK1J2+NGM1.eq.0) then
						K1min=K1maxv(Mv(j))
					else
						K1min=0
					end if
					do K1=K1min,K1max
						if (NGK1J2.eq.0) then
							K1J2max=0
						else
							K1J2max=NFAv(Mv(j))-K1
						end if
						if (NGM1.eq.0) then
							K1J2min=NFAv(Mv(j))-K1
						else
							K1J2min=0
						end if
						do K1J2=K1J2min,K1J2max
							M1=NFAv(Mv(J))-K1-K1J2
							if(Kv(Mv(j),1)+K1+K1J2.gt.0) then
								if(Kv(Mv(j),2)+2*K1+K1J2.gt.0) goto	807
								K1max=-(Kv(Mv(j),2)+2*K1+K1J2)/2
							else
								if(Kv(Mv(j),3)+K1.gt.0) goto 807
								K1max=-Kv(Mv(j),3)-K1
							end if
							if (K1max.gt.ntcp) then
								K1max=ntcp	! para completar las SMSs
							end if
							Mpre= Mpre + 1
							K1maxpv(Mpre)=K1max
							Kpv(Mpre,1) = Kv(Mv(j),1)+K1+K1J2
							Kpv(Mpre,2) = Kv(Mv(j),2)+2*K1+K1J2
							Kpv(Mpre,3) = Kv(Mv(j),3)+K1
! Ahora se forman todas las combinaciones concretas, correspondientes a la actual
!"meta-pre-terminaci�n" para la IMS concreta No. j, y v�lidas tambi�n para todas 
!	  las otras IMS con el mismo Mv(j) (o sea, hijas de la misma meta-IMS)
							iK1=NCOMBIN(K1,NGK1)
							iK1J2=NCOMBIN(K1J2,NGK1J2)
							do i1=1,iK1
								do i2=1,iK1J2
									ikst=ikst+1
									Mprv(ikst)=Mpre
									if (NGM1.ne.0) then
										itermp(ikst,MSV)=M1
									end if
									if (NGK1.ne.0) then
										do L=1,NGK1
											itermp(ikst,L)=ncomb(K1,i1,L)
										end do
									end if
									if (NGK1J2.ne.0) then
										do L=1,NGK1J2
											itermp(ikst,NGK1+L) = ncomb(K1J2,i2,L)
										end do
									end if
								end do
							end do
						end do
807						end do
				end if
			end if
!			ISF = Numero de subgrupos diferentes en la IMS
!			IST = Numero de subgrupos diferentes en cada SMS
!            
! El siguiente ciclo es para formar y evaluar IKST preSMS a partir de la IMS J
			DO 35 K=1,IKST
				Npre = Npre + 1 ! Number of pre-Final Solvents Generated
				ISF=0
				DO 33 I=1,MDV
					IF (MUD(I,J).NE.0) THEN
						ISF=ISF+1
						MS(NC,ISF,1) = NGDV(I)
						MS(NC,ISF,2) = MUD(I,J)
!						DO 34 K=1,5
!							NGSDV1(ISF,K)=NGSDV(I,K)
!  34						CONTINUE
					END IF
  33			CONTINUE
				if(sms) then
					contar=.TRUE.
					soldes=.FALSE.
					pmei=0
					goto 115
				end if
				do L=1,msv
					if (itermp(K,L).ne.0) then
						ISF=ISF+1
						MS(NC,isf,1) = NGSV(L)
						MS(NC,isf,2) = itermp(k,L)
					end if
				end do
!c
!c			ISF= Cantidad de subgrupos distintos que forman parte de la 
!C			pre-SMS k, correspondiente a la estructura intermedia J
!c
!c			Soldes = true para IMS que no satisfacen el criterio de factibilidad
!c			false para IMS factibles
 115			IF (.NOT.(SOLDES)) THEN
					CONTAR = .FALSE.
				    if(sms) goto 117
!						Hasta aqui corregi (Liliana)
!
!				    CALLING OF UNIFAC SUBROUTINES
!
!					COMPUTATION OF pre- solvent PROPERTIES
                    NPSWIPprev=NPSWIP
					IF(.NOT.PROP.and.ifam/=4.and.model/=3)THEN 
						
                        call gc_special_groups(ms(3,:,:),bulky,nbulky)
                        if (InputProblem01%mop==4)then  !PROVISORIO PARA REACCIONES!!!!
					       ! BPSolute = recorreSolutes%BoilingPoint
                            call Evaluate_Mixture(.True.,FinalStruct,MS,NC,temp1,&
							BPSolute,ifam,NPSWIP,NPSBI,contar,ListPerformance,it)                            
                        else
                        contar=.true.
                        recorreSolutes => InputProblem01%MixtureInput%Solutes	
                        it=1
			     		do while (contar)
					        ms(1,:,:) = recorreSolutes%Formula
					        BPSolute = recorreSolutes%BoilingPoint
                            call Evaluate_Mixture(.True.,FinalStruct,MS,NC,temp1,&
							BPSolute,ifam,NPSWIP,NPSBI,contar,ListPerformance,it)
                            recorreSolutes => recorreSolutes%next
                            if(.not.associated(recorreSolutes)) exit
                            it=it+1
                        enddo                                 
						endif
					else
					    CONTAR=.True.
				    ENDIF
  					IF(.NOT.PROP.and.model/=3)THEN
                        IF (.not.contar.and.ifam==3.and.NPSWIP.eq.NPSWIPprev) THEN
						    IF (IFAM.eq.3) THEN
!C
!C                 COMPUTATION OF SOLVENT LOSS AND SELECTIVITY WHITH
!C                 SINGLE VALENCY GROUPS
                                if(mop/=4)then !PROVISORIO para reacciones!!!!!
							    IST=ISF+1
								MS(NC,IST,2) = 1
								DO 41 I=1,MSV1 
								    MS(NC,IST,1) = NGSV1(I)
                                       
                                    call gc_special_groups(ms(3,:,:),bulky,nbulky)
                                    recorreSolutes => InputProblem01%MixtureInput%Solutes
                                    it=1
                                    nullify(ListPerformance)
			     		            do while (contar)
						                ms(1,:,:) = recorreSolutes%Formula
						                BPSolute = recorreSolutes%BoilingPoint
                                        call Evaluate_Mixture(.True.,FinalStruct,MS,NC,temp1,&
										BPSolute,ifam,NPSWIP,NPSBI,contar,ListPerformance,it)							                								    
                                        recorreSolutes => recorreSolutes%next
                                        if(.not.associated(recorreSolutes)) exit
            	                            it=it+1
                                    enddo                                             
                                       ! if(contar)call Incorporate_Solvent (mop,ms(3,:,:),ListProperties,ListPerformance,FMSs)

									MS(NC,IST,1) = 0
									MS(NC,IST,2) = 0
                                    if(contar)exit
41                              CONTINUE
                                endif
							END IF
						END IF
                    endif
!C  Continuo aqui
 117				IF (CONTAR) THEN
						JIS=JIS+1
						nIMS(JIS)=J
						Mprev(jis)=Mprv(k)
						if (sms) then
							ns=ns+1
							nsmsv(ns)=jis
						else if (sm1) then
							n1=n1+1
							nsm1v(n1)=jis
						end if
						DO 43 I=1,MDV
							NUD(I,JIS)=MUD(I,J)
  43						CONTINUE
						DO I=1,MSV
							NUD(MDV+I,JIS)=itermp(k,I)
						END DO
					END IF
					IF (JIS.EQ.NMAX) THEN
 900 					continue
						WRITE (6,880) NMAX,FAMILY(IFAM)
						READ (5,890) RESP
						IND1 = INDEX('CEce',RESP)
						IF (IND1.EQ.0) THEN
							GOTO 900
						ELSE IF ((RESP.EQ.'E').OR.(RESP.EQ.'e')) THEN
							SALIR=.TRUE.
							GOTO 1111
						ELSE
							GOTO 910
						END IF
					END IF
				END IF
  35		END DO
 205	CONTINUE
		PRINT*,CHAR(7)
		WRITE (6,1010)
		READ (5,890) SIGUE
!C
!c	 Quedan almacenadas en NUD, JIS pre-SMS para las cuales 
!c		NOACE  = false (todos los a necesarios est�n disponibles)
!c		SOLDES = false (factible)
!C		CONTAR = true (aceptable -en cuanto a SCLLI y SLSUP1 o SOLV)
!C	 y por lo tanto	siguen en carrera.
!C
!C        TERMINATION OF PROMISING pre-SMSs
!C
910		WRITE (6,2050)
!c			K = 0
!c			DO 510 I1=1,MSV
!c				DO 520 I2=I1,MSV
!c					K = K + 1
!c					NUS(I1,K) = 1
!c					NUS(I2,K)= NUS(I2,K) + 1
!c520				CONTINUE
!c510			CONTINUE
!c			KST = K
!c			IZ = 0
!c
!c	NUS = Vector que contiene todos los pares posibles entre los
!c		  MSV grupos terminales (se usa -cuando NATOT=2- para 
!c	      completar cada pre-SMS)				
!C ------------------------------------------------------------------
!C
!C        TERMINATION OF PROMISING STRUCTURES
!C
		JIST = 0
		NARES = 1
		IFAM1 = IALAR(1)
!c	Para incluir una "pre-SMS vacia" y permitir SMS de dos grupos:
		if (ifam.eq.3) then
			Npre = Npre + 1 
			NPSBI=NPSBI+1
			jis=jis+1 
			nIMS(JIS)= JEI + 1
			do 358 i=1,mdv+msv
				nud(i,jis)=0
358			continue
			Mpre=Mpre+1
			Mprev(jis)=Mpre
			K1maxpv(Mpre)=1
			Kpv(Mpre,1) =  1
			Kpv(Mpre,2) = -2
			Kpv(Mpre,3) =  0
		end if
!c
		NmFSG = 0
		FinalStruct=.True.
		NFSG = 0
		NSEW = 0
		NSBI = 0
		NSWIP = 0
		ns=1
		n1=1
601		DO 600 J=1,JIS
			ntcpa=ntcp
			DO 548 M=1,NSCM
				DO 548 K=1,2
					MS(NC,M,K) = 0
548			CONTINUE
			ISF = 0
			DO 605 I=1,MDV
				IF (NUD(I,J).NE.0) THEN
					ISF=ISF+1 !isf=n�mero de subgrupos distintos en la pre-SMS
					IMS(ISF,1) = NGDV(I)
					IMS(ISF,2) = NUD(I,J)
				END IF
605			CONTINUE
			isfi=isf
			DO 606 I=1,MSV
				IF (NUD(MDV+I,J).NE.0) THEN
					ISF=ISF+1 
					IMS(ISF,1) = NGSV(I)
					IMS(ISF,2) = NUD(MDV+I,J)
				END IF
606			CONTINUE
!C	El vector IMS se usa como intermediario entre NGDV-NGSV/NUD y MS(NC,_,_)
!C	para retener invariable la estructura intermedia, ya que MS es
!C	modificada por la subrutina INGRESAR_TERMINAL 
!C
			if (ifam.gt.3) then
				if (nsmsv(ns).eq.j) then
					ns=ns+1
					ikst=1
					ntcpa=0
					goto 118
				else if (nsm1v(n1).eq.j) then
					n1=n1+1
					ntcpa=1
				end if
			end if
			IF (IFAM.EQ.1) THEN !No deber�a pasar nunca por ac�...
				CALL AROMN (IPAREQ,ISF,NAC,NACCH,NACCH2,NATOT)
				IF (NATOT.EQ.0) THEN
					IKST = 1
				ELSE IF (NATOT.EQ.1) THEN
					IKST = MSV
				ELSE IF (NATOT.EQ.2) THEN
					IKST = KST
				END IF
			ELSE IF (IFAM.EQ.2) THEN ! ...ni por ac�.
				IKST = 1
				NATOT = 0
			ELSE ! ifam = 3, 4 � 5
				if(Mprev(j).ne.Mprev(j-1)) then	! para generar nuevas 
!c		terminaciones (meta y concretas) cuando cambia la meta pre-SMS
					ikst=0
					if (NGK1.eq.0) then
						K1max=0
					else
						K1max=K1maxpv(Mprev(j))
					end if
					if (NGK1J2+NGM1.eq.0) then
						K1min=K1maxpv(Mprev(j))
					else
						K1min=0
					end if
					do K1=K1min,K1max
						if (NGK1J2.eq.0) then
							K1J2max=0
						else
							K1J2max=ntcpa-K1
						end if
						if (NGM1.eq.0) then
							K1J2min=ntcpa-K1
						else
							K1J2min=0
						end if
						do K1J2=K1J2min,K1J2max
							M1=ntcpa-K1-K1J2
							if(Kpv(Mprev(j),1)+K1+K1J2.gt.0) then
								if(Kpv(Mprev(j),2)+2*K1+K1J2.gt.0) goto 907
							else
								if(Kpv(Mprev(j),3)+K1.gt.0) goto 907
							end if
							NmFSG=NmFSG+1 ! Number of metha 
!												Final Solvents Generated
!c Ahora se forman todas las combinaciones concretas, correspondientes a la actual
!c "meta-terminaci�n" para la preSMS concreta No. j, y v�lidas tambi�n para todas 
!c	  las otras IMS con el mismo Mv(j) (o sea, hijas de la misma meta-IMS)
							iK1=NCOMBIN(K1,NGK1)
							iK1J2=NCOMBIN(K1J2,NGK1J2)
							do i1=1,iK1
								do i2=1,iK1J2
									ikst=ikst+1
									i=0
									if (M1.ne.0) then
										i=1
										iterm(ikst,i,1)=NGSV(msv)
										iterm(ikst,i,2)=M1
									end if
									if (K1.ne.0) then
									 do L=1,NGK1
									  if (ncomb(K1,i1,L).ne.0) then
									    i=i+1
										iterm(ikst,i,1)=NGSV(L)
										iterm(ikst,i,2)=ncomb(K1,i1,L)
									  end if
									 end do
									end if
									if (K1J2.ne.0) then
									 do L=1,NGK1J2
									  if (ncomb(K1J2,i2,L).ne.0) then
									  i=i+1
									  iterm(ikst,i,1)=NGSV(NGK1+L)
									  iterm(ikst,i,2)=ncomb(K1J2,i2,L)
									  end if
									 end do
									end if
									niterm(ikst)=i
								end do
							end do
						end do
 907				end do
				end if
			END IF
!c	
!c			ISF = Numero de subgrupos diferentes en la pre-SMS
!c			IST = Numero de subgrupos diferentes en cada SMS
!c			
			if(nIMS(J).ne.nIMS(J-1)) then
				nIMSdist= nIMSdist+1
				n1SMS=NFSG+1
				nact=0
			end if
!c	n1SMS: orden de la 1er SMS correspondiente a la actual IMS
!c
!c	El siguiente ciclo es para formar IKST SMS a partir de la pre-SMS J
 118		DO 625 K=1,IKST
				NFSG=NFSG+1 ! Number of Final Solvents Generated
				rep = .FALSE.
				nact=nact+1
				SOLDES = .FALSE.
				IST = ISF
				DO 777 L=1,ISF
					MS(NC,L,1) = IMS(L,1)
					MS(NC,L,2) = IMS(L,2)
 777			CONTINUE
!c	la pre-SMS se copi� en MS(NC, , )
				if (ifam.ge.3) then
					if(nsmsv(ns-1).ne.j.and.nsm1v(n1-1).ne.j) then
						do L=1,niterm(k)
							do i=isfi+1,isf  ! isfi: No. de grupos dif. en la IMS
								if(IMS(i,1).eq.iterm(k,L,1)) then
									MS(NC,i,2)=MS(NC,i,2) + iterm(k,L,2)
									goto 778
								end if
							end do
							IST = IST + 1
							MS(NC,ist,1) = iterm(k,L,1)
							MS(NC,ist,2) = iterm(k,L,2)
 778					end do
						do i=1,isfi 
							TERM(nact,i,1) = 0
							TERM(nact,i,2) = 0
						end do	
						do i=isfi+1,ist 
							TERM(nact,i,1) = MS(NC,i,1)
							TERM(nact,i,2) = MS(NC,i,2)
						end do	
						TERM(nact,ist+1,1) = 0
						IF (nact.ne.1) then
							call eliminar_repetidos(nact,TERM,isfi,IST,rep)
							if(rep)then
								NFSG=NFSG-1
								nact=nact-1
								nSMSrep=nSMSrep+1
								GO TO 625
							end if
						END IF
					end if
				    IF (IFAM.EQ.4) THEN
					CALL AROMATICOS_ESPECIALES (IFAM1,IPAREQ,MS,IST)
				    END IF			
				else
				    IF (NATOT.EQ.0) THEN
					    IST = ISF
				    ELSE IF (NATOT.EQ.1) THEN
						CALL INGRESAR_TERMINAL (MS,ISF,NGSV(K),1)
						IST = ISF + 1

						IF (IFAM.EQ.1) CALL AROM1 (IPAREQ,IST,NACCH,NACCH2,SOLDES)

					ELSE IF (NATOT.EQ.2) THEN
						NSV = 0
						DO 49 I=1,MSV
							IF (NUS(I,K).NE.0) THEN
								NSV=NSV+1
								NSV1(NSV,K)=NGSV(I)	! n�mero de subgrupo
								NSV2(NSV,K)=NUS(I,K)  ! cantidad (1 o 2)
								NGSSV1(NSV,1)=NGSSV(I,1) !propiedades combinat.
								NGSSV1(NSV,3)=NGSSV(I,3)
							END	IF
  49					CONTINUE
						IF (NSV.NE.2) THEN ! cuando un mismo grupo esta 2 veces
							CALL INGRESAR_TERMINAL (MS,ISF,NSV1(1,K),2)
							MS(NC,ISF+2,1)=0
							MS(NC,ISF+2,2)=0
							IST=ISF+1
						ELSE  ! para dos grupos terminales diferentes
							CALL INGRESAR_TERMINAL (MS,ISF,NSV1(1,K),1)
							CALL INGRESAR_TERMINAL (MS,ISF+1,NSV1(2,K),1)
							IST=ISF+2
						END IF
!c	Ahora qued� la SMS completa, almacenada "provisoriamente" en MS(NC,_,_)
!C
!C					FEASIBILITY TEST OF FINAL STRUCTURES
						IF (IFAM.EQ.1) THEN
							CALL AROM1 (IPAREQ,IST,NACCH,NACCH2,SOLDES)
							CALL AROM2 (K,NACCH,NSV,NSV2,NGSSV1,SOLDES)
						ELSE IF (IFAM.EQ.3) THEN
							CALL ALIF2 (IPAREQ,SOLDES)
							IF (.NOT.(SOLDES)) THEN
								MSST=0
								KSST=0
								DO 52 I=1,NSV
									MSST=NSV2(I,K)*NGSSV1(I,1)+MSST
									KSST=NSV2(I,K)*NGSSV1(I,3)+KSST
  52							CONTINUE
								MS1=MSST+MTEST
								KS1=KSST+KTEST
								MJT=MS1+JTEST/2.0

								IF (KS1.GT.MJT) SOLDES = .TRUE.

							END IF
						END IF
					END IF
				end if
!					CALL MOSTRAR_COMPONENTE (IST)
				IF (.NOT.(SOLDES)) THEN
					DO L=IST+1,10
							MS(NC,L,1) = 0
							MS(NC,L,2) = 0
					END DO
                    call gc_special_groups(ms(3,:,:),bulky,nbulky)
                        
                    bulkyloop: &
                    do i = 1, nbulky 
                        if(mop==4)then !PROVISORIO para reacciones !!!!!!!
							call Evaluate_Pure (FinalStruct,ms(3,:,:),temp1,ifam,NSEW,solace,ListProperties)   
                        else
                        call Evaluate_Pure (FinalStruct,ms(3,:,:),temp1,ifam,NSEW,solace,ListProperties)   
                        recorreSolutes => InputProblem01%MixtureInput%Solutes
                        it=1
                        nullify(ListPerformance)
			     	    do while (solace)
					        ms(1,:,:) = recorreSolutes%Formula
					        BPSolute = recorreSolutes%BoilingPoint
                            call Evaluate_Mixture(.True.,FinalStruct,MS,NC,temp1,&
							BPSolute,ifam,NSWIP,NSBI,solace,ListPerformance,it)						    					
                            recorreSolutes => recorreSolutes%next
                            if(.not.associated(recorreSolutes)) exit
                            it=it+1
                        enddo     
                        endif
                       ! call llecalas(InputProblem01%T,InputProblem01%P,ms,3,xtemporal)

                        if(solace)then
                            call Incorporate_Solvent (mop,ms(3,:,:),ifam,ListProperties,ListPerformance)                        
                            jist = jist + 1
                        endif
                       ! nullify(ListPerformance)
                        
					    IF (JIST.GE.NMAX) THEN
 940						WRITE (6,930) NMAX,FAMILY(6)
						    READ (5,890) RESP
						    IND1 = INDEX('CEce',RESP)
						    IF (IND1.EQ.0) THEN
							    GOTO 940
						    ELSE IF ((RESP.EQ.'E').OR.(RESP.EQ.'e')) THEN
							    SALIR=.TRUE.
							    GOTO 1111
						    ELSE
							    GOTO 950
						    END IF
					    END IF
                    enddo bulkyloop
				END IF
 625		CONTINUE
600		CONTINUE
		IF (IFAM.EQ.4) THEN
			NARES = NARES + 1
			IF (IALAR(NARES).NE.0) THEN
				IFAM1 = IALAR(NARES)
				GOTO 601
			ELSE
				GOTO 950
			END IF
		END IF
	end if !ifam
!C
!C       SELECTION AND PRINTING OF MOST PROMISING SOLVENT STRUCTURES
950	CONTINUE
	PRINT*,CHAR(7)
	WRITE (6,1010)
	READ (5,890) SIGUE
	WRITE (2,955)
	WRITE (6,965)
	WRITE (2,965)
	IF (IFAM.NE.2) THEN
		IF (IFAM.NE.1) THEN
			WRITE (6,830) MEI
			WRITE (2,830) MEI
		END IF
		WRITE (6,820) JEI
		WRITE (2,820) JEI
		IF (IFAM.NE.1) THEN
		    WRITE (6,806) Mpre
		    WRITE (2,806) Mpre
		    WRITE (6,808) Npre
		    WRITE (2,808) Npre
		    WRITE (6,*) ' Number of pre-Solvents with high Mol. W.:',NPSEW
		    WRITE (2,*) ' Number of pre-Solvents with high Mol. W.:',NPSEW
		    WRITE (6,*) ' Number of pre-Solvents without int. par.:',NPSWIP
		    WRITE (2,*) ' Number of pre-Solvents without int. par.:',NPSWIP
		    WRITE (6,*) ' Number of pre-Solvents with bynary inf. :',NPSBI
		    WRITE (2,*) ' Number of pre-Solvents with bynary inf. :',NPSBI
		END IF
		WRITE (6,810) JIS
		WRITE (2,810) JIS
		IF (IFAM.NE.1) THEN
			WRITE (6,788) NmFSG
			WRITE (2,788) NmFSG
!c En realidad NmFSG est� sobredimensionado porque hay repetidos
!c Lo importante es que finalmente no haya solventes repetidos
!c y eso se solucion� con la subrutina eliminar_repetidos.
!c Entonces NFSG efectivamente indica un n�mero de solventes distintos
!c Debido a esto, siempre que no haya mas de un grupo terminal 
!c por metha grupo, NFSG ser� mayor que NmFSG.
		END IF
	END IF
	WRITE (6,789) NFSG
	WRITE (2,789) NFSG
	WRITE (6,*) ' Number of Solvents with high mol. weight:',NSEW
	WRITE (2,*) ' Number of Solvents with high mol. weight:',NSEW
	WRITE (6,*) ' Number of Solvents without int. par.    :',NSWIP
	WRITE (2,*) ' Number of Solvents without int. par.    :',NSWIP
	WRITE (6,*) ' Number of Solvents with bynary inf.     :',NSBI
	WRITE (2,*) ' Number of Solvents with bynary inf.     :',NSBI
	WRITE (6,790) JIST
	WRITE (2,790) JIST
	WRITE (6,1010)
	READ (5,890) SIGUE

	IF (KILOUT.EQ.1) WRITE (6,3290)	

1111 CONTINUE
!     FORMATS
788   FORMAT(' ',/,' NUMBER OF metha-FINAL SOLVENTS GENERATED:        ',I10)
789   FORMAT(' ',/,' NUMBER OF FINAL SOLVENTS GENERATED:              ',I10)
790   FORMAT(' ',/,' NUMBER OF FINAL SOLVENTS ACCEPTED:               ',I10)
806   FORMAT(' ',/,' NUMBER OF metha-pre-FINAL SOLVENTS GENERATED:    ',I10)
808   FORMAT(' ',/,' NUMBER OF pre-FINAL SOLVENTS GENERATED:          ',I10)
810   FORMAT(' ',/,' NUMBER OF ACCEPTED pre-FINAL SOLVENTS:           ',I10)
820   FORMAT (' NUMBER OF INTERMEDIATE STRUCTURES GENERATED:      ',I10)
830   FORMAT (' NUMBER OF METHA-INTERMEDIATE STRUCTURES GENERATED:',I10)
880   FORMAT (' ','*** Number of pre-final structures ',&
                'for termination too big ***',&
                //,10X,'Options:',//,19X,'Continue (C): Only ',&
                'the first ',i5,' pre-final structures',/,&
                33X,'will be taken into account.',//,19X,'Exit     ',&
                '(E): Choose fewer ',a21,'.',//,73x,'> ',$)
890   FORMAT (A)
930   FORMAT (' ','*** Number of final structures for termination ',&
                'too big ***',&
                //,10X,'Options:',//,19X,'Continue (C): Only ',&
                'the    first   ',i5,'    final   structures',/,&
                33X,'will be taken into account.',//,19X,'Exit     ',&
                '(E): Choose fewer ',a21,'.',//,76x,&
                '> ',$)
 955  FORMAT (' ',////)
 965  FORMAT (' ','*** MOLECULAR DESIGN RESULTS ***',//)
1010  FORMAT (' ',//,70X,'<RET>',$)
!2040  FORMAT (' ','***** GENERATING INTERMEDIATE STRUCTURES...')
2050  FORMAT (' ',/,' ***** GENERATING FINAL STRUCTURES...')
3290  FORMAT (' ','***** CREATING OUTPUT FILE...')
!4010  FORMAT (' ',/,' ***** EVALUATING SINGLE GROUP SOLVENTS...')     
      RETURN
endsubroutine structure_generator 

!================================================================================
subroutine Incorporate_Solvent (mop,compound,tsolv,ListProperties,ListPerformance)
    use StructuresDesign
    implicit none
!Variables de ENTRADA
    integer,intent(in)::mop,compound(DiffStructGroups,2),tsolv
    type(CalculatedProperties),pointer,intent(in)::ListProperties
    type(CalculatedPerformance),pointer,intent(in)::ListPerformance
!Variables de ENTRADA/SALIDA
   ! type(FinalStructure),pointer,intent(inout)::FMSs        
    if(mop==0)then
    elseif(mop==1)then
        call Incorporate_Structure (compound,tsolv,FMSs,ListProperties,ListPerformance)
        return
    elseif(mop==2)then
        call Incorporate_Structure (compound,tsolv,FMSs,ListProperties,ListPerformance)
        return
    elseif(mop==3)then
!	    FMSs%Selectivity = SelT
!	    FMSs%GroupsNumber = NMGRUP
!	    FMSs%Kow = KowT
!	    FMSs%FunctionalGroups = FuncGroupsT
	    return    
    elseif(mop==4)then !PROVISORIO para reacciones!!!
        call Incorporate_Structure (compound,tsolv,FMSs,ListProperties,ListPerformance)
        return        
    endif   
endsubroutine Incorporate_Solvent   

!=====================================================================
subroutine gc_special_groups(compuesto,bulky,nbulky)

use CONSTANTES
implicit none
integer,parameter::NESTM=7100
!Variables ENTRADA/SALIDA
integer,dimension(DiffStructGroups,2),intent(inout)::compuesto
integer,dimension(NESTM,DiffStructGroups,2),intent(inout)::bulky
integer,intent(out)::nbulky
!Variables INTERNAS
integer::i
integer::nC,nCH,nCH2,nCH3
logical::bulk


!
!bulk = .False.
!i=0
!do while(compuesto(i+1,1) /= 0)
!    i=i+1
!    if(compuesto(i,1) == 60) then
!        bulk = .True.
!        nC = compuesto(i,2)
!    endif
!enddo
!
!if(.not. bulk) return   !endsubroutine
!
!i=0
!do while(compuesto(i+1,1) /= 0)
!    i=i+1
!    if(compuesto(i,1) == 1)then 
!        compuesto(i,1) = 4
!        cycle
!    endif
!    if(compuesto(i,1) == 2)then 
!        compuesto(i,1) = 5
!        cycle
!    endif
!    if(compuesto(i,1) == 3)then 
!        compuesto(i,1) = 6
!        cycle
!    endif
!enddo
nbulky = 1
bulk = .False.
i=0
do while(compuesto(i+1,1) /= 0)
    i=i+1
    if(compuesto(i,1) == 60) bulk = .True.
enddo

if(.not. bulk) return   !endsubroutine

i=0
nbulky = 0 
do while(compuesto(i+1,1) /= 0)
    i=i+1
    if(compuesto(i,1) == 1)then 
        compuesto(i,1) = 4
        nbulky = nbulky + compuesto(i,2)
        if (nbulky==4)goto 900
        cycle
    endif
    if(compuesto(i,1) == 2)then 
        compuesto(i,1) = 5
        nbulky = nbulky + compuesto(i,2)
        if (nbulky==4)goto 900
        cycle
    endif
    if(compuesto(i,1) == 3)then 
        compuesto(i,1) = 6
        nbulky = nbulky + compuesto(i,2)
        if (nbulky==4)goto 900        
        cycle
    endif
enddo

900 nbulky = 1

endsubroutine gc_special_groups
   endmodule
   
