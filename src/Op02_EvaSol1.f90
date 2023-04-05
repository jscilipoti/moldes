SUBROUTINE evasol1 (nfunc)
!    use PropertiesData
!    use PureProp
!    use PROP_PDB, only:checkdatos,elegir_CXR,pvomegaybp
!    use input_data,only:ingresar_solventes, enter_problem
!    use SubGrupos
!    use StructuresDesign
!    use CONSTANTES
!    
!    PARAMETER (NMOD=3,NMSV=20,NA=150,NSVA=20,NSCM=10,NESTM=7100,IS1=5,NAA=4,ng=70)
!    implicit real*8 (A-H,O-Z)
!    EXTERNAL Leer_In,NCOMBIN
!    integer,intent(in)::nfunc
!    
!    
!    type(FinalStructure),pointer::FMSs
!    type(CalculatedProperties),pointer::ListProperties
!    type(CalculatedPerformance),pointer::ListPerformance
!    LOGICAL CONTAR,SOLACE,SOLDES,NOACE,ARCH,LICH2,LICH
!    LOGICAL USER,rep,sms,sm1,entro2,entro3
!	CHARACTER*8 FS(NA)
!    CHARACTER*70 TITULO
!    CHARACTER*70 TITL
!    CHARACTER*1 RESP,SIGUE,OPT
!    CHARACTER*35 NOMB, IFILE
!    CHARACTER*21 FAMILY(6)
!	integer UNIF(20),grupo1(na),maingr(na),Mprev(NMAX),Mprv(NMAX)
! 	CHARACTER*35 ANAME,FORM2
!	CHARACTER*20 FORM1
!	integer comp(20),ident(8),identsolv(NMAX,8),nisomsolv(NMAX)
!	character*35 nombre(8),formula(8) 
!	character*35 nombresolv(NMAX,8),formulasolv(NMAX,8) 
!	INTEGER NUD(NSVA,NMAX),NFAv(NESTM),K1maxv(NESTM)
!    integer K1maxpv(NMAX),ncomb(10,1001,5),MUD(NSVA,NMAX)
!    integer itermp(nmax,nscm),TERM(NMAX,NSCM,2)
!    real*8 TB(NMAX,8),DENST(NMAX,8),VLIQ(NMAX,8)
!    real*8 TMP
!    real*8 RDER(NMAX)
!    real*8 VREL(NMAX),VRELSOLV(NMAX,8)
!	dimension isv(NMG),isd(NMG),idv(NMG),icyc(NMG)
!    dimension iar(NMG),igrp(NMG),idrs(NMG,NSEL),nyrs(NMG,NSEL)
!    dimension grupos(nsel),kgrup(nsel),i3v(NMG),i4v(NMG)
!	dimension idr(nsel),nyr(nsel),idf(nsel),nyf(nsel)
!    DIMENSION ISOL(NMAX,2),NSV1(3,NMAX),NSV2(3,NMAX)
!    DIMENSION IGR0(NMAX,2),IGR1(NMAX,2),IGR2(NMAX,2),IALAR(NAA)
!    DIMENSION IGR3(NMAX,2)
!    dimension NPINT(NG),NINTT(NG)
!    DIMENSION NGDV(0:NA),NGSV(NSVA),NGSV1(NSVA),NGSDV(NSVA,5)
!    DIMENSION NGSSV(NSVA,5),NGSDV1(NSVA,5),NGSSV1(NSVA,5)
!    DIMENSION MS(NCOM,NSCM,2),IMS(NSCM,2)
!    DIMENSION NUS(NSVA,NMAX),MST(NMAX,NSCM,2),nicomp(nmax)
!    DIMENSION NGM1V(NMSV),nmst(nmax)
!    DIMENSION NPUNT(NA),NGRUP(NA),niterm(nmax)
!    DIMENSION MM(NMG),MJ(NMG),MK(NMG),MI(NMG),MH(NMG),PMG(NMG)
!    DIMENSION DELTC(NMG),DELV(NMG),DELPC(NMG),delvi(NMG),deln(NMG)
!    DIMENSION ICH2(NMOD),ICH(NMOD),IACOH(NMOD),IACCH2(NMOD)
!    DIMENSION IACCH(NMOD),IAC(NMOD),IACH(NMOD),IACCL(NMOD)
!    DIMENSION IACNH2(NMOD),IACNO2(NMOD),ICH3(NMOD),IOH(NMOD)
!    DIMENSION ICOOH(NMOD),Mv(NMAX),Kv(NESTM,3),iterm(nmax,10,2)
!    DIMENSION Kpv(NMAX,3),nIMS(NMAX),nsmsv(nmax),nsm1v(nmax)
!    common/estintram/MUD,ncomb,NFAv,K1maxv,Mv,Kv
!    common/int/ifin2,kgrup,maingr
!	COMMON/gru/GRUPO1,NGRMAX
!	COMMON/NOM/FS
!    COMMON/PUNSUB/NPUNT,NGRUP,NCANT
!    COMMON/PUNGRU/NPINT,NINTT,NUMINT
!    COMMON/US/MS
!    COMMON/IUS/IMS
!	COMMON/CONT/PMG,DELTC,DELV,DELPC,DELVI,DELN
!    COMMON/GRUPAR/IACH,IACCH2,IACCL,IACNH2,IACNO2
!    COMMON/GRUPAL/ICH3,ICH2
!    COMMON/GRUPES/ICH,IACOH,IACCH,IAC,IOH,ICOOH
!    COMMON/STAT/MM,MJ,MK,MI,MH
!    COMMON/EXTDIS/PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
!    COMMON/MOLCOM/ARCH
!	COMMON/SOLVISOM/NISOMSOLV,IDENTSOLV,NOMBRESOLV,FORMULASOLV,TB,DENST,VRELSOLV
!    COMMON/EVAL/USER
!    COMMON/PAREQ/IPAREQ
!    data family /'      aromatic groups','     molecular groups',&
!                    3*'  dual valence groups','single valence groups'/
!!...SENTENCIAS
!    idev = 6
!    idevr= 5
!!...LECTURA DE SOLUTOS
!    call enter_problem(nfunc)
!!    nopt = 0
!!    do while (nopt<1.or.nopt>2)
!!        write (idev,"(////' ','Component read from:',&
!!                              //,11x,' keyboard:           :  1',&
!!                               /,11x,' input file          :  2',&
!!     	                      //,60x,'> ',$)")
!!        read (idevr,*) nopt
!!    enddo
!!    if (nopt==1) then !por teclado
!!        write (idev,"(///' ','How many component to be recovered do you want evaluate?',//,60x,'> ',$)") 
!!        read(idevr,*) NSOLUT
!!    else !por archivo
!! 112    write (idev,*) "Give the file name (maximun 35 characters):"
!!        read (idevr,"(a)",err=112) IFILE
!!        open(unit=20,file=IFILE,status='old',access='direct', form='formatted',RECL=62)
!!        Do j=1, 100
!!            read  (20,"(20i3)",rec=j,err=113) (idrs(j,i),nyrs(j,i),i=1,10)
!!            write (6,"(20i3)") (idrs(j,i),nyrs(j,i),i=1,10)
!!        enddo
!! 113    NSOLUT = j-1   
!!    endif
!    MS(:,:,:) = 0    
!    Do 100 S=1, NSOLUT 
!        MS(1,:,:) = 0
!        MS(3,:,:) = 0
!        IS=5
!	    if (S.EQ.1) THEN 
!	        call ingresar_solutos(MOP,tem,entro2,entro3,&
!     						idrs,nyrs,msol,idf,nyf,mraf,bp1,a1,a2,a3,&
!     						dens2,b1,b2,b3,tazeo,xazeo,nsolut,nopt)
!	        do j=1,mraf !Carga el CPR
!	            MS(2,J,1)=idf(j)
!		        MS(2,J,2)=nyf(j)
!		        call cr_puntf(idf(j),NPUNT,NGRUP,NPINT,NINTT)
!	        end do     						
!	    ENDIF
!	    msol = 0
!	    do j=1,10 !Carga el CAR
!	        if (idrs(s,j).ne.0) then
!	            msol=msol+1
!	            MS(1,J,1)=idrs(s,j)
!		        MS(1,J,2)=nyrs(s,j)
!		        call cr_puntf(idrs(s,j),NPUNT,NGRUP,NPINT,NINTT)
!		    endif
!	    end do
!!c	MSOL = Cantidad de subgrupos diferentes en el CAR
!!c	MRAF = Cantidad de subgrupos diferentes en el CPR
!
!!c	USER indica si molde2eco fue llamada para evaluar solventes 
!!c	a ingresar por el usuario (true), o aquellos que se generarán 
!!c	a partir de un conjunto de grupos preseleccionados (false).
!
!	    kilout=0
!	    call car_car (ipareq,isv,isd,idv,i3v,i4v,icyc,iar,igrp,&
!                   imasv,imasd,imadv,n3v,n4v,imacy,imaar,imagr)
!	    if(ngrmax.eq.0)then
!	        ngrmax = Size_LSubGroups()
!		    do 19 i=1,ngrmax
!		        grupo1(i) = i
!			    maingr(i) = mainsg (i,ipareq) !función q devuelve n° de main group
!  19		continue
!	    end if
!	    if (s.eq.1) call ingresar_solventes (ipareq,jist,mst,nmst,nicomp)
!	    mdv=0
!	    do 147 i=1,jist
!		    do 148 j=1,nmst(i)
!			    mdv=mdv+1
!			    ngdv(mdv)=mst(i,j,1)
!  148		continue
!  147	continue
!	    CALL Store_Pr (IPAREQ)
!        CALL Store_In (IPAREQ)
!	    WRITE (6,1010)
!	    READ (5,890) SIGUE
!	    OPEN (UNIT=2,FILE='SOLVENTS.MDO',FORM='FORMATTED')
!!C
!!C         RAFFINATE AND SOLUTE PROPERTIES
!        PM2 = MolecularWeight (MS(2,:,:))
!        PM1 = MolecularWeight (MS(1,:,:))
!!c	PM1 and PM2 are the molecular weights of the CAR and CPR
!        IF (MOP.EQ.2) THEN !Destilación extractiva
!		    IF (TAZEO.EQ.1.0) THEN
!			    TEMP1 = BP1		!CAR normal boiling point
!		    ELSE
!			    TEMP1 = TAZEO	!normal azeotropic temperature
!		    END IF
!		    PSAT1 = EXP(A1 -A2/(TEMP1 + A3))
!		    PSAT2 = EXP(B1 -B2/(TEMP1 + B3))
!        ELSE ! Extracción líquido-líquido
!		    TEMP1 = tem   !operation temperature
!        END IF
!
!!C     PRINTING OF GROUP PARAMETERS FOR THE KEYS
!        MODEL=0			   ! ?
!       ! CALL ESCLIS (2,MOP,IPAREQ,MSOL,MRAF,NGDV,NGSDV,MDV,&
!     !		     NGSV,NGSSV,MSV,tem,BP1,TITL,A1,A2,&
!     !             A3,B1,B2,B3,TAZEO,X1AZEO,IFAM)
!        IF (KILOUT.EQ.0) CALL ESCLIS (6,MOP,IPAREQ,MSOL,MRAF,NGDV,NGSDV,MDV,&
!     	     		 NGSV,NGSSV,MSV,tem,BP1,TITL,A1,A2,&
!                      A3,B1,B2,B3,TAZEO,X1AZEO,IFAM)
!        LICH2 = .FALSE.
!        LICH = .FALSE.
!        CALL BUSCARAS (ICH2(IPAREQ),NGRUP,NMG,LICH2)
!        CALL BUSCARAS (ICH(IPAREQ),NGRUP,NMG,LICH)
!        if (LICH2) then
!	        TCCH2=DELTC(NPUNT(ICH2(IPAREQ)))
!	        PCCH2=DELPC(NPUNT(ICH2(IPAREQ)))
!	    endif
!	    if (LICH) PCCH= DELPC(NPUNT(ICH(IPAREQ)))
!        
!        nullify(FMSs)
!	    DO 200 J=1,JIST !bucle donde se evalúan los solventes
!	        if (nicomp(j).eq.3) then
!		        if (LICH2) DELTC(NPUNT(ICH2(IPAREQ)))=.013
!			    if (LICH2) DELPC(NPUNT(ICH2(IPAREQ)))=.184
!			    if (LICH) DELPC(NPUNT(ICH(IPAREQ)))=.192
!		    else
!			    if (LICH2) DELTC(NPUNT(ICH2(IPAREQ)))=TCCH2
!			    if (LICH2) DELPC(NPUNT(ICH2(IPAREQ)))=PCCH2
!			    if (LICH) DELPC(NPUNT(ICH(IPAREQ)))= PCCH
!		    end if
!		    DO 702 L=1,NMST(J)
!			    MS(3,L,1) = MST(J,L,1)
!			    MS(3,L,2) = MST(J,L,2)
!			    call cr_puntf(MST(J,L,1),NPUNT,NGRUP,NPINT,NINTT)
! 702		CONTINUE
!		    DO L=NMST(J)+1,10
!			    MS(3,L,1) = 0
!			    MS(3,L,2) = 0
!		    END DO
!
!            call Calc_Prop_Pure(ms(3,:,:),nicomp(j),temp1,Bpt=bpoint,HVAPT=Hvb,VISCT=visl,FuncGroup = nFuncGroupsT)
!            call Solvent_Properties (MS,IPAREQ,MOP,TEMP1,Selty,spower,sperd,distri,RELV,NOACE)
!            call Incorporate_Structure (MS(3,:,:),FMSs,ListProperties,ListPerformance)
!            FMSs%MolecularWeight = MolecularWeight (MS(3,:,:))
!            FMSs%BoilingPoint = BPOINT
!            FMSs%Hvap = HVB
!            FMSs%DifferenceBP = bpoint - BP1
!	        FMSs%Viscosity = VISL
!!		    FMSs%Selectivity = SelT
!!		    FMSs%SolventPower = SolPowT*100
!!		    FMSs%GroupsNumber = NMGRUP
!!		    FMSs%SolventLost = SlostT*100
!!		    FMSs%DistCoefficient = DistCoefT
!		    FMSs%FunctionalGroups = nFuncGroupsT
!            FMSs%RelVolatility =  RELV
!     	    IF (MOP.EQ.2) THEN
!			    SLV1 = SPOWER*PMOL/PM1
!			    CALL MI_SOLV_BREAK (IPAREQ,TEMP1,SLV1,SOLACE,SPERD)
!		    ELSE
!			    FMSs%RDER = DENSID/DENS2  
!		    END IF
!		    FMSs%Selectivity = SELTY
!            FMSs%SolventPower = SPOWER*100
!            FMSs%GroupsNumber = groups_number(MS(3,:,:))
!            FMSs%SolventLost = SPERD*100
!            IF (MOP.EQ.1) THEN
!		        FMSs%DistCoefficient = DISTRI
!		    ELSE
!			    FMSs%DistCoefficient = 100*Sel(J)/(SLost(J)*MW(J))
!            END IF
!	!	    DC2(J) = DistCoef(J)
!200	    CONTINUE
!        N=0
!	    if (ipareq.eq.2) then
!		    write(6,*) ' Checking solvent structures with data base'
!		    nsol=100
!            call score(FMSs,nsol)
!	    end if
!        call write_results (FMSs,mop,6,ifam,nsol)
!        call write_results (FMSs,mop,2,ifam,nsol)
!        WRITE (6,1010)
!        READ (5,890) SIGUE
! !  call ci_ban1
!100 continue
!1111 CLOSE (UNIT=1)
!	CLOSE (UNIT=2)
!!C	CLOSE (UNIT=4)
!2020IF (IERROR.NE.0) THEN
!        WRITE (6,*) '*****',nomb(1:lenfile)//'.mdi',': file not found *****'
!        WRITE (6,1010)
!        READ (5,890) SIGUE
!        ARCH = .FALSE.
!    END IF
!
!!C     FORMAToS
!955  FORMAT (' ',////)
! 965  FORMAT (' ','*** MOLECULAR DESIGN RESULTS ***',//)
!3290  FORMAT (' ','***** CREATING OUTPUT FILE...')
!2040  FORMAT (' ','***** GENERATING INTERMEDIATE STRUCTURES...')
!2050  FORMAT (' ',/,' ***** GENERATING FINAL SOLVENTS...')
!4010  FORMAT (' ',/,' ***** EVALUATING SINGLE GROUP SOLVENTS...')
!!2030  FORMAT (' ','*****',A6,'.MDI: FILE NOT FOUND *****')
!1010  FORMAT (' ',//,70X,'<RET>',$)
!11    FORMAT (A70)
!10    FORMAT (5I3)
!15    FORMAT (20I3)
!16    FORMAT (11I3)
!18    FORMAT (3I3)
!40    FORMAT (' ',/,10I3)
!46    FORMAT (' ','SELECTIVITY',1F10.4,/,' SOLV.POWER',1F10.4,/)
!694   FORMAT (' ',//,70X,'>',$)
!830   FORMAT (' NUMBER OF METHA-INTERMEDIATE STRUCTURES GENERATED:',I10)
!820   FORMAT (' NUMBER OF INTERMEDIATE STRUCTURES GENERATED:      ',I10)
!806   FORMAT(' ',/,' NUMBER OF metha-pre-FINAL SOLVENTS GENERATED:    ',I10)
!808   FORMAT(' ',/,' NUMBER OF pre-FINAL SOLVENTS GENERATED:          ',I10)
!810   FORMAT(' ',/,' NUMBER OF ACCEPTED pre-FINAL SOLVENTS:           ',I10)
!788   FORMAT(' ',/,' NUMBER OF metha-FINAL SOLVENTS GENERATED:        ',I10)
!789   FORMAT(' ',/,' NUMBER OF FINAL SOLVENTS GENERATED:              ',I10)
!790   FORMAT(' ',/,' NUMBER OF FINAL SOLVENTS ACCEPTED:               ',I10)
!880   FORMAT (' ','*** Number of pre-final structures ',&
!                'for termination too big ***',&
!                //,10X,'Options:',//,19X,'Continue (C): Only ',&
!                'the first ',i5,' pre-final structures',/,&
!                33X,'will be taken into account.',//,19X,'Exit     ',&
!                '(E): Choose fewer ',a21,'.',//,73x,&
!                '> ',$)
!930   FORMAT (' ','*** Number of final structures for termination ',&
!                'too big ***',&
!                //,10X,'Options:',//,19X,'Continue (C): Only ',&
!                'the    first   ',i5,'    final   structures',/,&
!                33X,'will be taken into account.',//,19X,'Exit     ',&
!                '(E): Choose fewer ',a21,'.',//,76x,&
!                '> ',$)
!890   FORMAT (A)
!	RETURN
      END

!====================================================================
subroutine ingresar_solutos(MOP,tem,entro2,entro3,&
     						idrs,nyrs,msol,idf,nyf,mraf,bp1,a1,a2,a3,&
     						dens2,b1,b2,b3,tazeo,xazeo,nsolut,nopt)
   ! use PROP_PDB, only:checkdatos,elegir_CXR,pvomegaybp,antoine,check_inter_param
    use input_data,only:Leer_Comp,ingresar_solventes,checkdatos,elegir_CXR,pvomegaybp,antoine,check_inter_param
    use CONSTANTES
    use SubGrupos
    use Input
    implicit real*8 (a-h,o-z)
    external max_sub,mainsg
    parameter (ins=20)
	dimension idr(nsel),nyr(nsel),idf(nsel),nyf(nsel),kgrup(nsel)
	dimension pmg(NMG),iar(NMG),&
               igrp(NMG),idv(NMG),icyc(NMG),isv(NMG),isd(NMG),&
     	i3v(NMG),i4v(NMG),idrs(NMG,NSEL),nyrs(NMG,NSEL)
    integer comp(10,2),grupo1(NMG),ident(8),maingr(NMG),grupos(nsel),h1
    character*8 fs(NMG)
	character*35 nombre(8),formula(8)
    character*37 nomco1,nomco2,nomco3,nomco4
    character*1 opt,iopt
    logical prob,entro2,entro3

	COMMON/NOM/FS
	COMMON/nr/nrcpr
    common/int/ifin2,kgrup,maingr
    COMMON/gru/GRUPO1,NGRMAX
    
    idev = 6
	idevr = 5
!c---Apertura del banco de datos BADAUN
   ! call ab_ban1 							 

!C---Seleccion de la operacion de separacion
!     call limp (idev)
60    write (idev,20)  ! Choose the separation operation...
      write (idev,30) !"> "
      read (idevr,*,err=60) mop                                        
      if ((mop.lt.1).or.(mop.gt.2)) go to 60
!
!
!c
!c	Ingreso modificaciones para incluir opciones 3 y 4 en ipareq (23/5)
!c
!c---Seleccion de Parametros Unifac
!c     if (mop.eq.1) then
      CALL TABLA_PARAMETROS(IPAREQ)
!c     else
!c        ipareq = 2
!c     end if
!C
!C---LECTURA DE LA TEMPERATURA DE OPERACION (SI ELL)
!C
      IF (MOP.EQ.1) THEN !Liquid-Liquid extraction
		WRITE(IDEV,5555)! Give operation temperature.
		READ(IDEVR,*) TEM
      ELSE IF (mop.eq.2) THEN !Extractive destilation
		tem=1.0
      END IF
      call SubGroups_Characterisation()

!c---Lectura de los datos segun la seleccion de parametros
      write (idev,240) !*** Loading data ***
      ngrmax = Size_LSubGroups()
      do 10 i=1,ngrmax
         grupo1(i) = i
         maingr(i) = mainsg (i,ipareq)
  10  continue
      call car_car (ipareq,isv,isd,idv,i3v,i4v,icyc,iar,igrp,&
                   imasv,imasd,imadv,ima3v,ima4v,imacy,imaar,imagr)

!c---Lectura del componente a recuperar (CAR)
      Do 113 S=1, NSOLUT
        entro2 = .true.
        entro3 = .true.
 960    if (mop.eq.1) then !OPERACIÓN DE SEPARACIÓN
           nomco1 = 'Component to be recovered:           '
        else
           nomco1 = 'Less  volatile  component:           '
        end if
        if (nopt.eq.2) then
            msol=0
            do j=1,10 !Carga el CAR
	          if (idrs(s,j).ne.0) then
	              msol=msol+1
	      	  endif
	        enddo
            goto 112
        endif
        write (idev,540) !Give the group....
        call leer_comp (ipareq,nomco1,comp)
        !idr: vector con los num ident de los subgrup seleccionados
        !nyr: vector con las cantidades de cada subgrup seleccionado
        do i=1, 10
         idrs(s,i)=idr(i)
        enddo
        do i=1, 10
         nyrs(s,i)=nyr(i)
        enddo
 112    continue
!c---Iniciacion del vector de grupos principales (VGP) para control
!c---de parametros de interaccion grupal (a(nm);a(mn))
 330    ifin0 = 0
      call asignar_grupos_principales (idrs(s,:),msol,grupos) !sale así y
!									(maingr,subgrup,nsub,grupos) llega de esta manera
        call armar_grupos_distintos (grupos,msol,kgrup,ifin0,ifin1) !sale así y
!								(grupos,ngrup,kgrup,ifin0,ifin) llega de esta manera
!c---Control de parametros de interaccion para el CAR
!c     call limp (idev)         
        if (msol.gt.1) then !msol es el número de grupos principales
         write (idev,1910) !*** Checking interaction parameters ***
       if(model /= 3)then  
         if (check_inter_param (ipareq,1,kgrup,1)) then
            write (idev,1000)
            write (idev,1010)
            read (idevr,510) return
            go to 960
         end if
        end if
       endif
!c---Confirmacion del CAR
 3000   write (idev,1210) nomco1,(Obtain_SubGroup_Name(idrs(s,i)),nyrs(s,i),i=1,msol) !nomco1="componente a ser recuperado
        if (nopt.ne.2) then
            write (idev,970)
            read (idevr,980,err=3000) icomp
        endif
        if (icomp.eq.1) then
         go to 960
        end if

	  nrcar=0
        if (ipareq.eq.2) then
!c	formación de comp (para el CAR) a partir de idr y nyr 
		    do 5 i=1,10
			    comp(i,1)=idrs(s,i) !grupos del componente
			    comp(i,2)=nyrs(s,i) !repetición de grupos de componente
   5	    continue
		    call checkdatos (comp,nisom,ident,nombre,formula)!está en"AporteDatos"
!															checkea en prop.pdb
		    if (nisom.ne.0) then !número de isómeros
			    call elegir_CXR (nisom,ident,nombre,formula,nrcar)
		    end if
	  end if
	  if (nrcar.ne.0) then
		 !   call pvomegaybp (nrcar,a1,a2,bp1,vliq)
	  else
!c
!c---Lectura de la temperatura de ebullicion del CAR
		    if (mop.eq.1) then
			    nomco4 = 'component to be recovered:           '
		    else
			    nomco3 = 'more  volatile  component:           '
			    nomco4 = 'less  volatile  component:           '
		    end if
3010		    write (idev,530) nomco4 !give the boiling point...
		    read (idevr,*,err=3010) bp1
!c
!c---Lectura de las constantes de Antoine del CAR
		    if (mop.eq.2.and.nrcar.eq.0) then
!c			    call limp (idev)
!			    call antoine (nomco4,'A1','A2','A3',a1,a2,a3)!carga las constantes 
!															de Antoine
		    end if
        end if
	!!!!!!NO ESTÁ DECLARADA ICANT!!!!!!
	  IF (ICANT.EQ.1) GOTO 350 !(cuando el cpr ya fue ingresado, 
!c	  pero se volvió a cambiar el CAR porque no estaban todos los amn)
113   continue
!c---Lectura del Componente Principal del Refinado (CPR)
1960  if (mop.eq.1) then
		nomco2 = 'Principal component in the raffinate:'
      else
		nomco2 = 'More  volatile  component:           '
      end if
      call leer_comp (ipareq,nomco2,comp)
 350  ifin0 = ifin1

!c---Ingreso del CPR al VGP
      call asignar_grupos_principales (idf,mraf,grupos)
      call armar_grupos_distintos (grupos,mraf,kgrup,ifin0,ifin2)

!c---Control de parametros de interaccion para el CPR
!c     call limp (idev)
      PROB=.FALSE.
      if (ifin2.gt.ifin0) then
		write (idev,1910) !*** Checking interaction parameters ***
		if (model /= 3)then
		    prob = check_inter_param (ipareq,1,kgrup,ifin0+1)
		else
		    prob = .False.
		endif
      end if
      if (prob) then
3020		write (idev,1000)
		write (idev,1020)
		read (idevr,980,err=3020) icomp
		if (icomp.eq.1) then
			ICANT=1
			go to 960
		else
			ICANT=0
			go to 1960
		end if
      end if
!c
!c---Confirmacion del CPR
3030  write (idev,1210) nomco2,(fs(idf(i)),nyf(i),i=1,mraf)
      write (idev,970)
      read (idevr,980,err=3030) icomp
      if (icomp.eq.1) then
         go to 1960
      end if
!c
	nrcpr=0
      if (ipareq.eq.2) then
!c	formación de comp (para el CPR) a partir de idf y nyf 
		do 6 i=1,10
			comp(i,1)=idf(i)
			comp(i,2)=nyf(i)
   6		continue
		call checkdatos (comp,nisom,ident,nombre,formula)
		if (nisom.ne.0) then
			call elegir_CXR (nisom,ident,nombre,formula,nrcpr)
		end if
	end if
	do 961 i=1,mraf
          CALL Leer_Pr (idf(I),IPAREQ,fs(idf(i)),M1,J1,K1,I1,H1,PMG(i),&
     	              RPAR,QPAR,DTC,DPC,DelV,dvi,dn,DTCA)   !en MANUNI1
 961  continue
	pmt=0
	DO 342 I=1,mraf
		PMT = PMT + PMG(I)*nyf(I)
 342  continue   
	if (nrcpr.ne.0) then
!		call pvomegaybp (nrcpr,b1,b2,bp2,vliq)
		Dens2=0.001*PMT/VLIQ
	else
 3011		write (idev,535) nomco2
		read (idevr,*,err=3011) dens2
!c---Lectura de las constantes de Antoine para el CPR
  85		if (mop.eq.2) then
!c		   call limp (idev)
		  ! call antoine (nomco2,'B1','B2','B3',b1,b2,b3)
		end if
      end if
!c
!c---Lectura de las composiciones azeotropicas 
!c        call limp (idev)
  	if (mop.eq.2) then
  80		write (idev,70)
		read (idevr,510,err=80) opt
		i1 = index ('12',opt)
		if (i1.eq.0) then
			go to 80
		end if
		if (opt.eq.'1') then
!c			call limp (idev)
 110			write (idev,90) 
			read (idevr,*,err=110)	tazeo
 120			write (idev,100)
			read (idevr,*,err=120) xazeo
!c			call limp (idev)
 160			write (idev,130) tazeo
			write (idev,140) xazeo
			write (idev,970)
			read (idevr,980,err=160) icomp
			if (icomp.eq.1) then
 170				write (idev,150) 
				read (idevr,510,err=160) iopt
				i1 = index ('12',iopt)
				if (i1.eq.0) then
!					call limp (idev)
					go to 170
				end if
!c				call limp (idev)
				if (iopt.eq.'1') then
 180					write (idev,90) 
					read (idevr,*,err=180) tazeo
					go to 160
				else if (iopt.eq.'2') then
 190					write (idev,100) 
					read (idevr,*,err=190) xazeo
					go to 160
				end if
			end if
		else
			tazeo = 1.0
			xazeo = 0.0
		end if
      end if
!	formatos
2900  format (' ',//,70x,'<ret>',$)
 510  format (a)
 240  format (/,'*** Loading data ***')
 20   format (' Choose the separation operation:',//,&
             29x,'Liquid-Liquid Extraction: 1',/,&
             29x,'Extractive  Distillation: 2')
 30   format (1x,/,60x,'> ',$)
 540  format (' Give the group composition of the two components',&
             /,' to be separated: id1,ny1,id2,ny2,etc.',/,&
             10x,'id1,id2: subgroup identification number',/,&
             10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
  70  format (1x,/,' Does the system form an azeotrope?',//,&
             51x,'Yes: 1',/,51x,'No : 2',//,60x,'> ',$)
  90  format (1x,/,' Give the azeotropic temperature (K)              ',' > ',$)
5555  FORMAT (1X,/,' Give operation temperature (K)                   ',' > ',$)
 100  format (1x,/,' Give the azeotropic composition',/,&
             ' (molar fraction of the less volatile component)   > ',$)
 130  format (' Azeotropic temperature:    ',f8.3,' K')
 140  format (' Azeotropic composition:    ',f8.3)
 150  format (' ',/,' Which variable do you want to change?',//,21x,&
            'tazeo: 1',/,21x,'xazeo: 2',//,42x,'> ',$)
1910  format (' ','*** Checking interaction parameters ***',/)
1000  format (' ',/,' ','* No interaction parameters for some of the ',&
             'subgroups',/,'   already entered.')
1010  format ('   Change the info. about the component:  <ret>. ',$)
1020  format ('   Do you want to change the components?',/,&
             '               Component to be recovered: 1',/,&
             '               Principal component in the raffinate: 2',&
             /,60x,'> ',$)                                                  
 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 980  format (i2)
1210  format (' ',a37,10(a8,i2))
 530  format (1x,/,' Give the boiling point (K) of the ',a37,/,64x,'> ',$)
 535  format (1x,/,' Give the density (gr/ml) of the ',a37,/,64x,'> ',$)
      return 
 501  end
!C ============================================================================
