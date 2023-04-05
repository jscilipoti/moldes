!    subroutine SolvRea ()
!!-------------------------------------------------
!!   Descripción
!!   - Variables de entrada
!!       vars:
!!   - Variables de salida
!!       vars:
!!-------------------------------------------------
!    use GRPS
!    use PROPI
!    parameter(NSCM=10,NSEL=20,NG=70,NMG=150,NMAX=60000)
!    implicit real*8(A-H,O-Z)
!!...Integers
!    integer::maingr(NMG),grupo1(NMG),ids(nsel),nys(nsel),MS(3,10,2),nr,np
!    integer,dimension(:,:,:),allocatable::reagents,products
!!...real*8
!    real*8:: MST(NMAX,NSCM,2),temp
!!...Characters
!    character*8 FS(NMG)
!    character*3 CAR
!    character*37 nomco
!!...Logical
!
!!...Commons
!    common/PAREQ/IPAREQ
!    common/NOM/FS
!    common/PUNSUB/NPUNT(NMG),NGRUP(NMG),NUM
!    common/PUNGRU/NPINT(NG),NINTT(NG),NUMINT
!!-----PRUEBAS
!    DIMENSION kgrup(nsel),ngr1(10),ngr2(10),ngr3(10),ngr4(10)
!    CHARACTER*21 family(8)
!
!!...Sentencias
!!   inicialización de variables GRPS
!    mop=3
!    call ab_ban1                          !Selección de modelo termodinámico y apertura de bases de datos 
!    call tabla_parametros(ipareq)         !Selección de tabla de parámetros	
!    ngrmax = max_sub (ipareq)             !Cantidad de grupos en la tabla de parámetros elegída
!    do J=1,ngrmax                         !Carga del vector fs (nombre de subgrupos)
!	  call carac (j,ipareq,fs(j),CAR)
!	  grupo1(j) = j
!	enddo
!
!!-----Ingresar REACTIVOS--------------------------------------
!	write(6,"('  Give number of reagents> ')")
!    write(6,"(' ',//,70X,'>',$)")
!    read (5,*)nr
!    allocate(reagents(nr,10,2))
!    reagents(:,:,:)=0
!	nomco = 'COMPONENT:                   '
!    do i=1,nr
!	    write(6,*)'Component structure number ',i
!		write(6,601)
!1960    call leer_comp (ipareq,fs,grupo1,ngrmax,nomco,ids,nys,msol)         !se presentan y seleccionan grupos
!3030    write (6,"(' ',a37,10(a8,i2))") nomco,(fs(ids(l)),nys(l),l=1,msol)  !Confirmacion del Componente
!	    write (6,970)
!	    read (5,"(i2)",err=3030) icomp
!	    if(icomp.eq.1)GOTO 1960
!	    do j=1,msol                         !carga del componente ingresado en REAGENTS
!		    reagents(i,j,1) = ids(j)
!			reagents(i,j,2) = nys(j)
!			call CR_PUNTF(IDS(J),NPUNT,NGRUP,NPINT,NINTT) !Carga grupo en NGRUP y NPINT
!  		enddo
!  		!CALL CR_PUNTF(REAGENTS(I,:,:),NPUNT,NGRUP,NPINT,NINTT)
!    enddo
!    write(6,"(/,'The following reagents are entered. Select a reagent for calculation of ,/,mixing properties.',/)")
!    call SeleccionarCompuesto(reagents,nr,ms(1,:,:))
!!-----Ingresar PRODUCTOS--------------------------------------
!	write(6,"('  Give number of products> ')")                          !Pide cantidad de productos
!    write(6,"(' ',//,70X,'>',$)")
!    read (5,*)np
!    allocate(products(np,10,2))
!    products(:,:,:)=0
!	nomco = 'COMPONENT:                   '
!    do i=1,NP
!	    write(6,*)'Component structure number ',i
!		write(6,601)
! 1961	call leer_comp (ipareq,fs,grupo1,ngrmax,nomco,ids,nys,msol) !se presentan y seleccionan grupos
! 3031	write(6,"(' ',a37,10(a8,i2))") nomco,(fs(ids(l)),nys(l),l=1,msol)           !Confirmacion del Componente
!	    write(6,970)
!	    read(5,"(i2)",err=3031) icomp
!	    if (icomp.eq.1)GOTO 1961
!	    do j=1,msol                         !carga del componente ingresado en REAGENTS
!		    products(i,j,1) = ids(j)
!			products(i,j,2) = nys(j)
!			call CR_PUNTF(IDS(J),NPUNT,NGRUP,NPINT,NINTT) !Carga grupo en NGRUP y NPINT
!  		enddo
!  		!CALL CR_PUNTF(PRODUCTS(I,:,:),NPUNT,NGRUP,NPINT,NINTT)
!    enddo
!    call SeleccionarCompuesto(products,np,ms(2,:,:))
!    write(6,"(/,' Give reaction temperature (K)',19X,' > ',$)")
!    read (5,*)temp
!    call Constraints_Selection(temp,3)
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!1200  write (6,590) !Give the family of solvents to be generated: ....
!      read (5,*,err=1200) ifam
!      if ((ifam.lt.1).or.(ifam.gt.4)) then
!		go to 1200
!      end if    
!3050  write (6,1600) !now you have to choose......screen(1);screen and output
!      read (5,"(i2)",err=3050) imprim   
!    CALL SELECCION_GRUPOS (fs,MAINGR,IMPRIM,KGRUP,IFIN2,NGR1,NGR2,NGR3,NGR4,MGR,family) !family
!!    CALL CR_PUNT (MS,MSOL,MRAF,NGDV,MDV,NGSV,MSV,NGSV1,MSV1,&
!!                      IPAREQ,IFAM,IALAR,NALAR)
!    call Store_Pr (ipareq)
!    call Store_In (ipareq)
!    CALL STRUCTURE_GENERATOR (JIST,MST,MS,3,NGSSV,TEMP1,NGSDV,SALIR,.FALSE.)    
!
!!------------------------------------------------    
!    call Store_Pr (ipareq)
!    call Store_In (ipareq)
!    call propi2(reagents(1,:,:),1,0,PMOL=PMOL)
!!-----------------------------------------------
!!...Formats
! 601format(' Give the group composition of the compound: ',&
!     	   ' id1,ny1,id2,ny2,etc.',/,                      &
!           10x,'id1,id2: subgroup identification number',/,&
!           10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
! 970format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
! 590  format (/,' Give the family of solvents to be generated:',///,&
!             10x,'-Aromatic solvents                  : 1',//&
!             10x,'-Single substance groups            : 2',//&
!             10x,'-No-aromatics with up to 12 groups in',/&
!             10x,' the final structure                : 3',//&
!             10x,'-Cyclic solvents                    : 4',///&
!     	      51x,' > ',$)
! 1600  format (' ',/,' * Now you have to choose the groups that take ',&
!             'part in the Molecular',/,&
!                   '   Design of Solvents.  The program performs a ',&
!             'preliminary selection',/,&
!                   '   in  order to  eliminate  those groups which ',&
!             'have  no  interaction',/,&
!                   '   parameters  with  respect  to  the two main ',&
!             'components.',/,&
!                   '   If you want to  see  the  pairings  without ',&
!             'interaction parameters',/,&
!                   '   enter the following code:',/,&
!              22x,'screen:                              1',/,&
!              22x,'screen and output file (nointer):    2',/,&
!              22x,'no printing:                      <ret>',/,&
!              63x,'> ',$)
!    endsubroutine



!    subroutine initialization()
!    use PropertiesData
!    use carteles
!    character*63::cartelt(2,18)
!!===Inicializazión
!    MW(:)      = 0.0
!    BP(:)      = 0.0
!    PC(:)      = 0.0
!    VC(:)      = 0.0
!    TC(:)      = 0.0
!    Visc(:)    = 0.0
!    Dens(:)    = 0.0
!    Hvap(:)    = 0.0
!    Kow(:)     = 0.0
!    SelInt(:)  = 0.0
!    Sel(:)     = 0.0
!    SolPow(:)  = 0.0
!    SLostInt(:)= 0.0
!    SLost(:)   = 0.0
!    DistCoef(:)= 0.0
!    VolatInt(:)= 0.0
!    Volat(:)   = 0.0
!    DTC2S(:)   = 0.0
!!
!    MWB      = (/0.0,0.0/)
!    BPB      = (/0.0,0.0/)
!    PCB      = (/0.0,0.0/)
!    VCB      = (/0.0,0.0/)
!    TCB      = (/0.0,0.0/)
!    ViscB    = (/0.0,0.0/)
!    DensB    = (/0.0,0.0/)
!    HvapB    = (/0.0,0.0/)
!    KowB     = (/0.0,0.0/)
!    SelIntB  = (/0.0,0.0/)
!    SelB     = (/0.0,0.0/)
!    SolPowB  = (/0.0,0.0/)
!    SLostIntB= (/0.0,0.0/)
!    SLostB   = (/0.0,0.0/)
!    DistCoefB= (/0.0,0.0/)
!    VolatIntB= (/0.0,0.0/)
!    VolatB   = (/0.0,0.0/)
!    DTC2SB   = (/0.0,0.0/)
!!
!    LMW      = (/.False.,.False./)
!    LBP      = (/.False.,.False./)
!    LPC      = (/.False.,.False./)
!    LVC      = (/.False.,.False./)
!    LTC      = (/.False.,.False./)
!    LVisc    = (/.False.,.False./)
!    LDens    = (/.False.,.False./)
!    LHvap    = (/.False.,.False./)
!    LKow     = (/.False.,.False./)
!    LSelInt  = (/.False.,.False./)
!    LSel     = (/.False.,.False./)
!    LSolPow  = (/.False.,.False./)
!    LSLostInt= (/.False.,.False./)
!    LSLost   = (/.False.,.False./)
!    LDistCoef= (/.False.,.False./)
!    LVolatInt= (/.False.,.False./)
!    LVolat   = (/.False.,.False./)
!    LDTC2S   = (/.False.,.False./)    
!    LDATO(:,:)=.False.
!    
!    cartelt=reshape((/'minimun molecular weight                                          ',& !1
!        'maximun molecular weight                                          ',&
!        'minimun boiling point                                         (K) ',& !2
!        'maximun boiling point                                         (K) ',& 
!        'minimun critical pressure                                      () ',&  !3
!        'maximun critical pressure                                      () ',&
!        'minimun critical volume                                        () ',&  !4
!        'maximun critical volume                                        () ',&
!        'minimun critical temperture                                    () ',&  !5
!        'maximun critical temperture                                    () ',&
!        'minimun viscosity                                              () ',&  !6
!        'maximun viscosity                                              () ',&
!        'minimun enthalpy of vaporization                               () ',&  !7
!        'maximun enthalpy of vaporization                               () ',&
!        'minimun octanol-water partition coeficient                     () ',&  !8
!        'maximun octanol-water partition coeficient                     () ',&
!        'minimun selectivity of an intermediate structure                  ',&  !9
!        'maximun selectivity of an intermediate structure                  ',&
!        'minimun selectivity of a final structure                     (wt) ',&  !10
!        'maximun selectivity of a final structure                     (wt) ',&
!        'minimun solvent power required                               (wt) ',&  !11
!        'maximun solvent power required                               (wt) ',&
!        'minimun solvent loss of an intermediate structure          (wt %) ',&  !12
!        'maximun solvent loss of an intermediate structure          (wt %) ',&
!        'minimun solvent loss of a final structure                  (wt %) ',&  !13
!        'maximun solvent loss of a final structure                  (wt %) ',&
!        'minimun distribution coefficient required                  (wt %) ',&  !14
!        'maximun distribution coefficient required                  (wt %) ',&
!        'minimun volatility with an intermediate structure                 ',&  !15
!        'maximun volatility with an intermediate structure                 ',&
!        'minimun volatility with a final structucture                      ',&  !16
!        'maximun volatility with a final structucture                      ',&
!        'minimun Minimun boiling temperature difference between solvent and less volatile component (K)',&  !17
!        'maximun Minimun boiling temperature difference between solvent and less volatile component (K)'/)&
!        ,shape(cartelt))
!        do i=1,size(cartelt(:,1))
!            do j=1,size(cartelt(1,:))
!                cartel(j,i)=cartelt(i,j)
!            enddo
!        enddo
!        do i=1,size(cartel(:,1))
!            do j=1,size(cartel(1,:))
!                write(*,*)i,j,cartel(i,j)
!            enddo
!        enddo
!    endsubroutine
!    subroutine ConstraintsSelectionIntStruct(mop)
!!-------------------------------------------------
!!   Descripción
!!   - Variables de entrada
!!       vars:
!!   - Variables de salida
!!       vars:
!!-------------------------------------------------
!    use PropertiesData
!    integer,intent(in)::mop
!    integer::nresp
!    
!    nresp=0
!    write(6,*)' What property do you want to evaluate pre-final structures?'    
!    if(mop==1)then
!        do while (nresp<1.or.nresp>3)
!            write(6,*)' Only selectivity                      1'
!		    write(6,*)' Selectivity and solvent loss          2'
!		    write(6,*)' Only solvent loss                     3'
!	        write(6,"(' ',60x,'<ret>',$)")
!	        read(5,*) nresp !Property to Evaluate Pre-Intermediate Structures
!        enddo
!        if(nresp==1)then
!            LSelInt(1)=.True.
!        elseif(nresp==2)then
!            LSelInt(1)=.True.
!            LSLostInt(2)=.True.
!        else
!            LSLostInt(2)=.True.
!        endif    
!    elseif(mop==2)then
!        do while (nresp<1.or.nresp>3)
!            write(6,*)' Only relative volatility              1'
!		    write(6,*)' Relative volatility and solvent power 2'
!		    write(6,*)' Only solvent power                    3'
!	        write(6,"(' ',60x,'<ret>',$)")
!	        read(5,*) nresp !Property to Evaluate Pre-Intermediate Structures
!        enddo
!        if(nresp==1)then
!            LSCLLI=.True.
!        elseif(nresp==2)then
!            LSCLLI=.True.
!            LSOLV=.True.
!        else
!            LSOLV=.True.
!        endif 
!    endif
!	endsubroutine
!      MODULE PROPI
!      !CONTAINS
!      INTERFACE
!        SUBROUTINE PROPI2 (COMPOUND,MOP,ICOMP,PMOL,TC,PC,VC,BPOINT,VISC,DENSID,HVAP,RELV)
!            IMPLICIT real*8 (A-H,O-Z)
!!           VARIABLES DE ENTRADA
!            INTEGER, INTENT(IN)::COMPOUND(10,2),MOP,ICOMP
!!           VARIABLES DE SALIDA (OPCIONALES)
!            real*8,INTENT(OUT),OPTIONAL::PMOL,TC,PC,VC,BPOINT,VISC,DENSID,HVAP,RELV
!        ENDSUBROUTINE
!      END INTERFACE
!      ENDMODULE
!      SUBROUTINE PROPI2 (COMPOUND,MOP,ICOMP,PMOL,TC,PC,VC,BPOINT,&
!                        VISC,DENSID,HVAP,RELV)
!-------------------------------------------------
!   Descripción: Calcula propiedades físicas de compuesto puro
!   - Variables de entrada
!       compound(10,2):componente
!       mop:    operación de separación (extracción L-L(1); destilación extractiva(2))
!       icomp:  structure type of the component
!		            Aliphatic		        :0
!		            Aromatic		        :1
!		            Single group solvent	:2
!   - Variables de salida
!     - OPCIONALES -
!       PMOL:   peso molecular
!       TC:     temperatura crítica
!       PC:     presión crítica
!       VC:     volumen crítico
!       BPOINT: punto de ebullición
!       VISC:   viscosidad
!       DENSID: densidad
!       HVAP:   calor de vaporización
!       RALV:   volatilidad relativa
!-------------------------------------------------

!      PARAMETER (NMG=150,NCOM=3,NA=150,NSCM=10)
!      IMPLICIT real*8 (A-H,O-Z)
!!     VARIABLES DE ENTRADA
!      INTEGER, INTENT(IN)::COMPOUND(10,2),MOP,ICOMP
!!     VARIABLES DE SALIDA (OPCIONALES)
!      real*8,INTENT(OUT),OPTIONAL::PMOL,TC,PC,VC,BPOINT,VISC,DENSID,HVAP,RELV
!!     NUEVO ARGUMENTO MOP
!
!      DIMENSION PMG(NMG),DELTC(NMG),DELV(NMG),DELPC(NMG),delvi(NMG),deln(NMG),DELTCA(NMG)    
!!     real*8 K11,K22
!      COMMON/PUNSUB/NPUNT(NA),NGRUP(NA),NCANT
!      COMMON/CONT/PMG,DELTC,DELV,DELPC,delvi,deln,DELTCA
!      ! COMMON/US/MS(NCOM,NSCM,2)
!
!      PMOLT = 0.0
!      SUMTC = 0.0
!      SUMPC = 0.0
!      SUMV = 0.0
!      SUMTCA = 0.0
!      i=1
!      do while (COMPOUND(I,1).ne.0)
!         PMOLT = PMOLT + PMG(NPUNT(COMPOUND(I,1)))*COMPOUND(I,2)
!         SUMTC = DELTC(NPUNT(COMPOUND(I,1)))*COMPOUND(I,2) + SUMTC
!         SUMTCA= DELTCA(NPUNT(COMPOUND(I,1)))*COMPOUND(I,2)+SUMTCA
!         SUMV = DELV(NPUNT(COMPOUND(I,1)))*COMPOUND(I,2) + SUMV
!         SUMPC = DELPC(NPUNT(COMPOUND(I,1)))*COMPOUND(I,2) + SUMPC
!         i=i+1
!      enddo
!      numgru=i-1
!      IF (ICOMP.EQ.2) THEN !Single group solvent
!         BPOINTT = SUMTC
!         DENSIDT = SUMV
!         HVAPT = SUMPC
!      ELSE
!         IF (ICOMP.EQ.3) THEN !Cyclic
!            TETA=1-(PMOLT**0.0268/(3.0+SUMTC))+((SUMTCA)/((PMOLT**1.3769)*(0.1219)))
!         ELSE !Alifatic
!            TETA=1-(PMOLT**0.0268/(1+SUMTC))+((SUMTCA)/((PMOLT**1.3769)*(0.1219)))
!         ENDIF           
!         PCT = PMOLT/(SUMPC + .34)**(2.0) !(Ec. 4.10)
!	   IF (ICOMP.EQ.1) THEN !Aromatic
!               SUMV = SUMV + 12
!         ELSE IF (ICOMP.EQ.3) THEN !Cyclic
!               SUMV = SUMV - 7
!         END IF
!         VCT =(SUMV/.285)**(1/1.048) !(Ec. 4.11)
!         ALFA = .9076*(1.0 + TETA*DLOG(PCT)/(1.0 - TETA)) !(Ec. 4.9)
!         TCT = VCT*PCT*(3.72 + .26*(ALFA - 7.0))/82.05 !(Ec. 4.8)
!         BPOINTT = TETA*TCT
!         HVAPT = 1.987*TCT*(TETA*LOG(PCT)/(1.0 - TETA))
!         IF (MOP.EQ.0) DENSIDT = PMOLT/SUMV 
!         IF (MOP.EQ.1) THEN
!            DENSIDT = PMOLT/SUMV
!            DENSC3 = PMOLT/VCT         
!            ZC = PCT*VCT/(82.05*TCT)
!            IF (ZC.LE..26) THEN
!               K22=-3.28257+13.6377*ZC+107.4844*ZC**2-384.211*ZC**3
!            ELSE
!               K22=60.2091-402.063*ZC+501.0*ZC**2+641.0*ZC**3
!            END IF
!            K11=17.4425-214.578*ZC+989.625*ZC**2-1522.06*ZC**3
!         !   ABC = T
!!           DENSIDT=DENSC3*(1+K11*(1-T/TCT)**(1/3)+K22*(1-T/TCT)**(2/3)
!!     *             +(.93-K22)*(1-T/TCT)**(4/3))
!!
!         END IF
!       END IF
!!	La viscosidad se calcula a 25 c
!!
!      
!	T25=25.+273.15
!	IF (T25.LE.BPOINTT) THEN
!	!  CALL Viscosity(3,NUMGRU,ICOMP,T25,BPOINTT,VISCT)
!	ELSE
!	  VISCT=0.0
!	ENDIF
!!     GENERACION DE VARIABLES DE SALIDA
!      IF(PRESENT(PMOL))PMOL=PMOLT
!      IF(PRESENT(TC))TC=TCT
!      IF(PRESENT(PC))PC=PCT
!      IF(PRESENT(VC))VC=VCT
!      IF(PRESENT(BPOINT))BPOINT=BPOINTT
!      IF(PRESENT(VISC))VISC=VISCT
!      IF(PRESENT(DENSID))DENSID=DENSIDT
!      IF(PRESENT(HVAP))HVAP=HVAPT 
!      !IF(PRESENT(RELV)) CALL RELVOLPRO (RELV,BPOINTT)
!      RETURN
!      END
!subroutine Moldes
!!-------------------------------------------------------------------------------
!!     Subrutinas y Funciones Utilizadas
!!     ---------------------------------
!!     Programa modificado para que calcule viscosidad
!!     Del paquete MANUNI:
!!                         NOM_GRU1
!!                         Size_LSubGroups()
!!                         MAINSG
!!                         AINT
!!-------------------------------------------------------------------------------
!!USE
!    use GRPS
!    use PropertiesData
!    use blockdatas
!    use CONSTANTES
!    use PROP_PDB, only:checkdatos,elegir_CXR,pvomegaybp
!    use input_data , only:seleccion_grupos,enter_mixture,rec_inf
!    use Evaluation, only:Constraints_Selection
!    use SubGrupos
!    use Input
!
!    implicit none
!!Variables INTERNAS    
!    !type(Compound),pointer::recorreSolutes    
!    integer::i,j,i0,i1,i2,qLimits,error
!    integer::idev,idevr,lenfile,imprim,imacy,imaar,imagr,ngrmax,mainsg,ifin0,ifin1,ifin2,imp,imasv
!    integer::imasd,imadv,imaev,ima3v,ima4v
!    integer,dimension(NSEL):: mgr,maux,lgrup,mgrup,igrupi,iaux,mdum,mst,ngr,jgrup,nsubgr,ngru1,ngru2,grupos
!    integer,dimension(NMG)::iar,igrp,idv,icyc,isv,isd,i3v,i4v,grupo1
!    integer,dimension(DiffStructGroups)::ngr1,ngr2,ngr3,ngr4,ngr5
!    real*8::temp,tazeo,dtc2s,solv,slsup1,slsupl,dist,pmma,sclli,select
!    real*8,dimension(20)::pmg
!    character*1::opt,iopt,optinv
!    character*8::nom_gru1
!    character*8,dimension(NMG)::fs
!    logical::prob,inter,esc,entro,entro1,entro2,entro3  
!    character*35::nomb       
!!COMMONS
!    common/EVAL/user,inv
!    logical::USER,INV      
!    common/PAREQ/IPAREQ
!    integer::ipareq
!    common/int/ifin2,kgrup,maingr
!    integer,dimension(NSEL)::kgrup
!    integer,dimension(NMG)::maingr
!    common /molcom/arch !arch=false cuando molde2 o rec_inf no encontraron el archivo .mdi
!    logical::arch
!!    common/NGX/NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1
!!    integer::NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1
!!EXTERNAL
!    external max_sub,mainsg,nom_gru1,Leer_In   
!
!!SENTENCIAS
!    idev = 6
!	idevr = 5
!
!!---Pregunta si se corre en modo MOLDES Invertido
!    i1 = 0
!    do while (i1==0)
!        write (6,"(////' ','Select one option :',/&
!                   //,11x,' Run a solvent molecular design problem     : 1',&
!                   //,11x,' Run a solute molecular design problem      : 2',&
!                   //,60x,'> ',$)")
!        read (5,'(a)') optinv
!        i1 = index ('12',optinv)            
!    enddo
!	InputProblem01%inv = .FALSE.
!    if (i1==2) InputProblem01%inv = .TRUE.
!
!!---Opciones de corrida
!    i1=0
!    do while(i1==0)
!        write (idev,"(' ','Options:',//,11x,'Create a new input file : 1',/,11x,'Run a created input file: 2',//,42x,'> ',$)") 
!        read (idevr,"(a)") opt
!        i1 = index ('12',opt)
!    enddo
! 
!    if (i1==1) then !Create a new input file : 1
!    
!! Pide el título del problema    
!2110    write (idev,"('Give problem title (maximum 70 characters)')") 
!        read (idevr,510,err=2110) titl
!        InputProblem01%ProblemTitle = titl
!        
!! Pide el nombre del archivo .mdi        
!3100    write (idev,"(/,' Give input file name (maximum 35 characters).')")
!        read (idevr,510,err=3100) nomb
!        InputProblem01%FileName = nomb
!        
!        entro2 = .false.
!        entro3 = .false.
!        arch = .true.
!  	    call enter_mixture()! (MOP,temp
!!       call car_car (ipareq,fs,isv,isd,idv,i3v,i4v,icyc,iar,igrp,imasv,imasd,imadv,ima3v,ima4v,imacy,imaar,imagr)
!105     entro1 = .false. 
!        
!        ifam=0 !Seleccion de la familia de componentes a generar
!        do while ((ifam.lt.1).or.(ifam.gt.4))
!            write (idev,590) !Give the family of solvents to be generated: ....
!            read (idevr,*) ifam
!        enddo
!        
!        imprim = 0
!        if(inv)then
!            write (idev,1601) !now you have to choose
!        else      
!3050        write (idev,1600) !now you have to choose......screen(1);screen and output
!        endif
!        read (idevr,980,err=3050) imprim
!        write (idev,1610) family(ifam)
!        if (imprim.eq.2) then
!            lenfile = index(nomb,' ')-1
!            open (unit=2,file=nomb(1:lenfile)//'.par',form='formatted')
!		    write (2,5000)
!		    write (2,1610) family(ifam)
!        end if
!        entro = .false.
!        
!1320    CALL SELECCION_GRUPOS(IMPRIM,KGRUP,IFIN2,NGR1,NGR2,NGR3,NGR4,InputProblem01%MainGroups,family)
!
!300     continue
!
!        call Constraints_Selection()
!
!        if (imprim.eq.2) close (unit=2)
!        
!        lenfile = index(nomb,' ')-1    
!        open(unit=1,file=nomb(1:lenfile)//'.mdi',form='formatted')
!
!!ESCRITURA INPUT FILE    
!        call Write_Input_File()
!						 
!    else !Run a created input file: 2
!
!        write (idev,2060)
!		read (idevr,510,err=375) nomb
!        lenfile = index(nomb,' ')-1
!		open(unit=1,file=nomb(1:lenfile)//'.mdi',status='old',form='formatted',iostat=error)
!		
!		do while (error/=0)
!		    write (6,*) '*ERROR* EN LA APERTURA DE "', nomb(1:lenfile)//'.mdi','"'
!375         write (idev,2060)
!		    read (idevr,510,err=375) nomb
!            lenfile = index(nomb,' ')-1
!		    open(unit=1,file=nomb(1:lenfile)//'.mdi',status='old',form='formatted',iostat=error)		
!		enddo
!	
!400	    call ab_ban1(InputProblem01%model)  !Apertura del banco de datos BADAUN 	        
!    endif
!    
!	USER = .FALSE.
!    arch = .true.
!    
!    call molde2 (nomb,opt)
!
!2080 if (arch) then
!2068    write(idev,2069)
!        read (5,510) opt
!		i0 = index ('12',opt)
!		if (i0.eq.1) then
!			user=.true.
!			call molde2 (nomb,opt)
!		else 
!			if (i0.eq.0) then 
!				goto 2068
!			end if
!		end if
!		if (mop.eq.2) then
!			write (idev,2070) nomb
!		else if (mop.eq.1) then
!			write(idev,2073) nomb
!		end if
!		read (5,510) opt
!		i1 = index ('123456789',opt)
!!c		call limp (idev)
!		if (i1.eq.0) then
!			go to 2080
!		else if (i1.eq.1) then
!			go to 400
!		else if ((i1.eq.2).or.(i1.eq.3).or.(i1.ge.7)) then
!			if (.not.(entro2)) then
!!c				call limp (idev)
!				write (idev,*) '*** Reading data from ',nomb(1:lenfile)//'.mdi'
!				call limp (idev)
!				call rec_inf (nomb,arch)
!				if (mop.eq.2) then
!					dato2(1) = SCLLI 
!					dato2(2) = SELECT
!					dato2(3) = SOLV  
!					dato2(4) = PMMA
!				end if
!			end if
!			if (arch) then
!				write (idev,340)
!  390				read (idevr,510) opt
!				i2 = index ('12',opt)
!				if (i2.eq.0) then
!					go to 390
!				else if (i2.eq.1) then
!  380					write (idev,360)
!					read (idevr,510,err=380) nomb
!				end if
!				if ((i1.eq.2).or.(i1.ge.8)) then
!					entro2 = .true.
!					go to 300
!				else if (i1.eq.7) then
!					if (.not.(entro2)) then
!						call car_car (ipareq,isv,isd,idv,i3v,i4v,icyc,iar,igrp,imasv,imasd,imadv,ima3v,ima4v,imacy,imaar,imagr)
!					end if
!					!goto 1860
!				else if (i1.eq.3) then
!					if (.not.(entro3)) then
!!c
!!c---Lectura de los datos segun la seleccion de parametros
!						write (idev,240) 
!						ngrmax = Size_LSubGroups()
!						do 315 i=1,ngrmax
!							grupo1(i) = i
!							maingr(i) = mainsg (i,ipareq)
! 315				    continue
!						call car_car (ipareq,isv,isd,idv,i3v,i4v,icyc,iar,igrp,imasv,imasd,imadv,ima3v,ima4v,imacy,imaar,imagr)
!						ifin0 = 0
!						call asignar_grupos_principales (idr,msol,grupos)
!						call armar_grupos_distintos (grupos,msol,kgrup,ifin0,ifin1)
!						ifin0 = ifin1
!						call asignar_grupos_principales (idf,mraf,grupos)
!						call armar_grupos_distintos (grupos,mraf,kgrup,ifin0,ifin2)
!					end if
!					entro2 = .true.
!					entro3 = .true.
!               		if (i1.eq.3) then
!						go to 105
!					else
!						goto 1320
!					end if
!				end if
!			else
! 314				write (idev,2060)
!				read (idevr,510,err=314) nomb
!				entro2 = .false.
!				entro3 = .false.
!				arch = .true.
!				go to 400
!			end if
!		else if (i1.eq.4) then
!			call ci_ban1
!			CLOSE(UNIT=IMP)
!			go to 2110
!		else if (i1.eq.5) then
! 310			write (idev,2060)
!			read (idevr,510,err=310) nomb
!			entro2 = .false.
!			entro3 = .false.
!			go to 400
!		end if
!	else
! 316		write (idev,2060)
!		read (idevr,510,err=316) nomb
!		entro2 = .false.
!		entro3 = .false.
!		arch = .true.
!		go to 400
!	end if
!10000 call ci_ban1
!      CLOSE(UNIT=IMP)
!
!!c---- formatos
!2069	format ('  Do you want to evaluate other specific ', 'solvents?',//,50x,'yes: 1',/,51x,'no: 2',/,60x,'>',$)
!   1	format ('  Choose the restriction level for the Molecular Design',' of Solvents:',//,&
!     		'  To use the original Brignole Feasibility Criteria ',&
!     		/,'  (better UNIFAC predictions):',50x,'1',/,&
!     		'  To perform a "less restricted" synthesis by using the',&
!     		' modified BFC (more ',&
!     		/,'  solvent structures but possibly reactive and',&
!     		' less reliable predictions):',6x,'2')
!
! 360  format (' ',//,' Enter the new input file name ','(maximum 35 characters): ',$)
! 340  format (' Do you want to keep the last results file?',//,&
!             '                                          Yes: 1',/,&
!             '                                          No : 2',/,&
!             60x,'> ',$)
!5000  format (' This file has UNIFAC pairs of groups which have no interaction parameters.',//)
!2070  format (' ','Options:',/&
!            /,11x,'Run the ',a6,' separation problem with',&
!            /,11x,'                               the same data: 1',&
!            /,11x,'            different restriction parameters: 2',&
!            /,11x,'    a new solvent family or different groups: 3',&
!            /,11x,'            different combination properties: 7',&
!            /,11x,'different evaluation of pre-final structures: 8',&
!            /,11x,'Create a new input file                     : 4',&
!            /,11x,'Run another file                            : 5',&
!            /,11x,'Exit                                        : 6',&
!    	      //,60x,'> ',$)
!2073  format (' ','Options:',/&
!             /,11x,'Run the ',a6,' separation problem with',&
!             /,11x,'                               the same data: 1',&
!             /,11x,'            different restriction parameters: 2',&
!             /,11x,'    a new solvent family or different groups: 3',&
!             /,11x,'            different combination properties: 7',&
!             /,11x,'different evaluation of pre-final structures: 8',&
!!c    *        /,11x,'            different restriction level: 8',&
!     	    /,11x,'                        other temperature   : 9',&
!             /,11x,'Create a new input file                     : 4',&
!             /,11x,'Run another file                            : 5',&
!             /,11x,'Exit                                        : 6',&
!     	      //,60x,'> ',$)
!
!2060  format (' ',//,' Enter the input file name (maximum 35 characters): ',$)
!1910  format (' ','*** Checking interaction parameters ***',/)
! 410  format (' ',/,' Do you want to change some group combination property ? (Y/N)  > ',$)
!1600  format (' ',/,' * Now you have to choose the groups that take ',&
!             'part in the Molecular',/,&
!                   '   Design of Solvents.  The program performs a ',&
!             'preliminary selection',/,&
!                   '   in  order to  eliminate  those groups which ',&
!             'have  no  interaction',/,&
!                   '   parameters  with  respect  to  the others ',&
!             'main components.',/,&
!                   '   If you want to  see  the  pairings  without ',&
!             'interaction parameters',/,&
!                   '   enter the following code:',/,&
!              22x,'screen:                              1',/,&
!              22x,'screen and output file (nointer):    2',/,&
!              22x,'no printing:                      <ret>',/,&
!              63x,'> ',$)
!1601  format (' ',/,' * Now you have to choose the groups that take ',&
!             'part in the Molecular',/,&
!                   '   Design of Solutes.  The program performs a ',&
!             'preliminary selection',/,&
!                   '   in  order to  eliminate  those groups which ',&
!             'have  no  interaction',/,&
!                   '   parameters  with  respect  to  the others ',&
!             'main components.',/,&
!                   '   If you want to  see  the  pairings  without ',&
!             'interaction parameters',/,&
!                   '   enter the following code:',/,&
!              22x,'screen:                              1',/,&
!              22x,'screen and output file (nointer):    2',/,&
!              22x,'no printing:                      <ret>',/,&
!              63x,'> ',$)
!1610  format (' ',/,' *** ',a21,' preliminary selection ***')
!1500  format (' ',/,' * No interaction parameters between the aromatic',&
!             ' part of the',/,'   structure of the solvent and the ',&
!             'dual valence groups already',/,'   entered.')
!1510  format (' ',/,' * No interaction parameters between the aromatic',&
!            ' part of the',/,'   structure of the solvent and the ',&
!            'single valence groups',/,'   already entered.')
!1490  format (' ',/,' * No interaction parameters available for ',&
!             'this family of solvents',/,'   in the ',a17,&
!             ' UNIFAC table')
!1040  format ('   Change the family of groups:  <ret>. ',$)
!1080  format (' ',/,' ','* No interaction parameters for some of the ',&
!            'subgroups',/,'   already entered and the two main ',&
!            'components.',/,'   Please change the subgroups: <ret>. '&
!            ,$)
!1090  format ('   Do you want to change some subgroups?',/,&
!             '                         First subgroups entered:  1',/,&
!             '                         Single valence subgroups: 2',&
!             /,60x,'> ',$)
! 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
! 800  format (' ',a63,f10.4)
! 801  format (' 12 Tem : Operation temperature                   ','          (K)',f8.2)
! 802  format (' ',/,' Give new operation temperature                 ','                    > ',$)
! 810  format (' ',a64,1x,i8)
! 990  format (' ',a63,2x,i8)
! 820  format ('  ')
! 830  format (' ','Code',3x,'Parameter',43x,'Default Value',/)
! 840  format (' ',/,' If you want to change a parameter value give code',' number.','  If not 0.',/,71x,'> ',$)
! 860  format (' ',a63,7x,'> ',$)
! 870  format (' ',a64,6x,'> ',$)
! 890  format (' ',a63,7x,'> ',$)
! 880  format (' ','Code',3x,'Parameter',47x,'New Value',/)
! !610  format (a70)
!! 620  format (20i3)
!! 630  format (6i3)
!! 625  format (11i3)
! !626  format (3i3)
! 500  format (2X,'****************************************************',&
!             '************'&
!             ,//,27X,'MOLDES 3.0',&
!             //,6X,' Based on the algorithm of Gani et al. (1983),',/,&
!             6X,' Brignole et al. (1986) and Pretel et al. (1994).',//,&
!             6x,' Present version prepared by Patricio Araya, ',&
!             'Eduardo Pretel',/,6X,&
!             ' and completed by Martin Cismondi from the',&
!     		' original version',/,6x,' of E. A. Brignole'&
!             ,//,6X,&
!             ' July, 2000.',//,6X,&
!             ' Planta Piloto de Ingenieria Quimica, 8000 Bahia Blanca,'&
!             ,/,6X,&
!             ' Argentina.',//,2X,&
!             '***************************************************',&
!             '************',//,&
!             '  MOLDES 3.0 is a Group Contribution Methodology to ',&
!             'Computer Aided',/,&
!     	      '  Molecular  Design.    It   selects  solvents  for ',&
!             ' Liquid-Liquid',/,& 
!             '  Extraction  or  Extractive  Distillation  Separati',&
!             'ons of  binary',/,& 
!             '  mixtures.   The selection  is  based on  UNIFAC   ',&
!             'prediction  of',/,&
!             '  solvent properties and on physical and molecular  ',&
!             'constraints.')
! 510  format (a) 
!2900  format (' ',60x,'<ret>',$)
! 240  format (/,'*** Loading data ***')
! 980  format (i2)
! 680  format (20i3)
! 590  format (/,' Give the family of solvents to be generated:',///,&
!             10x,'-Aromatic structures                : 1',//&
!             10x,'-Single substance groups            : 2',//&
!             10x,'-No-aromatics with up to 12 groups in',/&
!             10x,' the final structure                : 3',//&
!             10x,'-Cyclic structures                  : 4',///&
!     	      51x,' > ',$)
!
!    return
!    endsubroutine Moldes
!
!
!
!
!subroutine molde2 (NOMB,opt)
!!C	M O L E C U L A R    D E S I G N    O F   S O L V E N T S
!!C
!!C     SUBROUTINE MODIFICADA PARA EL CALCULO DE VISCOSIDAD         
!!C
!!Módulos
!      use GRPS
!      use PropertiesData
!      use PureProp
!      use StructuresDesign
!      use design
!      use CONSTANTES
!      use Input
!      use input_data, only:Leer_Comp,enter_mixture,ingresar_solventes,rec_inf
!      use Evaluation, only:Constraints_Selection
!      use SubGrupos
!
!!Parámetros
!      parameter (NMOD=3,NMSV=20,NA=150,NSVA=20,NSCM=10,NCOM=3,NESTM=7100,IS1=5,NAA=4)
!
!      implicit real*8 (A-H,O-Z)
!
!!Interfases 
!      interface
!        subroutine score(punt,nsol)
!            use StructuresDesign
!            implicit none
!        !..Variables de ENTRADA/SALIDA
!            type(FinalStructure),pointer,intent(inout)::punt
!            integer,intent(in)::nsol
!        !Variables INTERNAS
!            type(FinalStructure),pointer::recorre
!            integer::i
!        endsubroutine score
!      endinterface
!      
!      interface
!        subroutine write_results (FMSs,mop,idev,ifam,nsol)
!            use SubGrupos
!            use StructuresDesign
!            use constantes
!            implicit none
!        !Varibles de ENTRADA
!            type(FinalStructure),pointer,intent(in)::FMSs
!            integer,intent(in)::mop,idev,ifam,nsol
!        endsubroutine write_results 
!      endinterface
!
!!Tipo de datos
!      type(FinalStructure),pointer::FMSs
!      type(FinalStructure),pointer::recorre
!      type(Compound),pointer::recorreSolutes
!      type(Compound),pointer::ptr_pcp
!      
!      EXTERNAL Leer_In,NCOMBIN
!      LOGICAL SOLACE,ARCH,IFIND,SALIR,mezcla,puro,FinalStruct
!      LOGICAL USER,INV,entro2,entro3
!	  CHARACTER*8 FS(NA),FSF(NA)
!      CHARACTER*70 TITULO
!      CHARACTER*1 SIGUE,OPT
!      CHARACTER*35 NOMB
!      character*37 nomco1
!	  integer UNIF(20),grupo1(na),maingr(na)
! 	  CHARACTER*35 ANAME,FORM2
!	  CHARACTER*20 FORM1
!	  integer comp(20),ident(8),identsolv(NMAX,8),nisomsolv(NMAX)
!	  INTEGER IWANT(10,2),ISUS(10,2),H
!	  character*35 nombre(8),formula(8) 
!	  character*35 nombresolv(NMAX,8),formulasolv(NMAX,8) 
!	  INTEGER NUD(NSVA,NMAX)
!      real*8 TB(NMAX,8),DENST(NMAX,8),VLIQ(NMAX,8)
!      real*8 TMP, bpt
!      real*8 RDER(NMAX)
!      real*8 VREL(NMAX),VRELSOLV(NMAX,8)
!      real*8 MSMW(3)
!	  dimension isv(NMG),isd(NMG),idv(NMG),icyc(NMG),iar(NMG),igrp(NMG),grupos(nsel),kgrup(nsel),i3v(NMG),i4v(NMG)
!	  dimension idrf(nsel),nyrf(nsel)
!      DIMENSION ISOL(NMAX,2)
!      DIMENSION IGR0(NMAX,2),IGR1(NMAX,2),IGR2(NMAX,2)
!      DIMENSION IGR3(NMAX,2),IFSOL(NMAX,2)
!      DIMENSION NGSDV(NSVA,5)
!      DIMENSION NGSSV(NSVA,5)
!      DIMENSION MS(NCOM,NSCM,2)
!      DIMENSION MST(NMAX,NSCM,2),nicomp(nmax)
!      DIMENSION nmst(nmax),IACCOO(NMOD)
!      DIMENSION NPUNT(NA),NGRUP(NA)
!      DIMENSION MM(NMG),MJ(NMG),MK(NMG),MI(NMG),MH(NMG),PMG(NMG)
!      DIMENSION DELTC(NMG),DELV(NMG),DELPC(NMG),delvi(NMG),deln(NMG)
!      DIMENSION ICH2(NMOD),ICH(NMOD),IACOH(NMOD),IACCH2(NMOD)
!      DIMENSION IACCH(NMOD),IAC(NMOD),IACH(NMOD),IACCL(NMOD),DELTCA(NMG)
!      DIMENSION IACNH2(NMOD),IACNO2(NMOD),ICH3(NMOD),IOH(NMOD)
!      DIMENSION ICOOH(NMOD)
!      common/int/ifin2,kgrup,maingr
!	  COMMON/gru/GRUPO1,NGRMAX
!	  COMMON/NOM/FS
!      COMMON/PUNSUB/NPUNT,NGRUP,NCANT
!      COMMON/PUNGRU/NPINT(NINT),NINTT(NINT),NUMINT
!      COMMON/US/MS
!	  COMMON/CONT/PMG,DELTC,DELV,DELPC,DELVI,DELN,DELTCA
!      COMMON/GRUPAR/IACH,IACCH2,IACCL,IACNH2,IACNO2
!      COMMON/GRUPAL/ICH3,ICH2
!      COMMON/GRUPES/ICH,IACOH,IACCH,IAC,IOH,ICOOH,IACCOO
!      COMMON/STAT/MM,MJ,MK,MI,MH
!      COMMON/EXTDIS/PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
!      COMMON/MOLCOM/ARCH
!	  COMMON/SOLVISOM/NISOMSOLV,IDENTSOLV,NOMBRESOLV,FORMULASOLV,TB,DENST,VRELSOLV
!      COMMON/EVAL/USER,INV
!      COMMON/PAREQ/IPAREQ
!      common/NGX/NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1
!
!      IDEV = 6
!
!!C     INITIALIZATION OF MS MATRIX
!50  MS(:,:,:)=0
!!C
!    !lee inputfile
!    if(opt=='2')then
!        call rec_inf(nomb,arch)	    
!        if(.not.arch)then
!            goto 2020
!            write(6,*) '*****',nomb(1:lenfile)//'.mdi',': file not found *****'
!            write(6,1010)
!            read(5,890) SIGUE
!        endif
!    endif
!    
!
!
!	if(user) then !evaluar solventes 
!		InputProblem01%kilout=0
!		
!		call car_car (ipareq,isv,isd,idv,i3v,i4v,icyc,iar,igrp,imasv,imasd,imadv,n3v,n4v,imacy,imaar,imagr)
!		if(ngrmax.eq.0)then
!		ngrmax = Size_LSubGroups()
!			do 19 i=1,ngrmax
!				grupo1(i) = i
!				maingr(i) = mainsg (i,ipareq) !función q devuelve n° de main group
!  19			continue
!		end if
!		call ingresar_solventes (ipareq,jist,mst,nmst,nicomp)
!!---------Ahora los grupos que forman los solventes ingresados son 
!!		copiados en NGDV para que luego sean tenidos en cuenta por
!!		CR_PUNT y sus parámetros sean leídos por Store_Pr
!		mdv=0
!		do 147 i=1,jist
!			do 148 j=1,nmst(i)
!				mdv=mdv+1
!				ngdv(mdv)=mst(i,j,1)
!  148			continue
!  147		continue
!        CALL CR_PUNT (MS,MSOL,MRAF,NGDV,MDV,NGSV,MSV,NGSV1,MSV1,IPAREQ,IFAM,IALAR,NALAR)
!		CALL Store_Pr (IPAREQ)
!        CALL Store_In (IPAREQ)
!		WRITE (6,1010)
!		READ (5,890) SIGUE
!		OPEN (UNIT=2,FILE='SOLVENTS.MDO',FORM='FORMATTED')
!	else    
!        lenfile = index(nomb,' ')-1
!		OPEN (UNIT=2,file=nomb(1:lenfile)//'.MDO',FORM='FORMATTED')
!
!
!        i=1
!        do while(ngdv(i)/=0)
!            call CR_PUNTF(ngdv(i),NPUNT,NGRUP,NPINT,NINTT)
!            i=i+1
!		enddo
!		
!
!        i=1
!        do while(ngsv(i)/=0)
!            call CR_PUNTF(ngsv(i),NPUNT,NGRUP,NPINT,NINTT)
!            i=i+1
!		enddo
!		
!
!        i=1
!        do while(ngsv1(i)/=0)
!            call CR_PUNTF(ngsv1(i),NPUNT,NGRUP,NPINT,NINTT)
!            i=i+1
!		enddo
!
!		if(opt=='2')call Constraints_Selection()
!		ERROR = 1.E-3
!        WRITE (6,*) '*** MOLDES RUN ***'
!        WRITE (6,*) ' '
!        WRITE (6,*) '*** LOADING DATA ***'
!        !do i=1,2
!        recorreSolutes => InputProblem01%MixtureInput%Solutes
!        do while(associated(recorreSolutes))
!            k=1
!            do while(recorreSolutes%Formula(k,1)/=0)
!                call CR_PUNTF(recorreSolutes%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
!                k=k+1
!            enddo
!            recorreSolutes => recorreSolutes%next
!        enddo
!        k=1
!        do while (InputProblem01%MixtureInput%PCR%Formula(k,1)/=0)
!            call CR_PUNTF(InputProblem01%MixtureInput%PCR%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
!            k=k+1
!        enddo
!        !enddo
!		CALL Store_Pr (InputProblem01%ipareq)
!        CALL Store_In (InputProblem01%ipareq)
!        !Tipos de enlaces de cada subgrupo elegido para el diseño molecular
!		DO 6 I=1,MDV	!MDV: número de subgrupos intermedios
!			NGSDV(I,1)=Obtain_attM(NGDV(I))
!			NGSDV(I,2)=Obtain_attJ(NGDV(I))
!			NGSDV(I,3)=Obtain_attK(NGDV(I))
!			NGSDV(I,4)=Obtain_attI(NGDV(I))
!			NGSDV(I,5)=Obtain_attH(NGDV(I))
!   6    CONTINUE
!	    DO 7 I=1,MSV	!MSV número de grupos terminales
!			NGSSV(I,1)=Obtain_attM(NGSV(I))
!			NGSSV(I,2)=Obtain_attJ(NGSV(I))
!			NGSSV(I,3)=Obtain_attK(NGSV(I))
!			NGSSV(I,4)=Obtain_attI(NGSV(I))
!			NGSSV(I,5)=Obtain_attH(NGSV(I))
!   7    CONTINUE
!
!81      CONTINUE 
!        SOLV = SOLV/100.	 !SOLV  : Minimum solvent power required (wt)
!        SLSUP1 = SLSUP1/100. !SLSUP1: SolventLoss of an intermediate structure      
!        SLSUPL = SLSUPL/100. !SLSUPL: SolventLoss of a final structure 
!
!        IF(IFAM.EQ.5)THEN !IFAM = 5 --> Cyclic Solvents
!           DELTC(NPUNT(ICH2(IPAREQ)))=.013
!           DELPC(NPUNT(ICH2(IPAREQ)))=.184
!           DELPC(NPUNT(ICH(IPAREQ)))=.192
!        END IF
!        IF (IFAM.EQ.2) THEN
!           ICOMP = 2
!        ELSE IF ((IFAM.EQ.1).OR.(IFAM.EQ.4)) THEN
!           ICOMP = 1
!        ELSE IF (IFAM.EQ.5) THEN
!           ICOMP = 3
!        ELSE !IF (IFAM .EQ. 3) THEN
!           ICOMP = 0
!        END IF
!      end if
!
!!C         RAFFINATE AND SOLUTE PROPERTIES
!!      PM2=0
!!      i=1
!!      do while(MS(2,I,1)/=0) !CPR
!!        PM2 = PMG(NPUNT(MS(2,I,1)))*MS(2,I,2)+PM2
!!        i=i+1
!!      enddo
!!
!!      PM1=0
!!      i=1
!!      do while (MS(1,I,1)/=0)
!!		PM1  = PMG(NPUNT(MS(1,I,1)))*MS(1,I,2)+PM1
!!        i=i+1
!!      enddo
!
!!c	PM1 and PM2 are the molecular weights of the CAR and CPR
!      IF (MOP.EQ.2) THEN !Destilación extractiva
!		IF (TAZEO.EQ.1.0) THEN
!			TEMP1 = BP1		!CAR normal boiling point
!		ELSE
!			TEMP1 = TAZEO	!normal azeotropic temperature
!		END IF
!		PSAT1 = EXP(A1 -A2/(TEMP1 + A3))
!		PSAT2 = EXP(B1 -B2/(TEMP1 + B3))
!      ELSE ! Extracción líquido-líquido
!		TEMP1 = InputProblem01%T   !operation temperature
!      END IF
!
!!C     PRINTING OF GROUP PARAMETERS FOR THE KEYS
!      MODEL=0			   ! ?
!      CALL ESCLIS (2,InputProblem01%ipareq,MSOL,MRAF,NGDV,NGSDV,MDV,NGSV,NGSSV,MSV,Top,BP1,TITL,A1,A2,A3,B1,B2,B3,TAZEO,X1AZEO,IFAM)
!      IF (InputProblem01%kilout.EQ.0) THEN
!		CALL ESCLIS (6,InputProblem01%ipareq,MSOL,MRAF,NGDV,NGSDV,MDV,NGSV,NGSSV,MSV,Top,BP1,TITL,A1,A2,A3,B1,B2,B3,TAZEO,X1AZEO,IFAM)
!      END IF
!
!	if (.not.user) then	!generar solventes a partir de grupos preseleccionados
!	    ms(1,:,:)=InputProblem01%MixtureInput%Solutes%Formula(:,:) !Esto es para un solo soluto MODIFICAR!
!	    ms(2,:,:)=InputProblem01%MixtureInput%PCR%Formula(:,:)
!	    CALL STRUCTURE_GENERATOR (JIST,FMSs,MS,3,NGSSV,NGSDV,SALIR,.FALSE.)
!	    IF(SALIR)GOTO 1111
!	else ! para solventes ingresados por el usuario
!		TCCH2=DELTC(NPUNT(ICH2(IPAREQ)))
!		PCCH2=DELPC(NPUNT(ICH2(IPAREQ)))
!		PCCH= DELPC(NPUNT(ICH(IPAREQ)))
!		DO 701 J=1,JIST	
!			if (nicomp(j).eq.3) then
!				DELTC(NPUNT(ICH2(IPAREQ)))=.013
!				DELPC(NPUNT(ICH2(IPAREQ)))=.184
!				DELPC(NPUNT(ICH(IPAREQ)))=.192
!			else
!				DELTC(NPUNT(ICH2(IPAREQ)))=TCCH2
!				DELPC(NPUNT(ICH2(IPAREQ)))=PCCH2
!				DELPC(NPUNT(ICH(IPAREQ)))= PCCH
!			end if
!			MS(3,:,:)=0
!			DO 702 L=1,NMST(J)
!				MS(3,L,1) = MST(J,L,1)
!				MS(3,L,2) = MST(J,L,2)
! 702			CONTINUE
!			DO L=NMST(J)+1,10
!				MS(3,L,1) = 0
!				MS(3,L,2) = 0
!			END DO
!
!			mezcla=.True.
!			Puro=.True.
!            recorreSolutes => InputProblem01%MixtureInput%Solutes
!			do while (associated(recorreSolutes))
!			    ms(1,:,:) = recorreSolutes%Formula				
!!                call Evaluate(FinalStruct,Mezcla,Puro,MS,NC,TEMP1,ICOMP,
!!     *               NSEW,NSWIP,NSBI,SOLACE,FMSs)
!                if(.not.solace)exit
!                recorreSolutes => recorreSolutes%next
!            enddo              
! 
!            if(solace)jist = jist + 1
!701		CONTINUE
!	end if
!    N=0
!    DO 710 J=1,JIST
!		ISOL(J,1) = J
!		ISOL(J,2) = J
! 710	CONTINUE
!    PRINT*,CHAR(7)
!    WRITE (6,1010)
!    READ (5,890) SIGUE
!    DO 700 J=1,JIST-1
!	    if (nisomsolv(j).gt.0) then
!			do 802 n=1,nisomsolv(j)
!				DENST(J,N)=0.001*MW(J)/(VLIQ(j,n)*DENS2)
! 802			continue
!		end if					
!	  IF ((.not.user).and.MOP.EQ.1) THEN
!			ISOL(J,1)=J
!			ISOL(J,2)=J
!			MS(3,:,:)=0
!            DO 800 I=1,NGSL(J)
!		        MS(3,I,1)=MST(J,I,1)
!                MS(3,I,2)=MST(J,I,2)
! 800			CONTINUE
!		END IF
! 700	CONTINUE
!      ISIZE=J-1
!      call score(FMSs,nsol)
!803   CONTINUE
!      call write_results (FMSs,mop,6,ifam,nsol)
!      call write_results (FMSs,mop,2,ifam,nsol)
!
!1111    CLOSE (UNIT=2)
!!C	CLOSE (UNIT=4)
! 2020	continue
!! IF (IERROR.NE.0) THEN
!!        WRITE (6,*) '*****',nomb(1:lenfile)//'.mdi',': file not found *****'
!!        WRITE (6,1010)
!!        READ (5,890) SIGUE
!!        ARCH = .FALSE.
!!      END IF
!
!!C     FORMAToS
! 911  FORMAT (//,"*** THE COMPOUND IS WITHIN THE LIST ACCEPTED"," STRUCTURES ***",//,"    to view compound press <ret>",//)
! 912  FORMAT (//,"*** THE COMPOUND IS NOT WITHIN THE LIST ACCEPTED"," STRUCTURES ***",//)   
! 913  FORMAT (' Give the group composition of the component',&
!             /,' search: id1,ny1,id2,ny2,etc.',/,&
!             10x,'id1,id2: subgroup identification number',/,&
!             10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
!!2030  FORMAT (' ','*****',A6,'.MDI: FILE NOT FOUND *****')
!1010  FORMAT (' ',//,70X,'<RET>',$)
!11    FORMAT (A70)
!10    FORMAT (5I3)
!15    FORMAT (20I3)
!16    FORMAT (11I3)
!18    FORMAT (3I3)
!40    FORMAT (' ',/,10I3)
!46    FORMAT (' ','SELECTIVITY',1F10.4,/,' SOLV.POWER',1F10.4,/)
!890   FORMAT (A)
!	RETURN
!      END
