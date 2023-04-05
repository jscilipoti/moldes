MODULE PureProp
  use constantes
  CONTAINS
!===================================================================
  subroutine Calc_Prop_Pure(compuesto,icomp,Top,MWT,TCT,PCT,VCT,BPT,VISCT,DENST,HVAPT,PvapT,FuncGroup)
    use Input
    use Subgrupos
    implicit none
!   VARIABLES DE ENTRADA
    integer,dimension(DiffStructGroups,2),intent(in)::compuesto
    integer,intent(in)::ICOMP !tipo de estructura molecular
    real*8,intent(in)::Top
!   Variables INTERNAS
    real*8::BPTC,MW,TC,PC,VC,BP,VISC,DENS,HVAP,Volat,Pvap,Dc,R
    integer::NG,i
    integer,dimension(15,15)::nu
    real*8,dimension(15)::q, Tstr, gstr, gdot, gddot  ,Tb ,  TcGC, PcGC, BpGC, dcGC
    real*8::kstrT(15,15),kdotT(15,15),alphaT(15,15)
    integer::gruposgc(15)
    logical::fallo
    
!   VARIABLES DE SALIDA (OPCIONALES)
    real*8,intent(out),optional::MWT,TCT,PCT,VCT,BPT,VISCT,DENST,HVAPT,PvapT
    integer,intent(out),optional::FuncGroup
    COMMON/INTER/A,kstr,kdot,alpha
    real*8::A(NMG,NMG),kstr(NMG,NMG),kdot(NMG,NMG),alpha(NMG,NMG)
!Sentencias---------------------------------------------------------
   ! compuesto(:,1) = (/1,2,0,0,0,0,0,0,0,0/) !temporal
   ! compuesto(:,2) = (/2,1,0,0,0,0,0,0,0,0/) !temporal

    MW = MolecularWeight (compuesto)
    VC = CriticalVolume (compuesto,ICOMP)
    PC = CriticalPressure(compuesto,MW)
    BPTC= TbTcRatio(compuesto,icomp,MW)
    BP = NormalBoilingPoint (VC,PC,BPTC)
    i=1
    do while(compuesto(i,1)/=0)
        i=i+1
    enddo
    if(i==2 .and. compuesto(1,2)==1)then 
        tc= Obtain_dTC(compuesto(1,1))
    else
        TC = BP/BPTC
    endif
    if(Model == 3 .and. present(BPT))then
        R = 82.05D0 !atm·cm3/mol K
        Dc = (0.08943 * R * Tc / Pc)**(1./3.)
        DcGC(1) = Dc

        NG = 0
        nu = 0
        Tstr = 0
        gstr = 0
        gdot = 0
        gddot = 0
        kstrT = 0.
        kdotT = 0.
        alphaT = 0.        
        gruposgc(:) = 0
        
        call GC_Adapt(compuesto,1, NG, nu, q, Tstr, gstr, gdot, gddot,kstrT,kdotT,alphaT,gruposgc)
        TcGC(1) = Tc
        PcGC(1) = Pc
        BpGC(1) = Bp

     !   call boilingPoint (1,NG,0,TcGC(1),PcGC(1),BpGC(1),dcGC(1),nu,q,Tstr,gstr,gdot,gddot,kstr,kdot,alpha,Tb,fallo)
     !   if(fallo)Bp = 5000
     !   if(.not.isnan(Tb(1)))Bp = Tb(1)
        
    endif

    if(present(MWT))MWT = MW
    if(present(BPT))BPT = BP
    if(present(TCT))TCT = TC
    if(present(PCT)) PCT = PC
    if(present(DENST))DENST = MW/VC
    if(present(HVAPT))HVAPT = 1.987*TC*(BPTC*LOG(PC)/(1.0 - BPTC))
    if(present(PvapT))then
        PvapT = VaporPressure (Top,BP)
        if(Model /= 3)then
            if(present(ViscT))ViscT = Viscosity(compuesto,ICOMP,Top,BP,PvapT)
        endif
    endif
    if(present(FuncGroup).and.(icomp==1.or.icomp==3))FuncGroup = Grupos_Funcionales (compuesto,icomp)
  endsubroutine Calc_Prop_Pure
  
!===================================================================
  integer function groups_number (compuesto)
!------------------------------------------------------------------- 
!    Esta función devuelve la cantidad de grupos distintos
!    de la estructura ingresada en el vector "compuesto"
!-------------------------------------------------------------------
    implicit none
    integer,intent(in)::compuesto(10,2)
    integer::i
!Sentencias---------------------------------------------------------
    i=1
    do while(compuesto(i,1)/=0)
        i=i+1
    enddo
    groups_number = i-1
  endfunction groups_number
  
!===================================================================
  real*8 function MolecularWeight (compuesto)
!------------------------------------------------------------------- 
!    Esta función devuelve el peso molecular
!    de la estructura ingresada en el vector "compuesto"
!-------------------------------------------------------------------    
    use SubGrupos
    implicit none
!Parámetros
!Variables de ENTRADA
    integer,intent(in)::compuesto(10,2)
!Variables INTERNAS
    integer::i
    real*8::MWT    
!Sentencias--------------------------------------------------------- 
        MWT=0.0 
        i=1
        do while (compuesto(i,1)/=0)
            MWT = MWT + Obtain_MW(compuesto(i,1))*compuesto(i,2)
            i=i+1
        enddo
        MolecularWeight = MWT
  endfunction MolecularWeight    
  
!===================================================================
  real*8 function CriticalVolume(compuesto,icomp)
!------------------------------------------------------------------- 
!    Esta función devuelve el volumen crítico
!    de la estructura ingresada en el vector "compuesto".
!    La variable de entrada "icomp" indica el tipo de estructura
!    almacenada en "compuesto" (alifática, aromática, etc.)
!-------------------------------------------------------------------  
    use SubGrupos
    implicit none
!Parámetros
!Variables de ENTRADA
    integer,intent(in)::compuesto(10,2)
    integer::icomp
!Variables INTERNAS
    integer::i
    real*8::SUMV
!Sentencias---------------------------------------------------------    
    SUMV = 0.0
    i=1
    do while (compuesto(I,1).ne.0)
        SUMV = Obtain_dV(compuesto(i,1))*compuesto(I,2) + SUMV
        i=i+1
    enddo
    if(ICOMP.EQ.1)then !Aromatic
        SUMV = SUMV + 12
    elseif(ICOMP==2 .or. (i==2.and.compuesto(1,2)==1))then
        CriticalVolume=SUMV
        return
    elseif(ICOMP==4)then !Cyclic
        SUMV = SUMV - 7
    endif
    CriticalVolume =(SUMV/.285)**(1/1.048) !(Ec. 4.11)
  endfunction CriticalVolume

!===================================================================
  real*8 function CriticalPressure(compuesto,PMOLT)
!------------------------------------------------------------------- 
!    Esta función devuelve la presión crítica
!    de la estructura ingresada en el vector "compuesto".
!    Puesto que la ecuación para el cálculo de la presión crítica
!    requiere el valor del peso molecular del compuesto, esta 
!    información debe ingresarse en la variable de entrada "PMOLT"
!-------------------------------------------------------------------  
    use SubGrupos
    implicit none
!Variables de ENTRADA
    integer,intent(in)::compuesto(10,2)
    real*8,intent(in)::PMOLT
!Variables internas
    integer::i
    real*8::SUMPC
!Sentencias---------------------------------------------------------   
    SUMPC = 0.0
    i=1
    do while (compuesto(I,1).ne.0)
       SUMPC = Obtain_dPC(compuesto(I,1))*compuesto(I,2) + SUMPC
       i=i+1
    enddo
    CriticalPressure = PMOLT/(SUMPC + .34)**(2.0) !(Ec. 4.10)
  endfunction CriticalPressure
  
!===================================================================
  real*8 function TbTcRatio(compuesto,icomp,PMOLT)
!------------------------------------------------------------------- 
!    Esta función devuelve la presión crítica
!    de la estructura ingresada en el vector "compuesto".
!    La variable de entrada "icomp" indica el tipo de estructura
!    almacenada en "compuesto" (alifática, aromática, etc.)
!    Puesto que la ecuación para el cálculo de la presión crítica
!    requiere el valor del peso molecular del compuesto, esta 
!    información debe ingresarse en la variable de entrada "PMOLT"
!-------------------------------------------------------------------  
    use SubGrupos
    implicit none
!Variables de ENTRADA
    integer,intent(in)::compuesto(10,2),icomp
    real*8,intent(in)::PMOLT
!Variables INTERNAS
    integer::i
    real*8::a,b,c,d,SUMTC,SUMTCA
!Sentencias---------------------------------------------------------    
    a = 0.0268
    if(icomp==4) b = 3.0 !Cyclic compounds
    if(icomp/=4) b = 1.0 !Aliphatic compounds
    c = 1.3769
    d = 0.1219
    SUMTC = 0.0
    SUMTCA = 0.0
    i=1
    do while (compuesto(I,1).ne.0)
       SUMTC = Obtain_dTC(compuesto(I,1))*compuesto(I,2) + SUMTC
       SUMTCA= Obtain_dTCa(compuesto(I,1))*compuesto(I,2)+SUMTCA
       i=i+1
    enddo
    if(ICOMP.EQ.2.or. (i==2.and.compuesto(1,2)==1))then !Single group solvent
       TbTcRatio = SUMTC
    else
       TbTcRatio = 1-(PMOLT**a/(b+SUMTC))+((SUMTCA)/((PMOLT**c)*(d)))
    endif
  endfunction TbTcRatio
  
!===================================================================  
  real*8 function NormalBoilingPoint (Vc,Pc,TbTc)
!------------------------------------------------------------------- 
!    Esta función devuelve la temperatura normal de ebullición
!    a partir de los valores de volumen crítico (Vc), presión 
!    crítica (Pc) y de la relación Tb/Tc (TbTc)
!-------------------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    real*8,intent(in)::Vc,Pc,TbTc    
!Sentencias---------------------------------------------------------      

    NormalBoilingPoint = Vc*Pc*TbTc*(0.002876*Log(Pc)*TbTc/(1.0-TbTc)+0.02603)

  endfunction NormalBoilingPoint


!===================================================================  
  real*8 function VaporPressure (T,BP)
!------------------------------------------------------------------- 
!    Esta función devuelve la temperatura normal de ebullición
!    a partir de los valores de volumen crítico (Vc), presión 
!    crítica (Pc) y de la relación Tb/Tc (TbTc)
!-------------------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    real*8,intent(in)::T,BP
!Variables INTERNAS
    real*8::REL,AA,BB,CC
!Sentencias---------------------------------------------------------      

    REL  = T/BP
    AA   = 4.5398 + 1.0309*DLOG(BP)
    BB   = 1- ((3-2*REL)**0.19)/REL
	CC   = 0.38 * ((3-2*REL)**(-0.81)) * DLOG(REL)
	VaporPressure = DEXP(AA*(BB-CC))

  endfunction VaporPressure

!===================================================================  
  integer function Grupos_Funcionales (compuesto,icomp)
!------------------------------------------------------------------- 
!    Esta función devuelve la cantidad de grupos funcionales (en 
!    caso de que icomp==3) o la cantidad de carbonos sustituídos (en
!    caso de que icomp==1) de la estructura guardada en "compuesto"
!-------------------------------------------------------------------  
    implicit none
    integer,parameter::NMODEL=3
!Variables de ENTRADA
    integer,dimension(DiffStructGroups,2),intent(in)::compuesto
    integer,intent(in)::icomp
!Variables INTERNAS/COMMONS
    integer::j,GrupFunc
    common/GRUPAR/IACH,IACCH2,IACCL,IACNH2,IACNO2
        integer,dimension(NMODEL)::IACH,IACCH2,IACCL,IACNH2,IACNO2
    COMMON/PAREQ/ipareq
        integer::ipareq
!Sentencias---------------------------------------------------------  
    GrupFunc = 0
    if(icomp==1)then
        j=1
        do while (compuesto(j,1)/=0)
            if (compuesto(j,1) == IACH(ipareq)) GrupFunc = 6 - compuesto(j,2)
            j=j+1
        enddo    
    elseif(icomp==3.or.icomp==4)then
        j=1
        do while (compuesto(j,1)/=0)
			if(compuesto(j,1) > 4) GrupFunc = GrupFunc + compuesto(j,2)
            j=j+1
        enddo
    endif
    Grupos_Funcionales = GrupFunc
    endfunction Grupos_Funcionales


!=================================================================
  real*8 function Viscosity(MS,ICOMP,T,BPOINT,Pvap)
!C
!C----------------------------------------------------------------
!C      Esta rutina calcula la viscosidad de liquidos puros por
!C      contribucion grupal.El metodo se basa en el trabajo:
!C      "A new group contribution method predicting viscosity of 
!C      organic liquid" de S.R.S. Sastri & K.K. Kao,
!C      The Chemical Engineering Journal, 50(1992), 9-25.
!C      Se adaptan las contribuciones grupales de este metodo a 
!C      correspondientes al modelo UNIFAC.
!C----------------------------------------------------------------
   
   !   USE PROPI
   use SubGrupos
    PARAMETER(NSCM=10,NCOM=3)
	IMPLICIT real*8(A-H,O-Z)
	integer::ms(10,2)
	DIMENSION DELVI(NMG),DELN(NMG),DELVC(NMG),PMG(NMG),&
     		   DELNC(NMG),DELTC(NMG),DELPC(NMG),DELV(NMG),DELTCA(NMG)
      COMMON/PUNSUB/NPUNT(NMG),NGRUP(NMG),NCANT
	COMMON/CONT/PMG,DELTC,DELV,DELPC,DELVI,DELN,DELTCA
	!COMMON/US/MS(NCOM,NSCM,2)
     
     
     
	SUMVI=0.0
      SUMN=0.0
      SUMVC=0.0
	SUMNC=0.0
!C	
      DELVC(:)=0.0
      DELNC(:)=0.0      
!C
      INDAL1=0
      INDAL2=0
      INDAL3=0
      INDAL4=0
      INDAR=0
      INDETI=0
      INDNETI=0
      INDOH=0
      INDET=0
      INDNH=0
      INDKET=0
      INDNO2=0
      INDFOL=0
      INDALH=0
      INDNH2=0
      INDACL=0
      INDAN2=0
      NCARBAL=0
      INDNALQUI=0
      INDAC=0
!C	
      i=1
      DO while(MS(i,1)/=0)
	  IF(MS(I,1).EQ.1)THEN 
	      INDAL1=1
	      NCARBAL=NCARBAL+MS(I,2)
	  ENDIF
	  IF(MS(I,1).EQ.2)THEN
	      INDAL2=1
	      NCARBAL=NCARBAL+MS(I,2)
	  ENDIF
	  IF(MS(I,1).EQ.3)INDAL3=1
	  IF(MS(I,1).EQ.4)INDAL4=1
	  IF(MS(I,1).EQ.5)INDNETI=1
	  IF(MS(I,1).EQ.66)INDNALQUI=1
	  IF(MS(I,1).EQ.10)INDAR=1
	  IF(MS(I,1).EQ.43)INDAC=1

	  IF(MS(I,1).EQ.6)INDETI=1
	  IF(MS(I,1).EQ.15)INDOH=1
	  IF(MS(I,1).EQ.25.OR.MS(I,1).EQ.26.OR.MS(I,1).EQ.27)INDET=1
	  IF(MS(I,1).EQ.32.OR.MS(I,1).EQ.33.OR.MS(I,1).EQ.34)INDNH=1
	  IF(MS(I,1).EQ.19.OR.MS(I,1).EQ.20)INDKET=1
	  IF(MS(I,1).EQ.58)INDNO2=1
	  IF(MS(I,1).EQ.18)INDFOL=1
	  IF(MS(I,1).EQ.21)INDALH=1
	  IF(MS(I,1).EQ.29.OR.MS(I,1).EQ.30.OR.MS(I,1).EQ.31)INDNH2=1
	  IF(MS(I,1).EQ.54)INDACL=1
	  IF(MS(I,1).EQ.37)INDAN2=1
	  i=i+1
      enddo
      i=1
      DO while(MS(i,1)/=0)
        TDELVI=Obtain_dVi(MS(I,1))
        TDELN=Obtain_dN(MS(I,1))
        TDELVC=0.0 !DELVC(NPUNT(MS(I,1)))!correcciones
        TDELNC=0.0 !DELNC(NPUNT(MS(I,1)))!correcciones
        IF(NUMGRU.EQ.2.AND.INDAL1.EQ.1.AND.INDAL2.EQ.1.AND.& !Para compuestos alfifaticos con C>8
         (MS(1,2)+MS(2,2)).GT.8.AND.MS(I,1).EQ.1)  TDELN=0.05
        IF(NUMGRU.EQ.3.AND.INDAL1.EQ.1.AND.INDAL2.EQ.1.AND.INDNETI.EQ.1& !Para n-alquenos con C>8
         .AND.(MS(1,2)+MS(2,2)+MS(3,2)).GT.7.AND.MS(I,1).EQ.5)TDELN=0.05
        IF(NUMGRU.EQ.3.AND.INDAL1.EQ.1.AND.INDAL2.EQ.1.AND.& !Para n-alquinos con C>8
         INDNALQUI.EQ.1.AND.(MS(1,2)+MS(2,2)+MS(3,2)).GT.7.AND.MS(I,1).EQ.66)TDELN=0.05
!C.......PARA ALCOHOLES ALIFÁTICOS
        IF(NUMGRU.EQ.3.AND.INDAL1.EQ.1.AND.INDAL2.EQ.1.AND.& !Para n-alcoholes con c<11
         INDOH.EQ.1.AND.(MS(1,2)+MS(2,2)+MS(3,2)).LE.&
         11.AND.MS(I,1).EQ.15) TDELVI=0.615-0.092*NCARBAL+0.004*(NCARBAL**2)-10**(-0.58*NCARBAL)
!        IF(NUMGRU.EQ.3.AND.INDAL1.EQ.1.AND.INDAL2.EQ.1.AND. !Para n-alcoholes con c<11
!     *    INDOH.EQ.1.AND.(MS(1,2)+MS(2,2)+MS(3,2)).GT.11.AND.MS(I,1)
!     *    .EQ.15) TDELVI=0.2
        IF(NUMGRU.EQ.3.AND.INDAL1.EQ.1.AND.INDAL2.EQ.1.AND.& !Para n-alcoholes con c<11
         INDOH.EQ.1.AND.(MS(1,2)+MS(2,2)+MS(3,2)).LE.13.AND.MS(I,1).EQ.15) TDELN= 0.3
!        IF(NUMGRU.EQ.3.AND.INDAL1.EQ.1.AND.INDAL2.EQ.1.AND. !Para n-alcoholes con c<11
!     *    INDOH.EQ.1.AND.(MS(1,2)+MS(2,2)+MS(3,2)).GT.13.AND.MS(I,1)
!     *    .EQ.15) TDELN= 0.15
!C......Para ácidos alifáticos
        IF((NUMGRU.EQ.3.OR.NUMGRU.EQ.2).AND.INDAL1.EQ.1& !Para n-alcoholes con c<11
     .AND.INDAC.EQ.1.AND.(NCARBAL.EQ.1.OR.NCARBAL.EQ.2.OR.NCARBAL.EQ.3).AND.MS(I,1).EQ.43) THEN
        TDELVI=0.220
        TDELN= 0.05   
      endif
     
        !if (MS(I,1).eq.10) !call pausa
        IF(INDAL3.EQ.1.OR.INDAL4.EQ.1)THEN !si existe ramificación
	      IF((INDAL3.EQ.1).AND.(MS(I,1).EQ.4))THEN
		        TDELN=0.0
            ENDIF    
	      IF((MS(I,1).GE.9).AND.(MS(I,1).LT.45).AND.(MS(I,1).NE.32).OR.(MS(I,1).NE.33))THEN
		        TDELN=0.0
            END IF
        ENDIF
	  IF(MS(I,1).EQ.10)THEN
	      IF(MS(I,2).EQ.6)THEN
		        TDELN=0.1
		    ELSE IF (MS(I,2).LE.4)THEN
			    TDELN=0.025
			    TDELVC=0.070
		    ELSE IF((INDAL3.EQ.1.AND.INDAL4.EQ.1).AND.(MS(I,2).LE.5))THEN
		        TDELN=0.025
		    ELSE IF((INDAL3.EQ.1).AND.(MS(I,2)).EQ.5)THEN
			    TDELN=0.05
		    ENDIF
        ENDIF
        IF(INDOH.EQ.1)THEN
            IF(INDET.EQ.1)THEN
		        TDELN=0.1
            ELSE IF (INDNH.EQ.1)THEN
		        TDELN=0.3
            ELSE IF (INDKET.EQ.1)THEN
			    TDELN=0.125
            ENDIF
        ENDIF
	  IF(MS(I,1).EQ.18)THEN
            IF(INDFOL.EQ.1.OR.INDALH.EQ.1)THEN
		        TDELN=0.125
            ENDIF
        ENDIF
	  IF(MS(I,1).EQ.19.OR.MS(I,1).EQ.20)THEN
            IF(INDNH2.EQ.1)THEN
		        TDELVC=0.1
            ENDIF
        ENDIF
	  IF(MS(I,1).EQ.31.AND.(INDKET.EQ.1))THEN
		   TDELVC=0.1
        ENDIF
	  IF(MS(I,1).EQ.22.OR.MS(I,1).EQ.23.OR.MS(I,1).EQ.24)THEN
		   TDELN=0.05*MS(I,2)
        ENDIF     
	  IF(MS(I,1).EQ.58)THEN
	      IF(MS(I,2).GE.2)THEN
		        TDELN=0.1
            ENDIF
            IF(INDFOL.EQ.1)THEN
		        TDELN=0.125
            ENDIF
        ENDIF
	  IF(MS(I,1).EQ.18)THEN
            IF(INDNO2.EQ.1)THEN
		        TDELN=0.125
            END IF
	  ENDIF
	  IF(MS(I,1).EQ.54)THEN
	      IF(INDAN2.EQ.1)TDELNC=-0.1
	      IF(INDFOL.EQ.1)TDELNC=-0.05
	      IF(INDNO2.EQ.1)TDELNC=-0.05
        ENDIF
	  IF(MS(I,1).EQ.26)THEN
	      IF(MS(I,2).GT.1)THEN
		        TDELN=0.05*MS(I,2)
            ENDIF
	  ENDIF
	  IF(MS(I,1).EQ.64)THEN
	      IF(INDAL1.EQ.1.OR.INDAL2.EQ.1.OR.INDAL3.EQ.1)THEN
		        TDELVI=0.26
	      ENDIF
	      IF(INDAR.NE.0)THEN
		        TDELN=0.025
	      ENDIF
	  ENDIF
	  IF(MS(I,1).EQ.65)THEN
	      IF(INDAL1.EQ.1.OR.INDETI.EQ.1)THEN
		        TDELVI=0.24
	      ELSE
		        TDELVI=0.21
	      ENDIF
	      IF(INDAR.NE.0)THEN
		        TDELN=0.025
	      ENDIF
	  ENDIF
	  SUMVI = SUMVI + TDELVI*MS(I,2)
	  SUMN  = SUMN  + TDELN
	  SUMVC = SUMVC + TDELVC
	  SUMNC = SUMNC + TDELNC
	  i=i+1
      enddo

      IF(ICOMP.EQ.2)THEN !Single group solvent
        Viscosity=SUMVI              
        return
      ENDIF

      VIS0 = SUMVI + SUMVC
      ENE  = 0.2 + SUMN + SUMNC
!      if(sumvc.gt.0.or.sumnc.gt.0) then
!        write(*,*)"sumvc=", sumvc
!        write(*,*)"sumnc=", sumnc
!       ! call pausa
!      endif

!C
!C========VISCOSIDAD DINAMICA (mPa/seg)_
!C
	 Viscosity = VIS0*Pvap**(-ENE)
       endfunction Viscosity
       

       
       
endmodule