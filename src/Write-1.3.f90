subroutine write_results (mop,idev,ifam,nsol)
!--------------------------------------------------------------
!      
!--------------------------------------------------------------
  use SubGrupos
  use StructuresDesign
  use constantes
  implicit none
!INTERFACES
  interface
    subroutine write_compound(idev,mop,tazeo,recorre)
        use StructuresDesign
        type(FinalStructure),pointer,intent(in)::recorre
        integer,intent(in)::idev,mop
        real*8,intent(in)::tazeo
    endsubroutine write_compound
  endinterface
!Varibles de ENTRADA
  integer,intent(in)::mop,idev,ifam,nsol
!Variables INTERNAS/COMMONS
  character(len=70),dimension(4)::Titulo
  character::opt
  integer::CompXPag,i,j,k,iso,nGruposFunc,i1
  type(FinalStructure),pointer::recorre
  common/NOM/FS
      character(len=8),dimension(NMG)::FS
  common/PUNSUB/NPUNT,NGRUP,NCANT
      integer,dimension(NMG)::NPUNT,NGRUP
      integer::NCANT
  common/EXTDIS/PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
      real*8::psat1,psat2,tazeo,x1azeo,error

!SENTENCIAS
    if(idev==2)then
        CompXPag = COMPUESTOSXPAGINA_ARCHIVO
    else
        CompXPag = COMPUESTOSXPAGINA_CONSOLA
    endif
!Inicializa el vector TITULO
    if(ifam==1)then
        Titulo(1)= '**** AROMATIC RINGS WITHOUT SUBSTITUTED CARBON ****'
        Titulo(2)= '**** AROMATIC RINGS WITH ONE SUBSTITUTED CARBON ****'
        Titulo(3)= '**** AROMATIC RINGS WITH TWO ONE SUBSTITUTED CARBONS ****'
        Titulo(4)= '**** AROMATIC RINGS WITH THREE ONE SUBSTITUTED CARBONS ****'
    elseif(ifam==2)then
        Titulo(1)= ''
        Titulo(2)= ''
        Titulo(3)= ''
        Titulo(4)= ''
    elseif(ifam==3)then
        Titulo(1)= '**** SOLVENTS WITHOUT FUNCTIONAL GROUPS *****'
        Titulo(2)= '**** SOLVENTS WITH ONE FUNCTIONAL GROUP *****'
        Titulo(3)= '**** SOLVENTS WITH TWO FUNCTIONAL GROUPS *****'
        Titulo(4)= '**** SOLVENTS WITH MORE THAT TWO FUNCTIONAL GROUPS *****'
    elseif(ifam==4)then
        Titulo(1)= ''
        Titulo(2)= ''
        Titulo(3)= ''
        Titulo(4)= ''
    endif
    recorre=>FMSs
    k=1
    do while (k<=nsol)
        !Escribir encabezado
        if(mop == 0)then
            if(idev == 2)then
                write(idev,"")
            else
                write(idev,"")
            endif
        elseif(mop == 1)then
            if(idev == 2)then
                write(idev,"(/,' - -     SELECT.  ','SOL.POW.  ','SOL.LOS.    ','M.W.    ','B.P.   ','DENS.RAT.  ','VAP.HEAT   ','DIS.COEF.  ','VISC.SOLV.  ','Pow     ',' Conversion'&
                             /,' - -     [WT]     ','[WT %]    ','[WT %]      ',8X        ,'[K]    ',11X,          '[CAL/GMOL] ','[WT]       ','[mPa*seg]')")
                write(idev,"(110('*'))")
            else                                                                                                  
                write(idev,"(/,1X,' DIS.COEF.',1X,'   SELECT.',1X,'  SOL.POW.',1X,' SOL.LOST.',1X,'      B.P.',1X,'       Pow',1X,'      M.W.',/,&
                                  '    [ WT ]',1X,'    [ WT ]',1X,'  [ WT % ]',1X,'  [ WT % ]',1X,'       [K]',1X)")
                write(idev,"(80('*'))")
            endif
        elseif(mop == 2)then
            if(tazeo /= 1)then
                write(idev,"(' REL.VOL.',2X,'SOL.POW.',2X,' M.W. ',2X,'B.P. ',2X,&
                    'VAP.HEAT',4X,' MINX3  ',2X,'  P.I. (REL.VOL./M.W.)',&
                    'VISC. SOLV.',/,' ',&
                    10X,'[ WT % ]',10X,'[K]',3X,'[CAL/GMOL]',3X,'[MOL %]',9X,&
                    '(MINX3 M.W.))',4x,'[mPa*seg]')")
            else
                write(idev,"(' REL.VOL.',2X,'SOL.POW.',2X,' M.W. ',2X,'B.P. ',2X,&
                    'VAP.HEAT',11X,'  P.I. ',3x,'VISC.SOLV.',/,' ',10X,&
                    '[ WT % ]',10X,'[K]',3X,'[CAL/GMOL]',5X,'(REL.VOL./M.W.)',&
                    '[mPa*seg]')")
            endif       
        elseif(mop == 3)then
            if(idev == 2)then
                write(idev,"")
            else
                write(idev,"")
            endif
        elseif(mop == 4)then !PROVISORIO para reacciones!!!
                write(idev,"(/,1X,'      Conversion','     B.P.',1X,'      M.W.',/,'                      [K]',1X)")
                write(idev,"(80('*'))")            
        endif
        nGruposFunc = recorre%FunctionalGroups + 1
        if (nGruposFunc > 4) nGruposFunc = 4
        write(idev,"(A70)") Titulo(nGruposFunc) 
        i=1
        do while (i <= COMPXPAG) !bucle para imprimir páginas
            call write_compound(idev,mop,tazeo,recorre)
            if(.not.associated(recorre%next))then
                recorre => recorre%next
                exit
            endif
            if(recorre%FunctionalGroups /= recorre%next%FunctionalGroups)then
                recorre => recorre%next
                exit
            endif
            recorre => recorre%next
            i=i+1            
        enddo
        if(.not.associated(recorre))exit
        if (idev==6)then !opciones para impresión por consola
            write (IDEV,160)
            read (5,'(A)') OPT
            I1 = INDEX ('1234',OPT)
            if (i1==2)then
                k=0
                recorre=>FMSs
            elseif(i1==3)then
                exit
            endif
        endif 
        k=k+1
    enddo

!FORMATS
309	FORMAT (2X,I4,X,A27,A27,X,F5.1)
5020  format &
     (' The  list  of  solvents  is shown according to the following',&
      ' nomenclature:',//,&
      '                UNIFAC solvent structure.     ',/,&
      '                           |',/,&
      '               ---------------------------',/,&
      '      13**     (CH2) 1(CH2COO) 1   (CH3) 2 ',/,&
      '      ----',/,&
      '       |')
5030  FORMAT (&
      '    Solvent ranking based on distribution coefficient.',///,&
      ' SELECT.   : solvent selectivity [ wt ].',/,&
      ' SOL.POW.  : solvent power [ wt % ].',/,&
      ' SOL.LOS.  : solvent loss [ wt % ].',/,&
      ' M.W.      : solvent molecular weight.',/,&
      ' B.P.      : solvent boiling point [ K ].',/,&
      ' B.P.DIFF. : solvent and solute boiling points difference',&
      ' [ K ].',/,&
      ' DENS.RAT. : solvent density ratio.',/,&
      ' VAP.HEAT  : solvent latent heat of vaporization [cal/molgr].'&
      ,/,&
      ' DIST.COEF.: solvent distribution coefficient [ wt ].',/,&
      ' VISC.SOLV.: viscosidad del solvente [mPa*seg].' )
5040  FORMAT (&
      '    Solvent ranking based on PI index.',///,&
      ' REL.VOL.  : relative volatility [ wt ].',/,&
      ' SOL.POW.  : solvent power [ wt % ].',/,&
      ' M.W.      : solvent molecular weight.',/,&
      ' B.P.      : solvent boiling point [ K ].',/,&
      ' HEAT VAP. : solvent latent heat of vaporization [cal/molgr].'&
      ,/,&
      ' MINX3.    : solvent concentration to break the azeotrope [mol',&
      '%].',/,&
      ' P.I.      : performance index.',/,&
      ' VISC.SOLV.: viscosidad del solvente [mPa*seg].' )
5050  FORMAT (&
      '    Solvent ranking based on PI index.',///,&
      ' REL.VOL.  : relative volatility [ wt ].',/,&
      ' SOL.POW.  : solvent power [ wt % ].',/,&
      ' M.W.      : solvent molecular weight.',/,&
      ' B.P.      : solvent boiling point [ K ].',/,&
      ' HEAT VAP. : solvent latent heat of vaporization [cal/molgr].'&
      ,/,&
      ' P.I.      : performance index.',/,&
      ' VISC.SOLV.: viscosidad del solvente [mPa*seg].' )

160    FORMAT (' ',/,3X,'CONTINUE: 1',5X,'BEGIN AGAIN: 2',5X,&
              'EXIT: 3',5X,'HELP: 4',5X,'SEARCH: 5',5X,'> ',$)


1000   RETURN

endsubroutine write_results

!=======================================================
subroutine write_compound(idev,mop,tazeo,recorre)
  use SubGrupos
  use StructuresDesign
  implicit none
!INTERFACES
  interface
    subroutine write_isomers(idev,recorre)
        use StructuresDesign
        integer,intent(in)::idev
        type(FinalStructure),pointer,intent(in)::recorre 
    endsubroutine write_isomers
  endinterface
  interface
    subroutine write_performance(idev,recorre)
        use StructuresDesign
        type(FinalStructure),pointer,intent(in)::recorre
        integer,intent(in)::idev
    endsubroutine write_performance
  endinterface
!Variables de ENTRADA
  type(FinalStructure),pointer,intent(in)::recorre
  integer,intent(in)::idev,mop
  real*8,intent(in)::tazeo
!Variables INTERNAS
  integer::j

!   write (idev,"(1X, /, I4, '*', 10(1A8, X, 1I2, 3X))") recorre%position,(Obtain_SubGroup_Name(recorre%Formula(j,1)),&
!         recorre%Formula(j,2),j=1,recorre%GroupsNumber)
  if(mop == 1)then ! Propiedades estimadas
    
      if(idev == 6)then
  
        write (idev,"(1X, /, I4, '*', 10(1A8, X, 1I2, 3X))") recorre%position,(Obtain_SubGroup_Name(recorre%Formula(j,1)),&
            recorre%Formula(j,2),j=1,recorre%GroupsNumber)        
        write(idev,"(7(F11.2))")recorre%DistCoefficient,recorre%Selectivity,&
            recorre%SolventPower,recorre%SolventLost,&
            recorre%BoilingPoint,recorre%Kow,recorre%MolecularWeight
        if(recorre%numberIsomers/=0)call write_isomers(idev,recorre)
        if(recorre%SolutesNumber > 1)call write_performance(idev,recorre)
          
      else
 
        write (idev, "(1X, /, I4, ' º ', F8.2, F10.2, F10.3, F8.1, F8.1, F12.1, F11.2, X, G11.5, F11.3, F12.4, F10.2, 10(1A8, X, 1I2, 3X))") & 
       
          recorre%position, recorre%Selectivity, recorre%SolventPower, recorre%SolventLost, recorre%MolecularWeight,     &
          
          recorre%BoilingPoint, recorre%RDER, recorre%Hvap, recorre%DistCoefficient, recorre%Viscosity,                  &
          
          recorre%Kow, recorre%conversion,                                                                               &
          
          (Obtain_SubGroup_Name(recorre%Formula(j,1)), recorre%Formula(j,2),j=1,recorre%GroupsNumber)                   
        
!         write (idev,"(F8.2,F10.2,E10.3,F8.1,F8.1,F12.1&
!             ,F11.2,F11.2,F11.3,F12.4,F10.2)")recorre%Selectivity,recorre%SolventPower,&
!             recorre%SolventLost,recorre%MolecularWeight,recorre%BoilingPoint,&
!             recorre%RDER,recorre%Hvap,recorre%DistCoefficient,recorre%Viscosity,recorre%RelVolatility,recorre%Kow
        if (recorre%numberIsomers /= 0) call write_isomers(idev,recorre) 
        if (recorre%SolutesNumber > 1) call write_performance(idev,recorre)   
          
      endif
      
    elseif(mop == 2)then
      if (tazeo /= 1)then
            write(IDEV,"(' ',F8.2,2X,F8.2,2X,F6.1,2X,F5.1,2X,F8.2,2X,&
              F8.2,2X,F8.4,X,F8.4)")recorre%Selectivity,recorre%SolventPower,&
              recorre%MolecularWeight,recorre%BoilingPoint,recorre%Hvap,&
              recorre%SolventLost,recorre%DistCoefficient,recorre%Viscosity
          if(recorre%numberIsomers/=0)call write_isomers(idev,recorre)
      else 
            write(IDEV,"(' ',F8.2,2X,F8.2,2X,F6.1,2X,F5.1,2X,F8.2,2X,&
              F8.4,x,F8.4)")recorre%Selectivity,recorre%SolventPower,&
              recorre%MolecularWeight,recorre%BoilingPoint,recorre%Hvap,&
              recorre%DistCoefficient,recorre%Viscosity
          if(recorre%numberIsomers/=0)call write_isomers(idev,recorre)    
      endif
    elseif(mop == 4)then !PROVISORIO para reacciones!!!!
        write (idev,"(1X, /, I4, '*', 10(1A8, X, 1I2, 3X))") recorre%position,(Obtain_SubGroup_Name(recorre%Formula(j,1)),&
            recorre%Formula(j,2),j=1,recorre%GroupsNumber)        
        write(idev,"(F11.4,6X, 6(F11.2))")recorre%conversion, recorre%BoilingPoint,recorre%MolecularWeight
        if(recorre%numberIsomers/=0)call write_isomers(idev,recorre)   
    endif

endsubroutine write_compound

!====================================================
subroutine write_isomers (idev,recorre)
    use SubGrupos
    use StructuresDesign
    implicit none
!Variables de ENTRADA
    integer,intent(in)::idev
    type(FinalStructure),pointer,intent(in)::recorre
!Variables INTERNAS
    type(isomers),pointer::recorreIso    
    integer::iso
    
!SENTENCIAS   
     recorreIso => recorre%nextisomer
     do iso=1,recorre%numberIsomers ! propiedades prop.pdb
         if(iso==1)write(idev,"(2X,'PROP.PDB contains the following isomers for this ',&
             'solvent molecular strucure:',/,2X,'R.N.',x,'Name',23x,'Formula',20x,'B.P.(K)',3x)")
         write(idev,"(2X,I4,X,A27,A27,X,F5.1,4X,F8.2)	")recorreIso%index,&
     	        recorreIso%name,recorreIso%FormChem,recorreIso%BoilingPoint !,VRELSOLV(ISOL(K,1),N)              
         recorreIso => recorreIso%next
     enddo
endsubroutine write_isomers

!====================================================
subroutine write_performance(idev,recorre)
    use SubGrupos
    use StructuresDesign
    implicit none
!Variables de ENTRADA    
    type(FinalStructure),pointer,intent(in)::recorre
    integer,intent(in)::idev
!Variables INTERNAS
    type(CalculatedPerformance),pointer::recorrePerformance
    integer::i

!SENTENCIAS     
    write(idev,"(/,3X,75('-'))")
    write(idev,"('   COMPOUND    SELECT.[WT]    SOL.POW.[WT %]      DIS.COEF.[WT]      BPDIFF[K]')")  
    write(idev,"(3X,75('-'))")
    recorrePerformance => recorre%Performance
    i=1
    do while(associated(recorrePerformance))
        write(idev,"(I4,F22.2,F18.2,F19.3,F15.2)") i,recorrePerformance%Selectivity,&
             recorrePerformance%SolventPower,recorrePerformance%DistCoefficient,recorrePerformance%DifferenceBP
        recorrePerformance => recorrePerformance%next
        i=i+1
    enddo
    write(idev,"(3X,75('-'))")
endsubroutine write_performance
!====================================================
subroutine Mostrar_Componente(ISF)

PARAMETER (NMG=150,NCOM=3,NSCM=10,NA=150)
implicit real*8 (A-H,O-Z)
CHARACTER*8 FS(NA)
COMMON/US/MS(NCOM,NSCM,2)
COMMON/PUNSUB/NPUNT(NA),NGRUP(NA),NCANT
COMMON/NOM/FS
      WRITE (6,100)(FS(NPUNT(MS(3,I,1))),MS(3,I,2),I=1,ISF) ! muestra nombre y cantidad
      WRITE (2,100)(FS(NPUNT(MS(3,I,1))),MS(3,I,2),I=1,ISF) ! muestra nombre y cantidad
 100  FORMAT (1X,10(1A8,1I2))
endsubroutine Mostrar_Componente


!==================================================================================
SUBROUTINE ESCLIS (IDEV,MolecularDesign,MSOL,NGDV,NGSDV,MDV,&
     		       NGSV,NGSSV,MSV,A1,A2,A3,B1,B2,B3,TAZEO,X1AZEO,IFAM)
     		       
!----------------------------------------------------------------------------------
!		This subroutine writes the information about the case into the
!		output file (nomb).mdo (when idev=2) or screen (idev=6)
!----------------------------------------------------------------------------------

    use CONSTANTES
    use Input
    use PropertiesData
    use SubGrupos
    use StructuresDesign
    implicit none
    integer,parameter::NA=150,NSCM=10,NSVA=20

!Variables de ENTRADA
    integer,intent(in)::    idev,msol,ifam,mdv,msv
    integer,intent(in)::    NGDV(0:NA),NGSDV(NMG,5),NGSV(NSVA),NGSSV(NMG,5)
    real*8,intent(in)::     A1,A2,A3,B1,B2,B3,TAZEO,X1AZEO
    logical,intent(in)::    MolecularDesign
!Varaibles de SALIDA

!Variables INTERNAS
    integer::               i,j,l
    integer::               mop,nc
    real*8::                Top
    character*70::          titl
    character*1::           SALTO
    logical::               INV
    type(Compound),pointer::recorreSolute
    type(FinalStructure),pointer::ptr_FMS
!...Commons
    COMMON/PUNSUB/NPUNT,NGRUP,NCANT
    integer::NPUNT(NA),NGRUP(NA),ncant
    common/US/MS
    integer,dimension(NCOM,NSCM,2)::MS
    common/EVAL/INV	    

    mop = InputProblem01%mop
    Top = InputProblem01%T
    Titl = InputProblem01%ProblemTitle
    SALTO = CHAR(12)

    if(inv)then
        write(idev,362)
    else
        write(idev,361)
    endif
    IF (MOP.EQ.1) THEN
        WRITE (IDEV,110)
    ELSE
        WRITE (IDEV,120)
    END IF
    WRITE (IDEV,240)
    IF (IFAM.EQ.1) THEN
        WRITE (IDEV,250)
    ELSE IF (IFAM.EQ.2) THEN
        WRITE (IDEV,260)
    ELSE IF (IFAM.EQ.3) THEN
        WRITE (IDEV,270)
    ELSE IF (IFAM.EQ.4) THEN
        WRITE (IDEV,280)
    ELSE IF (IFAM.EQ.5) THEN
        WRITE (IDEV,290)
    END IF

    WRITE (IDEV,367)

    WRITE (IDEV,130)
	WRITE (IDEV,9) InputProblem01%ProblemTitle
	if (mop.eq.1) then
	    write (idev,1099) InputProblem01%T
	end if
    IF (IPAREQ.EQ.1) THEN
        WRITE(IDEV,357)
    ELSE IF (IPAREQ.EQ.2) THEN
        WRITE(IDEV,358)
    ELSE
        WRITE(IDEV,359)
	END IF
	if (mop.eq.1) then
	    if(inv)then
	        write(idev,"(' ',/,' SOLVENT                    : ',$)")
	    else
	        write(idev,"(' ',/,' COMPONENTS TO BE RECOVERED : ',$)")
	    endif    
	    recorreSolute => InputProblem01%MixtureInput%Solutes
	    i=1
	    do while(associated(recorreSolute)) 
	        if(i==1)then
	            write(idev,"(11X,I2,'-',10(A8,I2))")i,(Obtain_SubGroup_Name(recorreSolute%Formula(J,1)),&
	                                                  recorreSolute%Formula(J,2),J=1,recorreSOlute%GroupsNumber)
            else
                write(idev,"(41X,I2,'-',10(A8,I2))")i,(Obtain_SubGroup_Name(recorreSolute%Formula(J,1)),&
	                                                  recorreSolute%Formula(J,2),J=1,recorreSOlute%GroupsNumber)
            endif   
	        recorreSolute => recorreSolute%next
            i=i+1
	    enddo
    else if (mop.eq.2) then
	    write(idev,1362) (Obtain_SubGroup_Name(MS(1,J,1)),MS(1,J,2),J=1,MSOL)
	end if
	if(inv)then
	    write(idev,"(' ',/,' BOILING POINT (K)        : ',$)")
	else
        write(idev,"(' ',/,' BOILING POINT OF SOLUTES (K): ',$)")
	endif
	recorreSolute => InputProblem01%MixtureInput%Solutes
	i=1
	do while(associated(recorreSolute)) 
	    if(i==1)then
            write(IDEV," (10X,I2,'-',F9.2)") i,recorreSolute%BoilingPoint
        else
            write(IDEV," (41X,I2,'-',F9.2)") i,recorreSolute%BoilingPoint
        endif
        recorreSolute => recorreSolute%next
        i=i+1
	enddo	

    IF (MOP.EQ.2) THEN
        WRITE (IDEV,140) A1,A2,A3
        IF (TAZEO.NE.1) THEN
            WRITE (IDEV,160) X1AZEO
            WRITE (IDEV,170) TAZEO
        END IF
    END IF
	if (mop.eq.1) then
	    WRITE (IDEV,"(' ',/,' PRINCIPAL COMPONENT IN THE RAFFINATE: ',10(A8,I2))") &
	          (Obtain_SubGroup_Name(InputProblem01%MixtureInput%PCR%Formula(J,1)),&
	          InputProblem01%MixtureInput%PCR%Formula(J,2),J=1,InputProblem01%MixtureInput%PCR%GroupsNumber)
	          
	else if (mop.eq.2) then
	    write(idev,1363) &
	                (Obtain_SubGroup_Name(InputProblem01%MixtureInput%PCR%Formula(J,1)),&
	                InputProblem01%MixtureInput%PCR%Formula(J,2),J=1,InputProblem01%MixtureInput%PCR%GroupsNumber)
	end if
    IF (MOP.EQ.2) WRITE (IDEV,140) B1,B2,B3

!C

    IF(IDEV.EQ.6)CALL PAUSA 
    if(MolecularDesign)then
        WRITE (IDEV,"(' ',//,' ','** SELECTED GROUPS AND ATTACHMENT ',&
                'COMBINATION PROPERTIES **',/)")
        WRITE (IDEV,"(' ',/,5X,'GROUP  NO.    M    J    K    I    H'/)")
        DO 12 I=1,MDV
            L=NGDV(I)
            WRITE (IDEV,22) Obtain_SubGroup_Name(L),L,NGSDV(I,1),NGSDV(I,2),NGSDV(I,3),NGSDV(I,4),NGSDV(I,5)
  12    CONTINUE
        DO 16 I=1,MSV
            L=NGSV(I)
            WRITE (IDEV,22) Obtain_SubGroup_Name(L),L,NGSSV(I,1),NGSSV(I,2),NGSSV(I,3),NGSSV(I,4),NGSSV(I,5)
  16    CONTINUE
        IF(IDEV.EQ.6)CALL PAUSA
    
        WRITE (IDEV,"(' ',//,' ','** SELECTED PRIMARY SOLVENT PROPERTIES ','**')")
        WRITE (IDEV,"(' ',//,' * FOR THE INTERMEDIATE STRUCTURES:',/)") !
        IF (MOP.EQ.1) THEN
            WRITE (IDEV,"(' SELECTIVITY    SOVENT LOSS')")
            WRITE (IDEV,"('   [ WT ]         [ WT %]',/)")
            WRITE (IDEV,2060) SelIntB(1),SLostIntB(2) !*100
        ELSE
            WRITE (IDEV,180)
            WRITE (IDEV,190)
         !   WRITE (IDEV,220) SCLLI,SOLV*100
        END IF
        WRITE (IDEV,25)
        IF (MOP.EQ.1) THEN
            WRITE (IDEV,"(' SELECTIVITY    SOLVENT POWER    SOLVENT LOSS','    MOLC.WEIGHT    DIST.COEF.')")
            WRITE (IDEV,"('   [ WT ]          [ WT %]          [ WT %]  ','                     [ WT ]',/) ")
            WRITE (IDEV,26) SelB(1),SolPowB(1),SLostB(2),MWB(2),DistCoefB(1)
        ELSE
            WRITE (IDEV,200)
            WRITE (IDEV,210)
          !  WRITE (IDEV,230) SELECT,SOLV*100,DTC2S,PMMA
        END IF
        IF(IDEV.EQ.6)CALL PAUSA 
    endif !MolecularDesign
    
    IF (IDEV.EQ.2) WRITE (IDEV,100) SALTO
    WRITE (IDEV,19)
   ! MODEL = 0
    NC = 2
    ptr_FMS => FMSs
    do while (associated(ptr_FMS))
        ms(3,:,:)=ptr_FMS%formula
        CALL ESCPAR (NC,MODEL,IDEV)
        ptr_FMS => ptr_FMS%next
    enddo
    IF(IDEV.EQ.6)CALL PAUSA 

!...FORMATS
160   FORMAT (1X,'AZEOTROPIC COMPOSITION              : ',F9.2)
170	  FORMAT (1X,'AZEOTROPIC TEMPERATURE (K)          : ',F9.2)
140	  FORMAT (1X,'Antoine vapor pressure eqn.   : ','ln(Pv/mmHg) =',f8.4,' -',f9.2,'/(T/K + ',f7.2,')')
9	  FORMAT (' ',/,' PROBLEM TITLE: ',A70)
1099	  format (' Operation temperature (K) : ',10X,f9.2,/)
19        FORMAT ( ' ***  UNIFAC GROUP  DATA FOR THE KEYS ***',/)
357       FORMAT (' LIQUID-LIQUID UNIFAC PARAMETERS')
358       FORMAT (' LIQUID-VAPOR UNIFAC PARAMETERS')
359       FORMAT (' INFINITE DILUTION UNIFAC PARAMETERS')
367       FORMAT (&
     12X,' *******************************************************',/,&
     12X,' *                                                     *',/,&
     12X,' *    PROGRAM MOLDES: EVALUATION OF SPECIFIC SOLVENTS  *',/,&
     12X,' *                                                     *')
361       FORMAT (&
     12X,' *******************************************************',/,&
     12X,' *                                                     *',/,&
     12X,' *     PROGRAM MOLDES: MOLECULAR DESIGN OF SOLVENTS    *',/,&
     12X,' *                                                     *')
362       FORMAT (&
     12X,' *******************************************************',/,&
     12X,' *                                                     *',/,&
     12X,' *     PROGRAM MOLDES: MOLECULAR DESIGN OF SOLUTES    *',/,&
     12X,' *                                                     *')     
110       FORMAT (12X,' *               LIQUID-LIQUID EXTRACTION              *')
120       FORMAT (12X,' *               EXTRACTIVE  DISTILLATION              *')
240       FORMAT (12X,' *                                                     *')
250       FORMAT (12X,' *                / AROMATIC SOLVENTS /                *')
260       FORMAT (12X,' *              / SINGLE GROUP SOLVENTS /              *')
270       FORMAT (12X,' *                / ALIPHATIC SOLVENTS /               *')
280       FORMAT (12X,' *         / MIXED AROMATIC-ALIPHATIC SOLVENTS /       *')
290       FORMAT (12X,' *                  / CYCLIC SOLVENTS /                *')
130       FORMAT (&
     12X,' *                                                     *',/,&
     12X,' *******************************************************')
!362       FORMAT (' ',/,' COMPONENTS TO BE RECOVERED           : ',10(A8,I2))
1362	  format (' ',/,' LESS VOLATILE COMPONENT             : ',10(A8,I2))
365	  FORMAT (' Boiling point [K]                   : ',F9.2)
!363       FORMAT (' ',/,' PRINCIPAL COMPONENT IN THE RAFFINATE: ',10(A8,I2))
1363	  format (' ',/,' MORE VOLATILE COMPONENT             : ',10(A8,I2))
22        FORMAT (' ',2X,1A8,2X,I2,4X,I1,4X,I1,4X,I1,4X,I1,4X,I1)
!24        FORMAT (' ',//,' ','** SELECTED PRIMARY SOLVENT PROPERTIES ','**')
25        FORMAT (' ',//,' * FOR THE FINAL STRUCTUTRES:',/)
200       FORMAT (' VOLATILITY     SOLVENT POWER    MIN.DIFF.TEM.',&
                 '   MOLC.WEIGHT')
210       FORMAT ('                   [ WT %]           [ K ]  ',/)
!29        FORMAT (' SELECTIVITY    SOLVENT POWER    SOLVENT LOSS','    MOLC.WEIGHT    DIST.COEF.')
!2070      FORMAT ('   [ WT ]          [ WT %]          [ WT %]  ','                     [ WT ]',/)
220       FORMAT (1X,F8.2,7X,F8.2)
230       FORMAT (1X,F8.2,8X,F8.2,9X,F9.3,8X,F8.2)
26        FORMAT (1X,F8.2,8X,F8.2,9X,E9.3,8X,F6.1,5X,F8.2)
2060      FORMAT (1X,F8.2,7X,E9.3)
!27        FORMAT (' ',//,' * FOR THE INTERMEDIATE STRUCTURES:',/)
180       FORMAT ('  VOLATILITY   SOLVENT POWER')
190       FORMAT ('                  [ WT %]',/)
!28        FORMAT (' SELECTIVITY    SOVENT LOSS')
!2080      FORMAT ('   [ WT ]         [ WT %]',/)
100       FORMAT (A)
          RETURN
END


subroutine Groups_Present ()
!-----------------------------------------------------------------------------
!     Esta subrutina presenta por pantalla los subgrupos UNIFAC almacenados
!     en el vector grupos.
!
!     Variables de entrada:
!     --------------------
!                       fs: vector de nombres de identificacion de subgrupos
!                           UNIFAC.  Su dimension es NMG.
!                   grupos: vector  de  numeros  de  identificacion  de  los 
!                           subgrupos a presentar por pantalla.
!                    ngrup: numero de grupos en el vector ngrup.
!
!     La dimension de los vectores esta fijada por los siguientes parametros:
!
     use CONSTANTES
     use SubGrupos
!
!-----------------------------------------------------------------------------
      implicit none
!Variables de ENTRADA

!Variables INTERNAS
      type(Groups),pointer::recorre
      integer::i
!SENTENCIAS
      recorre=>LSubGroups
      do while (associated(recorre))
        do i=1,6
            if(associated(recorre)) then
                write (6,20) recorre%Name,recorre%Number
            else
                exit
            endif
            recorre=>recorre%next
        enddo
        write(*,*)
      enddo
!---- formatos
 20   format (1x,a8,'=',i2,$)
      return
      end