module evaluation
contains
!======================================================================================================
subroutine Evaluate_Mixture(MolecularDesign,FinalStruct,MS,NC,Top,BPSolute,icomp,NSWIP,NSBI,ACEPTABLE,&
           ListPerformance,ent)
!-------------------------------------------------
!   Descripción: Calcula propiedades solventes
!   - Variables de entrada
!       compound(10,2):componente
!       mop:    operación de separación (extracción L-L(1); destilación extractiva(2))
!       icomp:  structure type of the component
!		            Aliphatic		        :0
!		            Aromatic		        :1
!		            Single group solvent	:2
!   - Variables de salida
!-------------------------------------------------
    use PropertiesData
    use grps
    use PureProp
    use StructuresDesign
    use Input
    implicit none
!   VARIABLES DE ENTRADA
    real*8,intent(in)::Top,BPSolute
    logical,intent(in)::MolecularDesign,FinalStruct
    integer,intent(in)::NC
    integer,intent(in)::ent,icomp
!   VARIABLES DE SALIDA
    logical,intent(out)::ACEPTABLE
!   VARIABLES ENTRADA/SALIDA
    integer,intent(inout)::NSWIP,NSBI
    integer,intent(inout)::MS(NCOM,DiffStructGroups,2) !Para que pueda ser modificada por MoldesInvertido
    type(CalculatedPerformance),pointer,intent(inout)::ListPerformance    
!   VARIABLES INTERNAS
    real*8::MWT,BPT,PCT,VCT,TCT,BPTCT,ViscT,DensT,HvapT,DifBPT,KowT,PvapT,&
                      SelIntT,SelT,SolPowT,SLostIntT,SLostT,DistCoefT,VolatIntT,&
                      VolatT,SolPowTmop2,MSMW(NC),BPSolv,PvapSolv,TCSolv
    real*8::RDER !!Checkear esta variable
    integer::i,NMGRUP,FuncGroupsT
    type(CalculatedPerformance),pointer::puntero
    logical::NoAce
!   COMMONS
    common/EVAL/INV
    logical INV

    
!   SENTENCIAS
    mop=InputProblem01%mop
    nullify(puntero)
    allocate(puntero)
   ! write(*,*)ms(3,:,:)
    !nullify(ListPerformance)
    do i=1,NC-1
        MSMW(i) = MolecularWeight (MS(i,:,:))
    enddo
    ACEPTABLE=.False.
    do i=1,size(ms(3,:,2))
        if(ms(3,i,2)<0)return
    enddo
!   Cálculo y evaluación de peso molecular    
    MWT = MolecularWeight (MS(NC,:,:))    
    if (INV) call MoldesInvertido (MS)
!   Cálculo propiedades compuesto puro    
    call Calc_Prop_Pure(MS(3,:,:),icomp,Top,TCT = TCSolv, BPT=BPSolv, PvapT=PvapSolv)
    if(MolecularDesign .and. FinalStruct .and. (Top>TCSolv .or. Top>BPSolv))then
        NOACE = .True.
        return
    endif
!   Cálculo propiedades de mezcla    
    call Solvent_Properties (MS,MOP,Top,SelT,SolPowT,SLostT,DistCoefT,VolatT,NOACE)
    if(model == 3)SLostT = PvapSolv*MWT/100

    if (INV) call MoldesInvertido (MS)
    if (MolecularDesign .and. NOACE) then
        NSWIP=NSWIP+1 !número de solventes sin parámetros de interacción
        return
    endif
!   Evaluación de propiedades de mezcla    
    if(ent==1)NSBI=NSBI+1 !Number of solvents with binary inf.
    if(MOP.EQ.0)then ! Llamado desde Properties
	elseif(MOP.EQ.1) then !LIQUID-LIQUID EXTRACTION
        if(FinalStruct)then !Si es una estructura final
            if(inv)then
                puntero%Selectivity = SelT*MWT/MSMW(2)
                puntero%DistCoefficient = DistCoefT*MSMW(2)/MSMW(1)
                puntero%SolventLost = SLostT*MSMW(1)/MSMW(2)
                puntero%SolventPower = SolPowT*MWT/MSMW(1)
                puntero%DifferenceBP = BPSolute - BPSolv
            else
                puntero%Selectivity = SelT*MSMW(1)/MSMW(2)
                
                if(model/=3)then !UNIFAC o A-UNIFAC
                    puntero%SolventLost = SLostT*MWT/MSMW(2)
                    puntero%DistCoefficient = DistCoefT*MSMW(2)/MWT
                else !GC
                    puntero%SolventLost = SLostT
                    puntero%DistCoefficient = DistCoefT*MSMW(2)/MWT
                endif
                puntero%SolventPower = SolPowT*MSMW(1)/MWT
                puntero%DifferenceBP = BPSolute - BPSolv           
            endif
        else !Si es una estructura intermedia
            if(inv)then
		        puntero%SolventLost = SLostT*MSMW(1)/MSMW(2)
		        puntero%Selectivity = SelT*MWT/MSMW(2)
            else
		        puntero%SolventLost = SLostT*MWT/MSMW(2)
		        puntero%Selectivity = SelT*MSMW(1)/MSMW(2)
            endif
            if(Model == 3)then
                if((puntero%Selectivity >= SelIntB(1)))then
                    ACEPTABLE=.True.
                    return                
                endif
            else
                if((puntero%Selectivity >= SelIntB(1)).OR.(puntero%SolventLost <= SLostIntB(2))) then
                    ACEPTABLE=.True.
                    return
                endif
            endif
        endif
    elseif (MOP.EQ.2) then !EXTRACTIVE DISTILLATION
        if(FinalStruct)then
            puntero%SolventPower = SolPowT*MSMW(1)/MWT
            SolPowTmop2 = SolPowT*MWT/MSMW(1)
        else     !Para estructura intermedia
            if((puntero%RelVolatility >= VolatIntB(1)).OR.(puntero%SolventPower >= SolPowB(1)))then
                ACEPTABLE=.True.
                return
            endif   
        endif
    elseif(MOP==3)then
            puntero%Selectivity = SelT*MSMW(1)/MSMW(2)
    endif
!...Evaluación    
    puntero%SolventLost = puntero%SolventLost * 100
    puntero%SolventPower = puntero%SolventPower * 100
    if (MolecularDesign)call AcceptStruct_Mixture(FinalStruct,puntero,Aceptable)
    if(mop==2) then
        CALL MI_SOLV_BREAK (IPAREQ,Top,SolPowTmop2,ACEPTABLE,SlostT)
        puntero%DistCoefficient = 100*SelT/SLostT*MWT
    endif
    puntero%SoluteFormula = ms(1,:,:)
    call Incorporate_Performance(puntero,ListPerformance)
endsubroutine Evaluate_Mixture

!=======================================================================================
subroutine Evaluate_Pure (FinalStruct,component,Top,ICOMP,NSEW,ACEPTABLE,ListProperties)
!-------------------------------------------------
!   Descripción: Calcula propiedades físicas de compuesto puro
!   - Variables de entrada
!       component(10,2):componente
!       mop:    operación de separación (extracción L-L(1); destilación extractiva(2))
!       icomp:  structure type of the component
!		            Aliphatic		        :0
!		            Aromatic		        :1
!		            Single group solvent	:2
!   - Variables de salida
!-------------------------------------------------
    use PropertiesData
    use grps
    use PureProp
    use StructuresDesign
    use Input
    implicit none
!   VARIABLES DE ENTRADA
    integer,intent(in)::ICOMP !ICOMP debería sacarse cuando se elimine la evaluación de propiedades de compuesto puro
    real*8,intent(in)::Top
    logical,intent(in)::FinalStruct
!   VARIABLES DE SALIDA
    logical,intent(out)::ACEPTABLE
    type(CalculatedProperties),pointer,intent(out)::ListProperties       
!   VARIABLES ENTRADA/SALIDA
    integer,intent(inout)::NSEW
    integer,intent(inout)::component(DiffStructGroups,2) !Para que pueda ser modificada por MoldesInvertido
!   VARIABLES INTERNAS
    real*8::MWT,BPT,PCT,VCT,TCT,BPTCT,ViscT,DensT,HvapT,DifBPT,KowT,PvapT,&
                      SelIntT,SelT,SolPowT,SLostIntT,SLostT,DistCoefT,VolatIntT,&
                      VolatT,SolPowTmop2
    real*8::RDER !!Checkear esta variable
    integer::i,NMGRUP,FuncGroupsT
    logical::NoAce
!   COMMONS
    common/EVAL/INV
    logical INV

!-------------------------------------------------------------------------------------
!
!   SENTENCIAS
    nullify(ListProperties)
    allocate(ListProperties)
    ACEPTABLE=.False.
!   Cálculo y evaluación de peso molecular    
    MWT = MolecularWeight (component)
    if (LMW(2).and.(MWT>MWB(2)))then
        nsew=nsew+1 !número de estructuras con alto peso molecular
        return
    endif 
    if (LMW(1).and.(MWT<MWB(1))) return
    call Calc_Prop_Pure(component,icomp,Top,MWT,TCT,PCT,VCT,BPT,VISCT,DENST,HVAPT,PvapT,FuncGroupsT)
    call OW_Partition_Coefficient(KowT,component)
    ListProperties%MolecularWeight = MWT
    ListProperties%CriticalTemperature = TCT
    ListProperties%CriticalPressure = PCT
    ListProperties%CriticalVolume = VCT
    ListProperties%BoilingPoint = BPT
    ListProperties%Viscosity = VISCT
    ListProperties%Density = DENST
    ListProperties%Hvap = HVAPT
    ListProperties%VapourPressure = PvapT
    ListProperties%FunctionalGroups = FuncGroupsT
    ListProperties%Kow = KowT
    ListProperties%GroupsNumber = groups_number(component)
    
    if(mop==0)then
    elseif(mop==1)then
    elseif(mop==2)then
        DifBPT = BPT-BP1
    elseif(mop==3)then
    endif
    
    if (FinalStruct) call AcceptStruct_Pure(ListProperties,Aceptable)
endsubroutine Evaluate_Pure

!==========================================================================    
    subroutine EvaluationSolvRea (Solv)
!-------------------------------------------------
!   Descripción
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------    
    implicit none
!...VARIABLES DE ENTRADA    
    real*8::solv(10,2)
!...VARIABLES INTERNAS    
    real*8::compounds(20,10,2)
!...COMMONS
    common/reaction/reagents,products
    real*8::reagents(10,10,2),products(10,10,2)
!...Sentencias

    
    endsubroutine
    
!==========================================================================    
    subroutine Constraints_Selection()
!-------------------------------------------------
!   Descripción
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
    !use PropertiesData
    !use blockdatas
    use Input
    use Boundaries
    
    implicit real*8(A-H,O-Z)
!...VARIABLES DE ENTRADA y SALIDA
    !real*8,intent(inout)::Top
!...VARIABLES DE ENTRADA
  !  integer,intent(in)::mop
!...VARIABLES INTERNAS

!...Commons
!...Sentencias
    lugar(:,:)=0

    call DefinirLimites(InputProblem01%mop,InputProblem01%T)

    write (6,"(' ','Code',6x,'Parameter',52x,'Default',/,72x,'Value',/)")
    icod=99
    do while(icod/=0)
        call PresentarLimites(1)
 11     write (6,840)
		read (5,*,err=11) icod
		if(icod/=0) then
		    call ModificarLimites(icod)
		    write(6,"(' ','Code',6x,'Parameter',47x,'New Value',/)")
	    endif    
    enddo
    call GuardarLimites()
!---FORMATS
 840  format (' ',/,' If you want to change a parameter value give code',&
              ' number.','  If not 0.',/,71x,'> ',$)       
    
    endsubroutine

    subroutine AcceptStruct_Mixture (FinalStruct,ListPerformance,Aceptable)
    use PropertiesData
    use StructuresDesign
    implicit none
!   VARIABLES DE ENTRADA
    logical,intent(in)::FinalStruct
    type(CalculatedPerformance),pointer::ListPerformance
!   VARIABLES DE ENTRADA Y SALIDA
    logical,intent(inout)::Aceptable
    Aceptable=.False.
    if(FinalStruct)then
        if(LSel(1).and.ListPerformance%Selectivity<SelB(1))return
        if(LSel(2).and.ListPerformance%Selectivity>SelB(2))return
        if(LSolPow(1).and.ListPerformance%SolventPower<SolPowB(1))return
        if(LSolPow(2).and.ListPerformance%SolventPower>SolPowB(2))return
        if(LSLost(1).and.ListPerformance%SolventLost<SLostB(1))return
        if(LSLost(2).and.ListPerformance%SolventLost>SLostB(2))return  
        if(LDistCoef(1).and.ListPerformance%DistCoefficient<DistCoefB(1))return
        if(LDistCoef(2).and.ListPerformance%DistCoefficient>DistCoefB(2))return
        if(LVolat(1).and.ListPerformance%RelVolatility<VolatB(1))return
        if(LVolat(2).and.ListPerformance%RelVolatility>VolatB(2))return
        if(LDifBP(1).and.ListPerformance%DifferenceBP<DifBPB(1))return
        if(LDifBP(2).and.ListPerformance%DifferenceBP>DifBPB(2))return              
    else
        if(LSelInt(1).and.ListPerformance%Selectivity<SelIntB(1))return
        if(LSelInt(2).and.ListPerformance%Selectivity>SelIntB(2))return
        if(LSLostInt(1).and.ListPerformance%SolventLost<SLostIntB(1))return
        if(LSLostInt(2).and.ListPerformance%SolventLost>SLostIntB(2))return
        if(LVolatInt(1).and.ListPerformance%RelVolatility<VolatIntB(1))return
        if(LVolatInt(2).and.ListPerformance%RelVolatility<VolatIntB(2))return
    endif
    Aceptable=.True.
    return
    endsubroutine
    
!=================================================================================
    subroutine AcceptStruct_Pure (ListProperties,Aceptable)
    use PropertiesData
    use StructuresDesign
    implicit none
!   VARIABLES DE ENTRADA
    type(CalculatedProperties),pointer,intent(in)::ListProperties
!   VARIABLES DE ENTRADA Y SALIDA
    logical,intent(inout)::Aceptable
    Aceptable=.False.
    if(LBP(1).and.ListProperties%BoilingPoint<BPB(1))return
    if(LBP(2).and.ListProperties%BoilingPoint>BPB(2))return
    if(ListProperties%BoilingPoint==5000)return !por su falla la GC
    if(LPC(1).and.ListProperties%CriticalPressure<PCB(1))return
    if(LPC(2).and.ListProperties%CriticalPressure>PCB(2))return
    if(LVC(1).and.ListProperties%CriticalVolume<VCB(1))return
    if(LVC(2).and.ListProperties%CriticalVolume>VCB(2))return
    if(LTC(1).and.ListProperties%CriticalTemperature<TCB(1))return
    if(LTC(2).and.ListProperties%CriticalTemperature>TCB(2))return
    if(LVisc(1).and.ListProperties%Viscosity<ViscB(1))return
    if(LVisc(2).and.ListProperties%Viscosity>ViscB(2))return
    if(LDens(1).and.ListProperties%Density<DensB(1))return
    if(LDens(2).and.ListProperties%Density>DensB(2))return
    if(LHvap(1).and.ListProperties%Hvap<HvapB(1))return
    if(LHvap(2).and.ListProperties%Hvap>HvapB(2))return
    if(LKow(1).and.ListProperties%Kow<KowB(1))return
    if(LKow(2).and.ListProperties%Kow>KowB(2))return
    Aceptable=.True.
    return
    endsubroutine
    
    subroutine PresentarLimites(imp)
    use PropertiesData
    integer,intent(in)::imp !0=no imprimir, 1=imprimir

    
    j=1
    do i=1,size(limits)
        if(limits(i)%ExistRunValue) then
            if(imp==1)write(6,"(i2,1x,a8,': ',a50,a7,f8.2)")j, limits(i)%ID, limits(i)%cartel, limits(i)%unit, limits(i)%RunValue
            lugar(j,1)=i
            lugar(j,2)=0
            j=j+1
        endif
        if(limits(i)%ExistLowerBound) then
            if(imp==1)write(6,"(i2,1x,a8,': ','Minimun ',a42,a7,f8.2)")j, limits(i)%ID, limits(i)%cartel, limits(i)%unit, limits(i)%LowerBound
            lugar(j,1)=i
            lugar(j,2)=1
            j=j+1
        endif
        if(limits(i)%ExistUpperBound) then
            if(imp==1)write(6,"(i2,1x,a8,': ','Maximun ',a42,a7,f8.2)")j, limits(i)%ID, limits(i)%cartel, limits(i)%unit, limits(i)%UpperBound
            lugar(j,1)=i
            lugar(j,2)=2
            j=j+1
        endif
    enddo
    
    endsubroutine PresentarLimites
!
!
    subroutine ModificarLimites(icod)
    use PropertiesData
    integer,intent(in)::icod
    
    if(lugar(icod,2)==0)then
        write(6,*)limits(lugar(icod,1))%cartel
        read(5,*)limits(lugar(icod,1))%RunValue
    elseif(lugar(icod,2)==1)then
        write(6,*)'Minimun',limits(lugar(icod,1))%cartel
        read(5,*)limits(lugar(icod,1))%LowerBound    
    elseif(lugar(icod,2)==2)then
        write(6,*)'Maximun',limits(lugar(icod,1))%cartel
        read(5,*)limits(lugar(icod,1))%UpperBound
    endif
    
    
!"(' ',a63,7x,'> ',$)"
    
    endsubroutine ModificarLimites
  
!==========================================================================
    
    subroutine DefinirLimites(mop,Top)
    use PropertiesData
    integer,intent(in)::mop
    real*8, intent(in)::Top
    if(mop==0)then
    elseif(mop==1)then
        limits(1)%ExistUpperBound=.True.
        limits(1)%UpperBound=440. 
        limits(2)%ExistLowerBound=.True.
        limits(2)%LowerBound=300.         
        limits(8)%ExistLowerBound=.True.
        limits(8)%LowerBound=0.0
        limits(8)%ExistUpperBound=.True.
        limits(8)%UpperBound=150.        
        limits(9)%ExistLowerBound=.True.
        limits(9)%LowerBound=0.25
        limits(10)%ExistLowerBound=.True.
        limits(10)%LowerBound=1.
        limits(11)%ExistLowerBound=.True.
        limits(11)%LowerBound=0.00
       ! limits(12)%ExistUpperBound=.True.
       ! limits(12)%UpperBound=150.
        limits(13)%ExistUpperBound=.True.
        limits(13)%UpperBound=500.
        limits(14)%ExistLowerBound=.True.
        limits(14)%LowerBound=0.0
        limits(18)%ExistRunValue=.True.
        limits(18)%RunValue=0.9962
        limits(19)%ExistRunValue=.True.
        limits(19)%RunValue=100.
        limits(20)%ExistRunValue=.True.
        limits(20)%RunValue=21
        limits(21)%ExistRunValue=.True.
        limits(21)%RunValue=0
        limits(22)%ExistRunValue=.True.
        limits(22)%RunValue=Top
    elseif(mop==2)then
        limits(1)%ExistUpperBound=.True.
        limits(1)%UpperBound=240. 
        limits(15)%ExistLowerBound=.True.
        limits(15)%LowerBound=1.
        limits(16)%ExistLowerBound=.True.
        limits(16)%LowerBound=1.2
        limits(11)%ExistLowerBound=.True.
        limits(11)%LowerBound=15.
        limits(17)%ExistLowerBound=.True.
        limits(17)%LowerBound=30.
        limits(19)%ExistRunValue=.True.
        limits(19)%RunValue=100.
        limits(20)%ExistRunValue=.True.
        limits(20)%RunValue=21
        limits(21)%ExistRunValue=.True.
        limits(21)%RunValue=0
    elseif(mop==3)then
        limits(2)%ExistUpperBound=.True.
        limits(2)%UpperBound=Top
        limits(8)%ExistLowerBound=.True.
        limits(8)%LowerBound=3.
        limits(10)%ExistLowerBound=.True.
        limits(10)%LowerBound=0.6
    elseif(mop==4)then !PROVISORIO para reacciones!!!
        limits(2)%ExistLowerBound=.True.
        limits(2)%LowerBound=Top        
    endif
    endsubroutine

    subroutine GuardarLimites()
    use PropertiesData
    
    MWB(:)      =(/limits(1)%LowerBound ,limits(1)%UpperBound/)
    BPB(:)      =(/limits(2)%LowerBound ,limits(2)%UpperBound/)
    PCB(:)      =(/limits(3)%LowerBound ,limits(3)%UpperBound/)
    VCB(:)      =(/limits(4)%LowerBound ,limits(4)%UpperBound/)
    TCB(:)      =(/limits(5)%LowerBound ,limits(5)%UpperBound/)
    ViscB(:)    =(/limits(6)%LowerBound ,limits(6)%UpperBound/)
    HvapB(:)    =(/limits(7)%LowerBound ,limits(7)%UpperBound/)
    KowB(:)     =(/limits(8)%LowerBound ,limits(8)%UpperBound/)
    SelIntB(:)  =(/limits(9)%LowerBound ,limits(9)%UpperBound/) !intermediate
    SelB(:)     =(/limits(10)%LowerBound,limits(10)%UpperBound/)
    SolPowB(:)  =(/limits(11)%LowerBound,limits(11)%UpperBound/)
    SLostIntB(:)=(/limits(12)%LowerBound,limits(12)%UpperBound/) !intermediate
    SLostB(:)   =(/limits(13)%LowerBound,limits(13)%UpperBound/)
    DistCoefB(:)=(/limits(14)%LowerBound,limits(14)%UpperBound/)
    VolatIntB(:)=(/limits(15)%LowerBound,limits(15)%UpperBound/) !intermediate
    VolatB(:)   =(/limits(16)%LowerBound,limits(16)%UpperBound/)
    DifBPB(:)   =(/limits(17)%LowerBound,limits(17)%UpperBound/)

    LMW(:)      =(/limits( 1)%ExistLowerBound,limits( 1)%ExistUpperBound/)
    LBP(:)      =(/limits( 2)%ExistLowerBound,limits( 2)%ExistUpperBound/)
    LPC(:)      =(/limits( 3)%ExistLowerBound,limits( 3)%ExistUpperBound/)        
    LVC(:)      =(/limits( 4)%ExistLowerBound,limits( 4)%ExistUpperBound/)        
    LTC(:)      =(/limits( 5)%ExistLowerBound,limits( 5)%ExistUpperBound/)                    
    LVisc(:)    =(/limits( 6)%ExistLowerBound,limits( 6)%ExistUpperBound/)                        
    LHvap(:)    =(/limits( 7)%ExistLowerBound,limits( 7)%ExistUpperBound/)                
    LKow(:)     =(/limits( 8)%ExistLowerBound,limits( 8)%ExistUpperBound/)                    
    LSelInt(:)  =(/limits( 9)%ExistLowerBound,limits( 9)%ExistUpperBound/)  !intermediate
    LSel(:)     =(/limits(10)%ExistLowerBound,limits(10)%ExistUpperBound/)                
    LSolPow(:)  =(/limits(11)%ExistLowerBound,limits(11)%ExistUpperBound/)                            
    LSLostInt(:)=(/limits(12)%ExistLowerBound,limits(12)%ExistUpperBound/)  !intermediate     
    LSLost(:)   =(/limits(13)%ExistLowerBound,limits(13)%ExistUpperBound/)                
    LDistCoef(:)=(/limits(14)%ExistLowerBound,limits(14)%ExistUpperBound/)                        
    LVolatInt(:)=(/limits(15)%ExistLowerBound,limits(15)%ExistUpperBound/)  !intermediate                          
    LVolat(:)   =(/limits(16)%ExistLowerBound,limits(16)%ExistUpperBound/)                    
    LDifBP(:)   =(/limits(17)%ExistLowerBound,limits(17)%ExistUpperBound/)            

    endsubroutine

    subroutine initialization()
    use Boundaries
    use PropertiesData
!===Inicializazión
    limits(1) =bound('MW      ','Molecular Weight of the final structure              ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(2) =bound('BP      ','Boiling Boint                                        ','(K)    ',.False.,.False.,.False.,0.,0.,0.)
    limits(3) =bound('PC      ','Critical Pressure                                    ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(4) =bound('VC      ','Critical Volume                                      ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(5) =bound('TC      ','Critical Temperture                                  ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(6) =bound('Visc    ','Viscosity                                            ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(7) =bound('Hvap    ','Enthalpy of Vaporization                             ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(8) =bound('Kow     ','Octanol-Water Partition Coeficient                   ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(9) =bound('SelInt  ','Selectivity of an Intermediate Structure             ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(10)=bound('Sel     ','Selectivity of a Final Structure                     ','(wt)   ',.False.,.False.,.False.,0.,0.,0.)
    limits(11)=bound('SolPow  ','Solvent Power Required                               ','(wt)   ',.False.,.False.,.False.,0.,0.,0.)
    limits(12)=bound('SLostInt','Solvent Loss of an Intermediate Structure            ','(wt %) ',.False.,.False.,.False.,0.,0.,0.)
    limits(13)=bound('SLost   ','Solvent Loss of a Final Structure                    ','(wt %) ',.False.,.False.,.False.,0.,0.,0.)
    limits(14)=bound('DistCoef','Distribution Coefficient Required                    ','(wt %) ',.False.,.False.,.False.,0.,0.,0.)
    limits(15)=bound('VolatInt','Volatility with an Intermediate Structure            ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(16)=bound('Volat   ','Volatility with a Final Structucture                 ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(17)=bound('BPdiff  ','BP diff. between solvent and less volatile component ','(K)    ',.False.,.False.,.False.,0.,0.,0.)
    limits(18)=bound('DensRaff','Density of raffinate                                 ','(gr/ml)',.False.,.False.,.False.,0.,0.,0.)
    limits(19)=bound('Nsol    ','Maximum number of structures to be listed            ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(20)=bound('Is      ','Number of structures per page                        ','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(21)=bound('Out     ','Output type (0:Screen and output file, 1:Output file)','       ',.False.,.False.,.False.,0.,0.,0.)
    limits(22)=bound('Tem     ','Operation temperature                                ','(K)    ',.False.,.False.,.False.,0.,0.,0.)
   endsubroutine

endmodule evaluation