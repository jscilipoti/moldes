      MODULE GRPS
!      INTERFACE
        PARAMETER(NSEL1=150,NAA1=4,na1=150,nsva1=20)
        CHARACTER*70 titl
        INTEGER::MOP,NPEPIS,IDR(NSEL1),NYR(NSEL1),IDF(NSEL1),NYF(NSEL1),MSOL,MRAF,&
        MDV,MSV,MDV2,MDV1,MVAL(NSEL1),MSD(NSEL1),IFAM,NALAR,IALAR(NAA1),&
        NGM,NGM1V(NSEL1),JST1(NSEL1),KST1(NSEL1),IST1(NSEL1),HST(NSEL1),&
        MSV1,MJK,KILOUT,NGDV(0:NA1),NGSV(NSVA1),NGSV1(NSVA1),NSOL,IS
        real*8::TEM,BP1,XAZEO,A1,A2,A3,B1,B2,B3,DENS2
 !     ENDINTERFACE
      ENDMODULE
 !
    MODULE PropertiesData
        use Boundaries
        parameter(NMAX=60000)
        type (bound)::limits(22)
!===Número de grupos distintos en cada estructura generada
        integer::NGSL(NMAX)
!===Almacenamiento de estimaciones
!        1- MW:         molecular weight
!        2- BP:         boiling point
!        3- PC:         critical pressure
!        4- VC:         critical volume
!        5- TC:         critical temperture
!        6- Visc:       viscosity
!        7- Dens:       density
!        8- Hvap        enthalpy of vaporization
!        9- Kow:        octanol-water partition coeficient   
!	    10- SelInt:     selectivity of an intermediate structure
!		11- Sel:        selectivity of a final structure (wt)       
!		12- SolPow:     solvent power required (wt)                 
!		13- SLostInt:   solvent loss of an intermediate structure (wt %)      
!		14- SLost:      solvent loss of a final structure (wt %)              
!		15- DistCoef:   distribution coefficient required (wt %)
!       16- VolatInt:   volatility with an intermediate str.
!       17- Volat:      volatility with a final structucture
!       18- DifBPT:     boiling temperature difference between solvent and less volatile component (K)
        real*8::MW(NMAX),BP(NMAX),PC(NMAX),VC(NMAX),TC(NMAX),Visc(NMAX),Dens(NMAX),Hvap(NMAX),Kow(NMAX) !Puro
        real*8::SelInt(NMAX),Sel(NMAX),SolPow(NMAX),SLostInt(NMAX),SLost(NMAX),DistCoef(NMAX),&         !Mezcla
                          VolatInt(NMAX),Volat(NMAX),DifBP(NMAX),&
        DC2(NMAX) !Este último sirve para hacer el ranking de estructuras
!===LIMITES
!        1- MWB:        molecular weight
!        2- BPB:        boiling point
!        3- PCB:        critical pressure
!        4- VCB:        critical volume
!        5- TCB:        critical temperture
!        6- ViscB:      viscosity
!        7- DensB:      density
!        8- HvapB       enthalpy of vaporization
!        9- KowB:       octanol-water partition coeficient   
!	    10- SelIntB:    selectivity of an intermediate structure
!		11- SelB:       selectivity of a final structure (wt)       
!		12- SolPowB:    solvent power required (wt)                 
!		13- SLostIntB:  solvent loss of an intermediate structure (wt %)      
!		14- SLostB:     solvent loss of a final structure (wt %)              
!		15- DistCoefB:  distribution coefficient required (wt %)
!       16- VolatIntB:  volatility with an intermediate str.
!       17- VolatB:     volatility with a final structucture
!       18- DTC2SB:     Minimun boiling temperature difference between solvent and less volatile component (K)
        real*8::MWB(2),BPB(2),PCB(2),VCB(2),TCB(2),ViscB(2),DensB(2),HvapB(2),KowB(2),&
        SelIntB(2),SelB(2),SolPowB(2),SLostIntB(2),SLostB(2),DistCoefB(2),VolatIntB(2),VolatB(2),DifBPB(2)
!===Límites que serán tenidos en cuenta
!        1- LMW:        molecular weight
!        2- LBP:        boiling point
!        3- LPC:        critical pressure
!        4- LVC:        critical volume
!        5- LTC:        critical temperture
!        6- LVisc:      viscosity
!        7- LDens:      density
!        8- LHvap       enthalpy of vaporization
!        9- LKow:       octanol-water partition coeficient   
!	    10- LSelInt:    selectivity of an intermediate structure
!		11- LSel:       selectivity of a final structure (wt)       
!		12- LSolPow:    solvent power required (wt)                 
!		13- LSLostInt:  solvent loss of an intermediate structure (wt %)      
!		14- LSLost:     solvent loss of a final structure (wt %)              
!		15- LDistCoef:  distribution coefficient required (wt %)
!       16- LVolatInt:  volatility with an intermediate str.
!       17- LVolat:     volatility with a final structucture
!       18- LDTC2S:     Minimun boiling temperature difference between solvent and less volatile component (K)
        logical::LMW(2),LBP(2),LPC(2),LVC(2),LTC(2),LVisc(2),LDens(2),LHvap(2),LKow(2),&
        LSelInt(2),LSel(2),LSolPow(2),LSLostInt(2),LSLost(2),LDistCoef(2),LVolatInt(2),LVolat(2),LDifBP(2),&
        LDATO(18,2)

   ENDMODULE PropertiesData
    
    MODULE carteles
        character*63 cartel(18,2)

    ENDMODULE carteles
     
!
!       MODULE PropertiesData
!        parameter(NMAX=60000)
!!-------propiedades de compuesto puro
!!           MW: peso molecular; BP: punto de ebullición normal; PCF: presión crítica;
!!           VCF: volumen crítico; TCF: temperatura crítica; Visc: viscosidad; DENSIDF: densidad;
!!           Hvap: calor de vaporización; DENS: densidad; DLPOWF: Log(Pow)
!        real*8::MW(NMAX),BP(NMAX),PCF(NMAX),VCF(NMAX),TCF(NMAX),Visc(NMAX),&
!            DENSIDF(NMAX),Hvap(NMAX),DENS(NMAX),DLPOWF(NMAX)
!        integer::NGSL(NMAX)
!!-------propiedades de mezcla
!!           Sel: selectividad; SolPow: poder solvente; SLost: pérdida de solvente;  
!!           DistCoef y DC2: coeficiente de distribusión
!        real*8::Sel(NMAX),SolPow(NMAX),SLost(NMAX),DistCoef(NMAX),DC2(NMAX)
!!=======LIMITES
!!-------Límites propiedades compuesto puro (14)
!!           PMMALB y PMMA: peso molecular; BPLB y BPUB: pto de ebullición; PcLB y PcUB: p. crít.;
!!           VcLB y VcUB: v. crít.; TcLB y TcUB: T. crít.; VisUB y VisLB: viscosidad; 
!!           DensUB y DensLB: densidad; HVLB y HVUB: calor de vap; DLPOWLB y DLPOWUB: Log(Pow)
!!           DTC2S: mínima diferencia entre las TB del solvente y el componente menos volátil
!        real*8::PMMALB,PMMA=440.,BPLB,BPUB,PcLB,PcUB,VcLB,VcUB,TcLB,TcUBVisUB,VisLB,&
!            DensUB,DensLB,HVLB,HVUB,DLPOWLB,DLPOWUB,DTC2S
!!-------Límites propiedades mezcla (14)
!!	    	SCLLI : Minimum selectivity of an intermediate structure
!!	    	SELECT: Minimum selectivity of a final structure (wt)       
!!	    	SOLV  : Minimum solvent power required (wt)                 
!!	    	SLSUP1: Solvent loss of an intermediate structure (wt %)      
!!	    	SLSUPL: Solvent loss of a final structure (wt %)              
!!	    	DIST  : Minimum distribution coefficient required (wt %)      
!!	    	PMMA  : Maximum molecular weight of the final structure (wt)
!!	    	DENS2 : Density of raffinate (gr/ml)
!        real*8::SCLLI=.25,SELECT=1.,SOLV=.25,SLSUP1=150.,SLSUPL=150.,DIST=0.02,DENS2=1.0
!!-------Límites que serán tenidos en cuenta
!!       Propiedades compuesto puro
!        logical::LPMT=.False.,LBP=.False.,LPCF=.False.,LVCF=.False.,LTCF=.False.,LVISCF=.False.,&
!            LDENSIDF=.False.,LHVB=.False.,LDENS=.False.,LDLOWF=.False.
!!       Propiedades de mezcla
!        logical::LSCLLI=.False.,LSELECT=.False.,LSOLV=.False.,LSLSUP1=.False.,LSLSUPL=.False.,&
!            LDIST=.False.,LDENS2=.False.,LDTC2S=.False.
!      ENDMODULE