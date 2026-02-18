! ##############################################################################################
! ###### Hybrid Constitutive Material Model For Simultaneous Application of Hardwood     #######
! ###### (European Beech) and Softwood (Spruce) Incorporating All Deformation Mechanisms #######
! ###### Namely as: Elastic - Hygro-Expnasion - Visco-elastic - Mechano-sorptive -       #######
! ###### Multisurface Plasticity to Simulate Hybrid Adaptive Wood Elements               #######
! ##############################################################################################
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, 
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      CHARACTER*80 CMNAME
      INTEGER NTENS,NSTATV,NPROPS
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1  DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2  STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3  PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3) 
      DOUBLE PRECISION DTIME
	  
! Pick Wood Model
      Call UMAT_WOOD(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, 
     1 DDSDDT,DRPLDE,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,
     2 DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     3 DROT,DFGRD0,DFGRD1,NOEL,NPT,KSTEP,KINC)	 
	 
      RETURN
      END

! *********************************************************************************************
! ****** UMAT for Abaqus/Standard Incorporating Orthotropic Elasticity, Hygro-expansion, ****** 
! ****** Visco-elasticity, Mechano-sorptive creep, and Multi-surface Plasticity with     ******
! ****** Arbitrary Isotropic Hardening. Implicit Integration (Closest Point Projection   ******
! ****** Integration Scheme With Consistent Jacobian(Tangent Operator)is utilized        ******  
      SUBROUTINE UMAT_WOOD(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, 
     1  DDSDDT,DRPLDE,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,
     2  DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     3  DROT,DFGRD0,DFGRD1,NOEL,NPT,KSTEP,KINC)
      IMPLICIT NONE
	  
!      INCLUDE 'ABA_PARAM.INC'
	  
      CHARACTER*80 CMNAME
      INTEGER NTENS,NSTATV,NPROPS,KSTEP,KINC,NOEL,NPT,NDI,NSHR
      Double Precision STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3) 
      Double precision DTIME,TEMP,DTEMP,SSE,SPD,SCD
      Integer KW_CoSys
      Integer K1,K2,K3,Z1
      Integer NITER
      Integer, Dimension(3) :: KW_JAct_K,KW_Index	  
      Integer K4,J1,J2
      Integer NITER_EHP,KW_Countor,KW_num,KW_FCounter
      Integer KW_Counter2
      Integer KW_Res,KW_IA,KW_A,KW_depa,KW_Dealloc	
      Integer KW_PlasIdent, KW_ViscoIdent, KW_MechanoIdent
! #############################################################
      Double precision KW_Moist,KW_DMoist
      Double precision KW_J0L_1,KW_J0L_0
      Double precision KW_Refw,KW_ETT12
      Double precision KW_JOveral,KW_UFS,KW_UFSHE
      Double precision TOL1,TOL2,TOL3
      Double precision KW_T1,KW_T0
      Double precision KW_M1,KW_M0
      Double precision KW_norm1,KW_norm2,KW_norm3	  
      Double precision KW_n3
      Double Precision KW_MC0,KW_MC1
      Double Precision KW_MCJL0,KW_MCJL1	  
      Double precision KW_E10,KW_E20,KW_E30,KW_G120,KW_G130,KW_G230  
      Double precision KW_nu120,KW_nu130,KW_nu230
      Double precision KW_E11,KW_E21,KW_E31,KW_G121,KW_G131,KW_G231
      Double precision KW_nu121,KW_nu131,KW_nu231,KW_nu321
      Double precision KW_nu210,KW_nu310,KW_nu320,KW_nu211,KW_nu311
      Double precision KW_Coeff1,KW_Coeff2,KW_Coeff3
      Double precision Zero,One,Two,Three	  
! ------------------------------------------------------------	  
      Double precision KW_Hard1Out,KW_Hard2Out,KW_Hard3Out
      Double Precision KW_Q1Out,KW_Q2Out,KW_Q3Out
      Double precision KW_R_Norm
      Double Precision KW_fcr,KW_fct,KW_ftl,KW_fcl
      Double Precision KW_fvrt,KW_fvtl,KW_fvrl	  
! #############################################################	  
      Double precision, Dimension(6) :: KW_HygStr,KW_DHygStr
      Double precision, Dimension(4) :: KW_Tau,KW_Gamma_1,KW_Gamma_0
      Double precision, Dimension(6) :: KW_ElasStr_Old
      Double precision, Dimension(6) :: KW_Stress_Old,KW_Stran_Old
      Double precision, Dimension(4) :: KW_Jve
      Double precision, Dimension(6) :: KW_TotVEStr_New,KW_TotVEStr_Old
      Double precision, Dimension(6) :: KW_TotMSStr_New,KW_TotMSStr_Old
      Double precision, Dimension(3) :: KW_Mu
      Double precision, Dimension(4) :: KW_AllT1,KW_AllT0
      Double precision, Dimension(3) :: KW_AllM1,KW_AllM0	  
      Double precision, Dimension(48) :: KW_ResVec
! ------------------------------------------------------------	  
      Double precision, Dimension(6) :: KW_a1,KW_a2,KW_a3	  
      Double precision, Dimension(6) :: KW_Stran_Pl,KW_Stran_El
      Double precision, Dimension(6) :: KW_Stran_El_New
      Double precision, Dimension(6) :: KW_Stran_El_New2
      Double precision, Dimension(6) :: KW_Stran_Pl_New
      Double precision, Dimension(6) :: KW_DSTRESS_Elas,KW_STRESS2 
      Double precision, Dimension(6) :: KW_STRAN_New
      Double precision, Dimension(6) :: KW_EPS_Trial
      Double precision, Dimension(6) :: KW_dEPL
      Double precision, Dimension(18) :: KW_All_a
      Double precision, Dimension(6) :: KW_TotVisc_Old, KW_TotMech_Old
      Double precision, Dimension(6) :: KW_TotVisc_New, KW_TotMech_New	  
! #############################################################
      Double precision, Dimension(1,1) :: KW_DTE
      Double precision, Dimension(1,1) :: KW_DEE
      Double precision, Dimension(1,1) :: KW_DPE	  
      Double precision, Dimension(6,1) :: KW_R1
      Double precision, Dimension(6,1) :: KW_ElasStr_New
      Double precision, Dimension(6,1) :: KW_Stress_New
      Double precision, Dimension(6,1) :: KW_R2i
      Double precision, Dimension(6,1) :: KW_R2
      Double precision, Dimension(6,1) :: KW_R3j
      Double precision, Dimension(6,1) :: KW_R3
      Double precision, Dimension(6,1) :: KW_RSum
      Double precision, Dimension(6,1) :: KW_DStress
      Double precision, Dimension(6,1) :: KW_R31,KW_R32,KW_R33
      Double precision, Dimension(6,1) :: KW_DMSj1,KW_DMSj2,KW_DMSj3
      Double precision, Dimension(6,1) :: KW_DSTRHyg
      Double precision, Dimension(6,1) :: KW_DElasSTR_NO
      Double precision, Dimension(6,1) :: KW_DPlasSTR_NO	  
      Double precision, Dimension(1,4) :: KW_T1Gamma
      Double precision, Dimension(1,6) :: KW_DVEi
      Double precision, Dimension(1,6) :: KW_Stress_Avg
      Double precision, Dimension(6,4) :: KW_MatR2	  
      Double precision, Dimension(6,6) :: KW_Compl_0,KW_Compl_1
      Double precision, Dimension(6,6) :: KW_MOE_0,KW_MOE_1
      Double precision, Dimension(4,6) :: KW_VisStr_Old,KW_VisStr_New
      Double precision, Dimension(6,6) :: KW_MechCompl11,KW_MechCompl01
      Double precision, Dimension(6,6) :: KW_MechCompl12,KW_MechCompl02
      Double precision, Dimension(6,6) :: KW_MechCompl13,KW_MechCompl03
      Double precision, Dimension(6,6) :: KW_InvCjMech11
      Double precision, Dimension(6,6) :: KW_InvCjMech12
      Double precision, Dimension(6,6) :: KW_InvCjMech13
      Double precision, Dimension(6,6) :: KW_ElasTanOpt	  
      Double precision, Dimension(6,6) :: KW_ViscoTanOpt
      Double precision, Dimension(6,6) :: KW_MechTanOpt
      Double precision, Dimension(6,6) :: KW_InvCT,KW_CT
      Double precision, Dimension(4,6) :: KW_DVEStr,KW_DVEStr1	  
      Double precision, Dimension(3,6) :: KW_MechStr_Old,KW_MechStr_New
      Double precision, Dimension(3,6) :: KW_DMSStr,KW_DMSStr1
      Double precision, Dimension(6,3) :: KW_MatR3
      Double precision, Dimension(3,2) :: KW_mTmL
      Double precision, Dimension(6,6) :: KW_EPTanOpt
      Double precision, Dimension(6,6) :: KW_InvEPTanOpt		  
! ------------------------------------------------------------	  
      Double precision, Dimension(6,6) :: KW_b1,KW_b2,KW_b3	  
      Double precision, Dimension(6,6) :: KW_D2FDS2Out
      Double precision, Dimension(6,6) :: KW_DDSDDEOutMat
      Double precision, Dimension(6,6) :: KW_Hessain
      Double precision, Dimension(6,1) :: KW_H61, KW_Ee_Trial
      Double precision, Dimension(6,1) :: KW_ElasticE	  
      Double precision, Dimension(6,1) :: KW_Str_Trial
      Double precision, Dimension(6,1) :: KW_Stress_NewEHP
      Double precision, Dimension(6,1) :: KW_DFDS
      Double precision, Dimension(3,1) :: KW_Eq_Eplas,KW_Eq_Eplas_New
      Double precision, Dimension(3,1) :: KW_dAlpha
      Double precision, Dimension(3,1) :: KW_Plas_Mult,KW_Plas_Mult_Bar 
      Double precision, Dimension(3,1) :: KW_YFuncs_Trial,KW_YFuncs
      Double precision, Dimension(3,1) :: KW_Cur_YFuncs
      Double precision, Dimension(3,1) :: KW_AllHarden
      Double precision, Dimension(3,1) :: KW_DFDQ
      Double precision, Dimension(3,1) :: KW_All_Q,KW_All_fc1
      Double precision, Dimension(3,1) :: KW_All_fcq
      Double precision, Dimension(3,1) :: KW_dPlas_Mult
      Double precision, Dimension(3,1) :: KW_dQdw0	  
      Double precision, Dimension(1,1) :: KW_YFValTr,KW_YFValK
      Double precision, Dimension(6,3) :: KW_All_DFDS
      Double precision, Dimension(6,18) :: KW_All_b
      Double precision, Dimension(6,18) :: KW_All_D2FDS2
! #############################################################	  
      Double precision, Allocatable :: KW_Residual(:,:)
      Double precision, Allocatable :: KW_InvA(:,:)
      Double precision, Allocatable :: KW_AMatrix(:,:)
      Double precision, Allocatable :: KW_dEPL_Alpha(:,:)
      Double precision KW_PlasCoeff, KW_ViscoCoeff, KW_MechanoCoeff
! #############################################################	  
      Parameter (Zero = 0.D0, One = 1.D0, Two = 2.D0, Three = 3.D0)
      !Parameter (KW_UFS = 1.D0)    ! ###### For Mechano-sorption Calculation default (20.D0) 	
      Parameter (KW_UFSHE = 0.3D0) ! Fiber Saturation Point ###### For Hygro-expansion Calculation ######

! ###########################################################################################
! ######################### Switch On/Off different deformation modes #######################
! ###########################################################################################
      !KW_ViscoIdent = 0 ! If =1 viscoelasticity is switched on and =0 viscoelasticity is off
      INCLUDE 'activate_ve_creep.inc' ! If KW_ViscoIdent=1 viscoelasticity is switched on and KW_ViscoIdent=0 viscoelasticity is off
      KW_MechanoIdent = 0   ! If =1 mechanosorption is switched on and =0 mechanosorption is off       
      ! Set viscoelastic coefficient
      If (KW_ViscoIdent .EQ. 1) Then 
        KW_ViscoCoeff = (1.D0) 
      Else
        KW_ViscoCoeff = (10.D60) 
      End If
      ! Set mechanosorption coefficient
      If (KW_MechanoIdent .EQ. 1) Then 
        KW_MechanoCoeff = (1.D0) 
      Else 
        KW_MechanoCoeff = (10.D30) 
      End If

! ###########################################################################################
! ##############################   Hygro-expansion Calculation   ############################
! ###########################################################################################
      ! Make analogy between Temperature and Moisture at t(n)
      KW_Moist = TEMP
      KW_DMoist = DTEMP      
      ! Store initial Moisture Content as STATEV(1) at t(0)
      IF (KSTEP .EQ. 1) THEN
        IF (KINC .EQ. 1) THEN
    	    IF (KW_Moist .LE. (KW_UFSHE)) THEN
            STATEV(1) = KW_Moist
          ELSE IF (KW_Moist .GT. (KW_UFSHE)) THEN
    	      STATEV(1) = (KW_UFSHE)
          END IF
        END IF
      END IF
      ! Set hygro-expansion coefficients (local coord. system)
      ! 1 = Fiber direction
      ! 2 = Tangential direction
      ! 3 = Longitudinal direction
      IF (KW_Moist .LE. (KW_UFSHE)) THEN  ! if below FSP
        IF (CMNAME == 'ML') THEN
	        KW_Coeff1 = (0.19D0)
          KW_Coeff2 = (0.19D0)
          KW_Coeff3 = (0.32D0)
        ELSE IF (CMNAME == 'P') THEN
          KW_Coeff1 = (0.19D0)
          KW_Coeff2 = (0.19D0)
          KW_Coeff3 = (0.32D0)          
        ELSE IF (CMNAME == 'S1') THEN
	        KW_Coeff1 = (0.006D0)
          KW_Coeff2 = (0.479D0)
          KW_Coeff3 = (0.479D0)
        ELSE IF (CMNAME == 'S2') THEN
	        KW_Coeff1 = (0.003D0)
          KW_Coeff2 = (0.367D0)
          KW_Coeff3 = (0.367D0)
        ELSE IF (CMNAME == 'S3') THEN
	        KW_Coeff1 = (0.003D0)
          KW_Coeff2 = (0.455D0)
          KW_Coeff3 = (0.455D0)
        ELSE IF (CMNAME == 'TEST') THEN
	        KW_Coeff1 = (1.D0)
          KW_Coeff2 = (1.D0)
          KW_Coeff3 = (1.D0)

        ELSE
          PRINT *, 'Error: Material "', CMNAME, '" not supported.'
          RETURN
        END IF ! End material selection	
		 
      ELSE IF (KW_Moist .GT. (KW_UFSHE)) THEN
        KW_Coeff1 = (0.D0)          
        KW_Coeff2 = (0.D0)          
        KW_Coeff3 = (0.D0)
      END IF
      ! Calculate hygro-expansion strain tensor at t(n+1)
      KW_HygStr = (Zero)
      KW_HygStr(1) = KW_Coeff1*((KW_Moist+KW_DMoist)-STATEV(1))
      KW_HygStr(2) = KW_Coeff2*((KW_Moist+KW_DMoist)-STATEV(1))
      KW_HygStr(3) = KW_Coeff3*((KW_Moist+KW_DMoist)-STATEV(1))
      KW_HygStr(4) = (0.D0)
      KW_HygStr(5) = (0.D0)
      KW_HygStr(6) = (0.D0)
      ! Calculate variation of hygro-expansion strain tensor
      KW_DHygStr = (Zero)
      KW_DHygStr(1) = KW_Coeff1*(KW_DMoist)
      KW_DHygStr(2) = KW_Coeff2*(KW_DMoist)
      KW_DHygStr(3) = KW_Coeff3*(KW_DMoist)
      KW_DHygStr(4) = (0.D0)
      KW_DHygStr(5) = (0.D0)
      KW_DHygStr(6) = (0.D0)	  
! ############################## END hygro-expansion Calculation ############################


! ###########################################################################################
! #############################   Elastic Strain Calculation    #############################
! ###########################################################################################
      ! Set Moisture Content at t(n)
      KW_MC0 = (0.D0)
      KW_MC0 = (KW_Moist)
      ! Compute elastic compliance tensor as a function of Moisture Content at t(n)
      KW_Compl_0 = (0.D0)
      Call KW_ElastEngConst(KW_MC0, CMNAME, KW_Compl_0) ! KW_Compl_0 = Compliance tensor at t(n)
      ! Compute elastic stiffness tensor from compliance tensor at t(n) 
      KW_MOE_0 = (0.D0)
      Call KW_Inv6(KW_Compl_0, KW_MOE_0) ! KW_MOE_0 = Stiffness tensor at t(n)

      ! Set Moisture Content at t(n+1)
      KW_MC1 = (0.D0)
      KW_MC1 = (KW_Moist+KW_DMoist)
      ! Compute elastic compliance tensor as a function of Moisture Content at t(n+1)
      KW_Compl_1 = (0.D0)
      Call KW_ElastEngConst(KW_MC1, CMNAME, KW_Compl_1) ! KW_Compl_1 = Compliance tensor at t(n+1)
      ! Compute elastic stiffness tensor from compliance tensor at t(n+1) 
      KW_MOE_1 = (0.D0)
      Call KW_Inv6(KW_Compl_1, KW_MOE_1) ! KW_MOE_1 = Stiffness tensor at t(n+1)

! ###########################################################################################
! ###########################   Viscoelastic Strain Calculation    ##########################
! ###########################################################################################
      ! Store elastic compliance in fiber direction at t(n)
      KW_J0L_0 = KW_Compl_0(1,1)
      ! Set Moisture Content at t(n)
      KW_MCJL0 = (0.D0)
      KW_MCJL0 = (KW_Moist)
      ! Compute viscoelastic parameters as a function of Moisture Content at t(n)
      Call KW_ViscoelastParam(KW_MCJL0, CMNAME, KW_J0L_0, KW_ViscoCoeff,
     1 KW_Gamma_0, KW_Tau)

      ! Store elastic compliance in fiber direction at t(n+1)
      KW_J0L_1 = 1. !KW_Compl_1(1,1)
      ! Set Moisture Content at t(n+1)	 
      KW_MCJL1 = (0.D0)
      KW_MCJL1 = (KW_Moist+KW_DMoist)
      ! Compute viscoelastic parameters as a function of Moisture Content at t(n+1)
      Call KW_ViscoelastParam(KW_MCJL1, CMNAME, KW_J0L_1, KW_ViscoCoeff,
     1 KW_Gamma_1, KW_Tau)
      ! Compute viscoelastic compliance tensor as a function of Moisture Content at t(n+1)
      Do K1=1, 4	     
	      KW_Jve(K1) = (KW_J0L_1/KW_Gamma_1(K1))*(One-Exp(-(Time(1)+DTIME)/KW_Tau(K1)))   
      End Do	  
      KW_JOveral = Sum(KW_Jve(1:4)) ! KW_JOveral = Viscoelastic compliance tensor at t(n+1)

! #############################################################################################
! Set initial values of allowed number of iterations and acceptable tolerances
      NITER_EHP = 100
      NITER = 50
      TOL1 = (1.0D-6)
      TOL2 = (1.0D-6)
      TOL3 = (1.0D-6)
! #############################################################################################
! ############## To Store old values of following components as State Variables: ##############
! ############## 1) Initial moisture             ---> STATEV(1)                  ##############
! ############## 2) Elastic strain               ---> STATEV(2:7)                ##############
! ############## 3) Plastic strain               ---> STATEV(8:13)               ##############
! ############## 4) Hardening variables          ---> STATEV(14:16)              ##############
! ############## 5) Visco-elastic strain         ---> STATEV(17:40)              ##############
! ############## 6) Mechano-sorption strain      ---> STATEV(41:58)              ##############
! ############## 7) Total viscoelastic           ---> STATEV(59:64)              ##############
! ############## 8) Total mechanosorption        ---> STATEV(65:70)              ##############
! #############################################################################################
      KW_ElasStr_Old(1:6) = STATEV(2:7) ! Store Elastic Strain at t(n)
      !KW_Stran_Pl(1:6)      = STATEV(8:13)
      !KW_Eq_Eplas(1:3,1)    = STATEV(14:16)		 
      KW_VisStr_Old(1,1:6)  = STATEV(17:22) ! Viscoelastic strain of 1st KV elem.
      KW_VisStr_Old(2,1:6)  = STATEV(23:28) ! Viscoelastic strain of 2nd KV elem.
      KW_VisStr_Old(3,1:6)  = STATEV(29:34) ! Viscoelastic strain of 3rd KV elem.
      KW_VisStr_Old(4,1:6)  = STATEV(35:40)	! Viscoelastic strain of 4th KV elem.
      !KW_MechStr_Old(1,1:6) = STATEV(41:46) ! Mechanosorptive strain of 1st KV elem.
      !KW_MechStr_Old(2,1:6) = STATEV(47:52) ! Mechanosorptive strain of 2nd KV elem.
      !KW_MechStr_Old(3,1:6) = STATEV(53:58) ! Mechanosorptive strain of 3rd KV elem.
      KW_TotVEStr_Old(1:6)  = STATEV(59:64)	! Total viscoelastic strain
      !KW_TotMSStr_Old(1:6)  = STATEV(65:70) ! Total mechanosorptive strain

      ! Store values of stress and strain at the start of the time increment t(n)
      KW_Stress_Old = STRESS ! Store Total Stress at t(n)	
	    KW_Stran_Old = STRAN ! Store Total Strain at t(n)
      ! Store values of strain at the end of the time increment t(n+1)
      KW_STRAN_New = STRAN + DSTRAN ! Store total strain at t(n+1)
      KW_VisStr_New = KW_VisStr_Old ! Initialize viscoelastic strain
      !KW_MechStr_New = KW_MechStr_Old ! Initialize mechanosorptive strain

! ###########################################################################################
! ############################   Viscoelastic Tangent Operator    ###########################
! ###########################################################################################
      ! Calculate time functions at t(n) and t(n+1)
      KW_AllT1 = (Zero) ! Time functions at t(n+1)
      KW_AllT0 = (Zero) ! Time functions at t(n)
      Do K3=1, 4                                       			
        KW_T1 = (Zero)
    	  KW_T0 = (Zero)
        Call KW_TimeFun(KW_T1,KW_T0,DTIME,KW_Tau(K3)) 
	      KW_AllT1(K3) = KW_T1
		    KW_AllT0(K3) = KW_T0
      End Do 		 
      ! Calculate tangent operator	 
      KW_T1Gamma = (Zero)		 
      Do K3=1, 4                                       			
        KW_T1Gamma(1,K3) = (KW_AllT1(K3)/KW_Gamma_1(K3))
      End Do	  
      KW_ViscoTanOpt = (Zero)        		 
      KW_ViscoTanOpt = (Sum(KW_T1Gamma(1,1:4)))*KW_Compl_1 

! #############################################################################################
! ##############################  Computation of Stress Tensor   ##############################
! #############################################################################################
      Z1=0 
      Do Z1=1, NITER
        If (Z1 .GT. NITER) Then
     	    If ((NOEL .EQ. 1) .AND. (NPT .EQ. 1)) Then
		        Write(*,*) '--------------------------------------------'
            Write(*,*) '### Maximum number of iterations reached ###'
            Write(*,*) '--------------------------------------------'
          End If
          Exit
        End If
      
      ! Compute elastic tangent operator = Elastic compliance tensor (inverse of Jacobian matrix)
        KW_ElasTanOpt = (Zero)
        KW_ElasTanOpt = KW_Compl_1
      ! Compute new stress at t(n+1)
        KW_ElasStr_New = (Zero)
        Do K2=1, NTENS
          KW_ElasStr_New(K2,1) = KW_STRAN_New(K2)-KW_HygStr(K2)-(Sum(KW_VisStr_New(1:4,K2)))
        End Do	 
        KW_Stress_New = Matmul(KW_MOE_1,KW_ElasStr_New) ! Stress = elastic stiffness * elastic strain

      ! Compute Total Tangent Operator
        KW_InvCT = (Zero)
        KW_InvCT = KW_ElasTanOpt + KW_ViscoTanOpt
        KW_CT = (Zero)
        Call KW_Inv6(KW_InvCT, KW_CT)

      ! Calculate elastic residual vector = R1
        KW_R1 = (Zero)
        Call KW_ResR1(KW_R1,KW_Stress_New,KW_STRAN_New,KW_HygStr,KW_VisStr_New,KW_Compl_1)
      ! Calculate viscoelastic residual vector = R2
        ! Calculate R2 components
        KW_MatR2 = (Zero)
        Do K1=1, 4
	        KW_R2i = (Zero)
          Call KW_ResR2(KW_R2i, KW_VisStr_Old(K1,1:6), KW_VisStr_New(K1,1:6),
     1      KW_Stress_Old, KW_Stress_New, KW_Compl_0, KW_Compl_1, KW_Gamma_1(K1),
     2      KW_Gamma_0(K1), KW_AllT1(K1), KW_AllT0(K1), DTIME, KW_Tau(K1))
	        KW_MatR2(1:6,K1) = KW_R2i(1:6,1)
	        KW_MatR2(1:6,K1) = KW_R2i(1:6,1)
        End Do		 
        ! Calculate total R2
        KW_R2 = (Zero)
        Do K1=1, 6
          KW_R2(K1,1) = (Sum(KW_MatR2(K1,1:4)))
        End Do

      ! Calculate change of stress
        KW_RSum = (Zero)
        KW_RSum = (KW_R1+KW_R2) ! Total residuals
        KW_DStress = (Zero)
        KW_DStress = -Matmul(KW_CT, KW_RSum)	  
      ! Calculate change of viscoelastic strain
        KW_DVEStr1 = (Zero)
        KW_DVEStr = (Zero)
        Do K1=1, 4
    	    KW_DVEi = (Zero)
        Call KW_DeltaVE(KW_DVEi, KW_MatR2(1:6,K1), KW_Compl_1,
     1	    KW_DStress, KW_Gamma_1(K1), KW_AllT1(K1))
          KW_DVEStr1(K1,1:6) = KW_DVEi(1,1:6)
        End Do	  
        ! Consider the definition of engineering shear strain
        KW_DVEStr(1,1:3) = KW_DVEStr1(1,1:3)	  
        KW_DVEStr(1,4:6) = KW_DVEStr1(1,4:6)		 
        KW_DVEStr(2,1:3) = KW_DVEStr1(2,1:3)	  
        KW_DVEStr(2,4:6) = KW_DVEStr1(2,4:6)
        KW_DVEStr(3,1:3) = KW_DVEStr1(3,1:3)	  
        KW_DVEStr(3,4:6) = KW_DVEStr1(3,4:6)
        KW_DVEStr(4,1:3) = KW_DVEStr1(4,1:3)	  
        KW_DVEStr(4,4:6) = KW_DVEStr1(4,4:6)

      ! Update viscoelastic and mechanosorptive strain
        KW_VisStr_New = KW_VisStr_New + KW_DVEStr
        !KW_MechStr_New = KW_MechStr_New + KW_DMSStr

      ! Update stress based on the updated values of viscoelastic and mechanosorptive strain
        KW_ElasStr_New = (Zero)
        Do K2=1, 6
          KW_ElasStr_New(K2,1) = KW_STRAN_New(K2)-KW_HygStr(K2)-(Sum(KW_VisStr_New(1:4,K2)))
        End Do	 
        KW_Stress_New = Matmul(KW_MOE_1,KW_ElasStr_New)

      ! Re-compute residual vectors
        ! Calculate elastic residual vector = R1
        KW_R1 = (Zero)
        Call KW_ResR1(KW_R1, KW_Stress_New, KW_STRAN_New, KW_HygStr,
     1    KW_VisStr_New, KW_Compl_1)  
        ! Calculate viscoelastic residual vector = R2
        KW_MatR2 = (Zero)
        Do K1=1, 4
	        KW_R2i = (Zero)
          Call KW_ResR2(KW_R2i,KW_VisStr_Old(K1,1:6),
     1	    KW_VisStr_New(K1,1:6),KW_Stress_Old,KW_Stress_New,
     2	    KW_Compl_0,KW_Compl_1,KW_Gamma_1(K1),KW_Gamma_0(K1),
     3	    KW_AllT1(K1),KW_AllT0(K1),DTIME,KW_Tau(K1))
	        KW_MatR2(1:6,K1) = KW_R2i(1:6,1)
        End Do     
      ! Construct residual vector
        KW_ResVec = (Zero)
        KW_ResVec(1:6) = KW_R1(1:6,1)
        KW_ResVec(7:12)  = KW_MatR2(1:6,1)
        KW_ResVec(13:18) = KW_MatR2(1:6,2)
        KW_ResVec(19:24) = KW_MatR2(1:6,3)
        KW_ResVec(25:30) = KW_MatR2(1:6,4)
      ! Compute Norm of the residual vector
        KW_n3 = (Zero)
        Do K2=1, 30 !48
	        KW_n3 = KW_n3 + (KW_ResVec(K2))**2
        End Do
        KW_norm3 = DSqrt(KW_n3)

        If (KW_norm3 .LT. TOL3) Then

          ! Define the Jacobian and the Stress at t(n+1)
          STRESS(1:6) = KW_Stress_New(1:6,1)  
          DDSDDE = KW_CT

          ! Calculate total viscoelastic strain
          Do K1=1, 6      
            KW_TotVEStr_New(K1) = Sum(KW_VisStr_New(1:4,K1))
            !KW_TotMSStr_New(K1) = Sum(KW_MechStr_New(1:3,K1))
          End Do
        
          ! Update State Variables to the values at t(n+1)
          STATEV(2:7)   = KW_ElasStr_New(1:6,1)
          STATEV(17:22) = KW_VisStr_New(1,1:6)
          STATEV(23:28) = KW_VisStr_New(2,1:6)
          STATEV(29:34) = KW_VisStr_New(3,1:6) 
          STATEV(35:40) = KW_VisStr_New(4,1:6)
          !STATEV(41:46) = KW_MechStr_New(1,1:6) 
          !STATEV(47:52) = KW_MechStr_New(2,1:6) 
          !STATEV(53:58) = KW_MechStr_New(3,1:6)		 
          STATEV(59:64) = KW_TotVEStr_New(1:6)
          !STATEV(65:70) = KW_TotMSStr_New(1:6)

          ! Calculate specific energy
          ! Calculate change in total specific energy
          KW_DTE = (Zero)
          Do K3=1, 6
            ! Hygroexpansion strain has no contribution in the change of specific energy
            KW_DSTRHyg(K3,1) = (DSTRAN(K3)-KW_DHygStr(K3))
            KW_Stress_Avg(1,K3) = (0.5D0)*(KW_Stress_New(K3,1)+KW_Stress_Old(K3))
          End Do
          KW_DTE = Matmul(KW_Stress_Avg,KW_DSTRHyg)	   	 
          ! Change in elastic specific energy
          KW_DEE = (Zero)
          Do K1=1, 6 
            KW_DElasSTR_NO(K1,1)=(KW_ElasStr_New(K1,1)-KW_ElasStr_Old(K1)) 
          End Do		 
              KW_DEE = Matmul(KW_Stress_Avg,KW_DElasSTR_NO)
          ! Calculate elastic strain energy
          SSE = SSE + KW_DEE(1,1)
          ! Calculate creep dissipation energy
          SCD = SCD +(KW_DTE(1,1)-KW_DEE(1,1))

          Exit		 
        End If!(KW_norm3 .LT. TOL3)
      


      End Do ! Z1=1, NITER
	  
      RETURN
      END




! ===========================================================================================

! ###########################################################################################
! ##############################     Inverse of a 6x6 matrix     ############################
! ###########################################################################################
      SUBROUTINE KW_Inv6(KW_DI,KW_CO)
      Implicit none	
	  
      Double precision KW_DDET,KW_DDET1,KW_D1
      Double precision, Dimension(6,6), Intent(IN) :: KW_DI
      Double precision, Dimension(6,6), Intent(OUT) :: KW_CO
	  
      KW_DDET1 = KW_DI(1,1)*KW_DI(2,2)*KW_DI(3,3) 
     1  -KW_DI(1,1)*KW_DI(2,3)*KW_DI(2,3)
     2  -KW_DI(2,2)*KW_DI(1,3)*KW_DI(1,3) 
      KW_DDET = KW_DDET1-KW_DI(3,3)*KW_DI(1,2)*KW_DI(1,2)
     1  +(2.D0)*KW_DI(1,2)*KW_DI(2,3)*KW_DI(1,3) 
      KW_D1 = (1.D0)/KW_DDET
      KW_CO(1,1) = (KW_DI(2,2)*KW_DI(3,3)-KW_DI(2,3)*KW_DI(2,3))*KW_D1 
      KW_CO(1,2) = (KW_DI(1,3)*KW_DI(2,3)-KW_DI(1,2)*KW_DI(3,3))*KW_D1
      KW_CO(1,3) = (KW_DI(1,2)*KW_DI(2,3)-KW_DI(1,3)*KW_DI(2,2))*KW_D1
      KW_CO(2,1) = (KW_DI(1,3)*KW_DI(2,3)-KW_DI(1,2)*KW_DI(3,3))*KW_D1
      KW_CO(3,1) = (KW_DI(1,2)*KW_DI(2,3)-KW_DI(1,3)*KW_DI(2,2))*KW_D1
      KW_CO(2,2) = (KW_DI(3,3)*KW_DI(1,1)-KW_DI(1,3)*KW_DI(1,3))*KW_D1
      KW_CO(2,3) = (KW_DI(1,2)*KW_DI(1,3)-KW_DI(2,3)*KW_DI(1,1))*KW_D1
      KW_CO(3,2) = (KW_DI(1,2)*KW_DI(1,3)-KW_DI(2,3)*KW_DI(1,1))*KW_D1
      KW_CO(3,3) = (KW_DI(1,1)*KW_DI(2,2)-KW_DI(1,2)*KW_DI(1,2))*KW_D1
      KW_CO(4,4) = (1.D0)/KW_DI(4,4)
      KW_CO(5,5) = (1.D0)/KW_DI(5,5)
      KW_CO(6,6) = (1.D0)/KW_DI(6,6)

      RETURN
      End
! ###########################################################################################

! ###########################################################################################
! ##############################    Elastic compliance tensor    ############################
! ###########################################################################################
      SUBROUTINE KW_ElastEngConst(KW_MC, MATNAME, KW_Compl)
      Implicit none	
	  
      ! Define a parameter for unity
      DOUBLE PRECISION, PARAMETER :: One = 1.D0
      ! Input parameters
      DOUBLE PRECISION, INTENT(IN) :: KW_MC
      CHARACTER*(*)      , INTENT(IN) :: MATNAME
      ! Output: Elastic compliance tensor (6x6 matrix)
      Double precision, Dimension(6,6), Intent(OUT) :: KW_Compl
      ! Local variables for engineering constants
      DOUBLE PRECISION :: KW_E1, KW_E2, KW_E3
      DOUBLE PRECISION :: KW_G12, KW_G13, KW_G23
      DOUBLE PRECISION :: KW_nu12, KW_nu13, KW_nu23

	    IF (MATNAME == 'ML') THEN
	      KW_E1 = (2806146977.943D0)*(KW_MC)**10 + (-4641648853.314D0)*(KW_MC)**9 + (3303116949.83D0)*(KW_MC)**8 +
     1  (-1317562279.701D0)*(KW_MC)**7 + (320958933.22D0)*(KW_MC)**6 + (-48569827.335D0)*(KW_MC)**5 +
     2  (4399632.672D0)*(KW_MC)**4 + (-210611.458D0)*(KW_MC)**3 + (3478.562D0)*(KW_MC)**2 + (14.55D0)*(KW_MC) + (5.689D0) 
	      KW_E2 = (2806146977.943D0)*(KW_MC)**10 + (-4641648853.314D0)*(KW_MC)**9 + (3303116949.83D0)*(KW_MC)**8 +
     1  (-1317562279.701D0)*(KW_MC)**7 + (320958933.22D0)*(KW_MC)**6 + (-48569827.335D0)*(KW_MC)**5 +
     2  (4399632.672D0)*(KW_MC)**4 + (-210611.458D0)*(KW_MC)**3 + (3478.562D0)*(KW_MC)**2 + (14.55D0)*(KW_MC) + (5.689D0) 
	      KW_E3 = (2806146977.943D0)*(KW_MC)**10 + (-4641648853.314D0)*(KW_MC)**9 + (3303116949.83D0)*(KW_MC)**8 +
     1  (-1317562279.701D0)*(KW_MC)**7 + (320958933.22D0)*(KW_MC)**6 + (-48569827.335D0)*(KW_MC)**5 +
     2  (4399632.672D0)*(KW_MC)**4 + (-210611.458D0)*(KW_MC)**3 + (3478.562D0)*(KW_MC)**2 + (14.55D0)*(KW_MC) + (5.689D0) 
	      KW_nu23 = (14033992.731D0)*(KW_MC)**10 + (-26688932.99D0)*(KW_MC)**9 + (21728888.938D0)*(KW_MC)**8 +
     1  (-9865563.698D0)*(KW_MC)**7 + (2720793.638D0)*(KW_MC)**6 + (-463665.676D0)*(KW_MC)**5 +
     2  (47121.73D0)*(KW_MC)**4 + (-2540.929D0)*(KW_MC)**3 + (50.626D0)*(KW_MC)**2 + (0.037D0)*(KW_MC) + (0.377D0) 
	      KW_nu13 = (14033992.731D0)*(KW_MC)**10 + (-26688932.99D0)*(KW_MC)**9 + (21728888.938D0)*(KW_MC)**8 +
     1  (-9865563.698D0)*(KW_MC)**7 + (2720793.638D0)*(KW_MC)**6 + (-463665.676D0)*(KW_MC)**5 +
     2  (47121.73D0)*(KW_MC)**4 + (-2540.929D0)*(KW_MC)**3 + (50.626D0)*(KW_MC)**2 + (0.037D0)*(KW_MC) + (0.377D0) 
	      KW_nu12 = (14033992.731D0)*(KW_MC)**10 + (-26688932.99D0)*(KW_MC)**9 + (21728888.938D0)*(KW_MC)**8 +
     1  (-9865563.698D0)*(KW_MC)**7 + (2720793.638D0)*(KW_MC)**6 + (-463665.676D0)*(KW_MC)**5 +
     2  (47121.73D0)*(KW_MC)**4 + (-2540.929D0)*(KW_MC)**3 + (50.626D0)*(KW_MC)**2 + (0.037D0)*(KW_MC) + (0.377D0) 
	      KW_G23 = (1032489675.073D0)*(KW_MC)**10 + (-1713861472.187D0)*(KW_MC)**9 + (1224162397.86D0)*(KW_MC)**8 +
     1  (-490188214.182D0)*(KW_MC)**7 + (119882175.493D0)*(KW_MC)**6 + (-18213047.977D0)*(KW_MC)**5 +
     2  (1656132.113D0)*(KW_MC)**4 + (-79582.653D0)*(KW_MC)**3 + (1327.146D0)*(KW_MC)**2 + (3.859D0)*(KW_MC) + (2.24D0) 
	      KW_G13 = (1032489675.073D0)*(KW_MC)**10 + (-1713861472.187D0)*(KW_MC)**9 + (1224162397.86D0)*(KW_MC)**8 +
     1  (-490188214.182D0)*(KW_MC)**7 + (119882175.493D0)*(KW_MC)**6 + (-18213047.977D0)*(KW_MC)**5 +
     2  (1656132.113D0)*(KW_MC)**4 + (-79582.653D0)*(KW_MC)**3 + (1327.146D0)*(KW_MC)**2 + (3.859D0)*(KW_MC) + (2.24D0) 
	      KW_G12 = (1032489675.073D0)*(KW_MC)**10 + (-1713861472.187D0)*(KW_MC)**9 + (1224162397.86D0)*(KW_MC)**8 +
     1  (-490188214.182D0)*(KW_MC)**7 + (119882175.493D0)*(KW_MC)**6 + (-18213047.977D0)*(KW_MC)**5 +
     2  (1656132.113D0)*(KW_MC)**4 + (-79582.653D0)*(KW_MC)**3 + (1327.146D0)*(KW_MC)**2 + (3.859D0)*(KW_MC) + (2.24D0) 

	    ELSE IF (MATNAME == 'P') THEN
	      KW_E1 = (2806146977.943D0)*(KW_MC)**10 + (-4641648853.314D0)*(KW_MC)**9 + (3303116949.83D0)*(KW_MC)**8 +
     1  (-1317562279.701D0)*(KW_MC)**7 + (320958933.22D0)*(KW_MC)**6 + (-48569827.335D0)*(KW_MC)**5 +
     2  (4399632.672D0)*(KW_MC)**4 + (-210611.458D0)*(KW_MC)**3 + (3478.562D0)*(KW_MC)**2 + (14.55D0)*(KW_MC) + (5.689D0) 
	      KW_E2 = (2806146977.943D0)*(KW_MC)**10 + (-4641648853.314D0)*(KW_MC)**9 + (3303116949.83D0)*(KW_MC)**8 +
     1  (-1317562279.701D0)*(KW_MC)**7 + (320958933.22D0)*(KW_MC)**6 + (-48569827.335D0)*(KW_MC)**5 +
     2  (4399632.672D0)*(KW_MC)**4 + (-210611.458D0)*(KW_MC)**3 + (3478.562D0)*(KW_MC)**2 + (14.55D0)*(KW_MC) + (5.689D0) 
	      KW_E3 = (2806146977.943D0)*(KW_MC)**10 + (-4641648853.314D0)*(KW_MC)**9 + (3303116949.83D0)*(KW_MC)**8 +
     1  (-1317562279.701D0)*(KW_MC)**7 + (320958933.22D0)*(KW_MC)**6 + (-48569827.335D0)*(KW_MC)**5 +
     2  (4399632.672D0)*(KW_MC)**4 + (-210611.458D0)*(KW_MC)**3 + (3478.562D0)*(KW_MC)**2 + (14.55D0)*(KW_MC) + (5.689D0) 
	      KW_nu23 = (14033992.731D0)*(KW_MC)**10 + (-26688932.99D0)*(KW_MC)**9 + (21728888.938D0)*(KW_MC)**8 +
     1  (-9865563.698D0)*(KW_MC)**7 + (2720793.638D0)*(KW_MC)**6 + (-463665.676D0)*(KW_MC)**5 +
     2  (47121.73D0)*(KW_MC)**4 + (-2540.929D0)*(KW_MC)**3 + (50.626D0)*(KW_MC)**2 + (0.037D0)*(KW_MC) + (0.377D0) 
	    KW_nu13 = (14033992.731D0)*(KW_MC)**10 + (-26688932.99D0)*(KW_MC)**9 + (21728888.938D0)*(KW_MC)**8 +
     1  (-9865563.698D0)*(KW_MC)**7 + (2720793.638D0)*(KW_MC)**6 + (-463665.676D0)*(KW_MC)**5 +
     2  (47121.73D0)*(KW_MC)**4 + (-2540.929D0)*(KW_MC)**3 + (50.626D0)*(KW_MC)**2 + (0.037D0)*(KW_MC) + (0.377D0) 
	      KW_nu12 = (14033992.731D0)*(KW_MC)**10 + (-26688932.99D0)*(KW_MC)**9 + (21728888.938D0)*(KW_MC)**8 +
     1  (-9865563.698D0)*(KW_MC)**7 + (2720793.638D0)*(KW_MC)**6 + (-463665.676D0)*(KW_MC)**5 +
     2  (47121.73D0)*(KW_MC)**4 + (-2540.929D0)*(KW_MC)**3 + (50.626D0)*(KW_MC)**2 + (0.037D0)*(KW_MC) + (0.377D0) 
	      KW_G23 = (1032489675.073D0)*(KW_MC)**10 + (-1713861472.187D0)*(KW_MC)**9 + (1224162397.86D0)*(KW_MC)**8 +
     1  (-490188214.182D0)*(KW_MC)**7 + (119882175.493D0)*(KW_MC)**6 + (-18213047.977D0)*(KW_MC)**5 +
     2  (1656132.113D0)*(KW_MC)**4 + (-79582.653D0)*(KW_MC)**3 + (1327.146D0)*(KW_MC)**2 + (3.859D0)*(KW_MC) + (2.24D0) 
	      KW_G13 = (1032489675.073D0)*(KW_MC)**10 + (-1713861472.187D0)*(KW_MC)**9 + (1224162397.86D0)*(KW_MC)**8 +
     1  (-490188214.182D0)*(KW_MC)**7 + (119882175.493D0)*(KW_MC)**6 + (-18213047.977D0)*(KW_MC)**5 +
     2  (1656132.113D0)*(KW_MC)**4 + (-79582.653D0)*(KW_MC)**3 + (1327.146D0)*(KW_MC)**2 + (3.859D0)*(KW_MC) + (2.24D0) 
	      KW_G12 = (1032489675.073D0)*(KW_MC)**10 + (-1713861472.187D0)*(KW_MC)**9 + (1224162397.86D0)*(KW_MC)**8 +
     1  (-490188214.182D0)*(KW_MC)**7 + (119882175.493D0)*(KW_MC)**6 + (-18213047.977D0)*(KW_MC)**5 +
     2  (1656132.113D0)*(KW_MC)**4 + (-79582.653D0)*(KW_MC)**3 + (1327.146D0)*(KW_MC)**2 + (3.859D0)*(KW_MC) + (2.24D0) 

	    ELSE IF (MATNAME == 'S1') THEN
	      KW_E1 = (1314753273.457D0)*(KW_MC)**10 + (-2091851107.862D0)*(KW_MC)**9 + (1426366268.91D0)*(KW_MC)**8 +
     1  (-543074654.128D0)*(KW_MC)**7 + (125825894.75D0)*(KW_MC)**6 + (-18057356.105D0)*(KW_MC)**5 +
     2  (1548616.603D0)*(KW_MC)**4 + (-70240.791D0)*(KW_MC)**3 + (1134.894D0)*(KW_MC)**2 + (-7.22D0)*(KW_MC) + (50.564D0) 
	      KW_E2 = (1370908815.978D0)*(KW_MC)**10 + (-2344731684.774D0)*(KW_MC)**9 + (1725024864.19D0)*(KW_MC)**8 +
     1  (-710865284.323D0)*(KW_MC)**7 + (178658365.599D0)*(KW_MC)**6 + (-27831966.76D0)*(KW_MC)**5 +
     2  (2586450.015D0)*(KW_MC)**4 + (-126387.678D0)*(KW_MC)**3 + (2142.232D0)*(KW_MC)**2 + (-9.809D0)*(KW_MC) + (6.875D0) 
	      KW_E3 = (1370908815.978D0)*(KW_MC)**10 + (-2344731684.774D0)*(KW_MC)**9 + (1725024864.19D0)*(KW_MC)**8 +
     1  (-710865284.323D0)*(KW_MC)**7 + (178658365.599D0)*(KW_MC)**6 + (-27831966.76D0)*(KW_MC)**5 +
     2  (2586450.015D0)*(KW_MC)**4 + (-126387.678D0)*(KW_MC)**3 + (2142.232D0)*(KW_MC)**2 + (-9.809D0)*(KW_MC) + (6.875D0) 
	      KW_nu23 = (480473.676D0)*(KW_MC)**10 + (974453.341D0)*(KW_MC)**9 + (-2077836.474D0)*(KW_MC)**8 +
     1  (1429832.568D0)*(KW_MC)**7 + (-505698.335D0)*(KW_MC)**6 + (101825.505D0)*(KW_MC)**5 +
     2  (-11640.393D0)*(KW_MC)**4 + (680.969D0)*(KW_MC)**3 + (-14.22D0)*(KW_MC)**2 + (-0.264D0)*(KW_MC) + (0.422D0) 
	      KW_nu13 = (-0.0D0)*(KW_MC)**10 + (0.0D0)*(KW_MC)**9 + (-0.0D0)*(KW_MC)**8 +
     1  (0.0D0)*(KW_MC)**7 + (-0.0D0)*(KW_MC)**6 + (0.0D0)*(KW_MC)**5 +
     2  (-0.0D0)*(KW_MC)**4 + (0.0D0)*(KW_MC)**3 + (-0.0D0)*(KW_MC)**2 + (0.0D0)*(KW_MC) + (0.336D0) 
	      KW_nu12 = (-0.0D0)*(KW_MC)**10 + (0.0D0)*(KW_MC)**9 + (-0.0D0)*(KW_MC)**8 +
     1  (0.0D0)*(KW_MC)**7 + (-0.0D0)*(KW_MC)**6 + (0.0D0)*(KW_MC)**5 +
     2  (0.0D0)*(KW_MC)**4 + (-0.0D0)*(KW_MC)**3 + (0.0D0)*(KW_MC)**2 + (-0.0D0)*(KW_MC) + (0.336D0) 
	      KW_G23 = (493724296.826D0)*(KW_MC)**10 + (-845385862.904D0)*(KW_MC)**9 + (622712286.323D0)*(KW_MC)**8 +
     1  (-256955323.975D0)*(KW_MC)**7 + (64673860.53D0)*(KW_MC)**6 + (-10091699.317D0)*(KW_MC)**5 +
     2  (939678.456D0)*(KW_MC)**4 + (-46043.224D0)*(KW_MC)**3 + (784.311D0)*(KW_MC)**2 + (-3.02D0)*(KW_MC) + (2.417D0) 
	      KW_G13 = (475752172.593D0)*(KW_MC)**10 + (-811311794.433D0)*(KW_MC)**9 + (595243632.545D0)*(KW_MC)**8 +
     1  (-244680035.171D0)*(KW_MC)**7 + (61359552.969D0)*(KW_MC)**6 + (-9541719.822D0)*(KW_MC)**5 +
     2  (885649.958D0)*(KW_MC)**4 + (-43267.487D0)*(KW_MC)**3 + (734.538D0)*(KW_MC)**2 + (-3.403D0)*(KW_MC) + (2.441D0) 
	      KW_G12 = (475752172.593D0)*(KW_MC)**10 + (-811311794.433D0)*(KW_MC)**9 + (595243632.545D0)*(KW_MC)**8 +
     1  (-244680035.171D0)*(KW_MC)**7 + (61359552.969D0)*(KW_MC)**6 + (-9541719.822D0)*(KW_MC)**5 +
     2  (885649.958D0)*(KW_MC)**4 + (-43267.487D0)*(KW_MC)**3 + (734.538D0)*(KW_MC)**2 + (-3.403D0)*(KW_MC) + (2.441D0) 

	    ELSE IF (MATNAME == 'S2') THEN
	      KW_E1 = (820639024.188D0)*(KW_MC)**10 + (-1282100252.018D0)*(KW_MC)**9 + (855941378.171D0)*(KW_MC)**8 +
     1  (-318018894.205D0)*(KW_MC)**7 + (71637849.191D0)*(KW_MC)**6 + (-9956761.269D0)*(KW_MC)**5 +
     2  (823944.255D0)*(KW_MC)**4 + (-35966.175D0)*(KW_MC)**3 + (570.289D0)*(KW_MC)**2 + (-5.544D0)*(KW_MC) + (69.25D0) 
	      KW_E2 = (964440757.596D0)*(KW_MC)**10 + (-1713541733.989D0)*(KW_MC)**9 + (1307299767.767D0)*(KW_MC)**8 +
     1  (-557587067.944D0)*(KW_MC)**7 + (144742259.071D0)*(KW_MC)**6 + (-23239299.516D0)*(KW_MC)**5 +
     2  (2221215.722D0)*(KW_MC)**4 + (-111515.112D0)*(KW_MC)**3 + (1957.231D0)*(KW_MC)**2 + (-15.323D0)*(KW_MC) + (8.27D0) 
	      KW_E3 = (964440757.596D0)*(KW_MC)**10 + (-1713541733.989D0)*(KW_MC)**9 + (1307299767.767D0)*(KW_MC)**8 +
     1  (-557587067.944D0)*(KW_MC)**7 + (144742259.071D0)*(KW_MC)**6 + (-23239299.516D0)*(KW_MC)**5 +
     2  (2221215.722D0)*(KW_MC)**4 + (-111515.112D0)*(KW_MC)**3 + (1957.231D0)*(KW_MC)**2 + (-15.323D0)*(KW_MC) + (8.27D0) 
	      KW_nu23 = (4992849.429D0)*(KW_MC)**10 + (-6178932.317D0)*(KW_MC)**9 + (2771602.935D0)*(KW_MC)**8 +
     1  (-400709.983D0)*(KW_MC)**7 + (-86657.739D0)*(KW_MC)**6 + (42643.068D0)*(KW_MC)**5 +
     2  (-6664.508D0)*(KW_MC)**4 + (459.846D0)*(KW_MC)**3 + (-10.811D0)*(KW_MC)**2 + (-0.33D0)*(KW_MC) + (0.414D0) 
	      KW_nu13 = (0.0D0)*(KW_MC)**10 + (-0.0D0)*(KW_MC)**9 + (0.0D0)*(KW_MC)**8 +
     1  (-0.0D0)*(KW_MC)**7 + (0.0D0)*(KW_MC)**6 + (-0.0D0)*(KW_MC)**5 +
     2  (0.0D0)*(KW_MC)**4 + (-0.0D0)*(KW_MC)**3 + (0.0D0)*(KW_MC)**2 + (-0.0D0)*(KW_MC) + (0.337D0) 
	      KW_nu12 = (0.0D0)*(KW_MC)**10 + (-0.0D0)*(KW_MC)**9 + (0.0D0)*(KW_MC)**8 +
     1  (-0.0D0)*(KW_MC)**7 + (0.0D0)*(KW_MC)**6 + (-0.0D0)*(KW_MC)**5 +
     2  (0.0D0)*(KW_MC)**4 + (-0.0D0)*(KW_MC)**3 + (0.0D0)*(KW_MC)**2 + (-0.0D0)*(KW_MC) + (0.337D0) 
	      KW_G23 = (340223079.787D0)*(KW_MC)**10 + (-608326087.475D0)*(KW_MC)**9 + (466903734.293D0)*(KW_MC)**8 +
     1  (-200282066.319D0)*(KW_MC)**7 + (52274138.904D0)*(KW_MC)**6 + (-8437239.339D0)*(KW_MC)**5 +
     2  (810711.624D0)*(KW_MC)**4 + (-40944.552D0)*(KW_MC)**3 + (724.811D0)*(KW_MC)**2 + (-4.798D0)*(KW_MC) + (2.925D0) 
	      KW_G13 = (317061401.426D0)*(KW_MC)**10 + (-564343032.187D0)*(KW_MC)**9 + (431370844.084D0)*(KW_MC)**8 +
     1  (-184359872.867D0)*(KW_MC)**7 + (47960877.016D0)*(KW_MC)**6 + (-7718561.344D0)*(KW_MC)**5 +
     2  (739728.797D0)*(KW_MC)**4 + (-37266.458D0)*(KW_MC)**3 + (657.017D0)*(KW_MC)**2 + (-5.206D0)*(KW_MC) + (2.962D0) 
	      KW_G12 = (317061401.426D0)*(KW_MC)**10 + (-564343032.187D0)*(KW_MC)**9 + (431370844.084D0)*(KW_MC)**8 +
     1  (-184359872.867D0)*(KW_MC)**7 + (47960877.016D0)*(KW_MC)**6 + (-7718561.344D0)*(KW_MC)**5 +
     2  (739728.797D0)*(KW_MC)**4 + (-37266.458D0)*(KW_MC)**3 + (657.017D0)*(KW_MC)**2 + (-5.206D0)*(KW_MC) + (2.962D0) 

	    ELSE IF (MATNAME == 'S3') THEN
	      KW_E1 = (765388246.723D0)*(KW_MC)**10 + (-1206068079.689D0)*(KW_MC)**9 + (813351866.983D0)*(KW_MC)**8 +
     1  (-305815168.817D0)*(KW_MC)**7 + (69862013.789D0)*(KW_MC)**6 + (-9871013.467D0)*(KW_MC)**5 +
     2  (832645.868D0)*(KW_MC)**4 + (-37184.143D0)*(KW_MC)**3 + (615.281D0)*(KW_MC)**2 + (-10.061D0)*(KW_MC) + (63.481D0) 
	      KW_E2 = (736377654.629D0)*(KW_MC)**10 + (-1304919113.349D0)*(KW_MC)**9 + (992881340.658D0)*(KW_MC)**8 +
     1  (-422304448.354D0)*(KW_MC)**7 + (109302190.873D0)*(KW_MC)**6 + (-17492207.363D0)*(KW_MC)**5 +
     2  (1665427.285D0)*(KW_MC)**4 + (-83156.524D0)*(KW_MC)**3 + (1449.809D0)*(KW_MC)**2 + (-18.548D0)*(KW_MC) + (7.691D0) 
	      KW_E3 = (736377654.629D0)*(KW_MC)**10 + (-1304919113.349D0)*(KW_MC)**9 + (992881340.658D0)*(KW_MC)**8 +
     1  (-422304448.354D0)*(KW_MC)**7 + (109302190.873D0)*(KW_MC)**6 + (-17492207.363D0)*(KW_MC)**5 +
     2  (1665427.285D0)*(KW_MC)**4 + (-83156.524D0)*(KW_MC)**3 + (1449.809D0)*(KW_MC)**2 + (-18.548D0)*(KW_MC) + (7.691D0) 
	      KW_nu23 = (3901032.851D0)*(KW_MC)**10 + (-4646139.793D0)*(KW_MC)**9 + (1903401.341D0)*(KW_MC)**8 +
     1  (-154647.045D0)*(KW_MC)**7 + (-119482.527D0)*(KW_MC)**6 + (43126.293D0)*(KW_MC)**5 +
     2  (-6268.39D0)*(KW_MC)**4 + (418.217D0)*(KW_MC)**3 + (-9.578D0)*(KW_MC)**2 + (-0.339D0)*(KW_MC) + (0.409D0) 
	      KW_nu13 = (0.0D0)*(KW_MC)**10 + (-0.0D0)*(KW_MC)**9 + (0.0D0)*(KW_MC)**8 +
     1  (-0.0D0)*(KW_MC)**7 + (0.0D0)*(KW_MC)**6 + (-0.0D0)*(KW_MC)**5 +
     2  (0.0D0)*(KW_MC)**4 + (-0.0D0)*(KW_MC)**3 + (0.0D0)*(KW_MC)**2 + (-0.0D0)*(KW_MC) + (0.337D0) 
	      KW_nu12 = (0.0D0)*(KW_MC)**10 + (-0.0D0)*(KW_MC)**9 + (0.0D0)*(KW_MC)**8 +
     1  (-0.0D0)*(KW_MC)**7 + (0.0D0)*(KW_MC)**6 + (-0.0D0)*(KW_MC)**5 +
     2  (0.0D0)*(KW_MC)**4 + (-0.0D0)*(KW_MC)**3 + (0.0D0)*(KW_MC)**2 + (-0.0D0)*(KW_MC) + (0.337D0) 
	      KW_G23 = (261137032.97D0)*(KW_MC)**10 + (-465766582.896D0)*(KW_MC)**9 + (356572320.097D0)*(KW_MC)**8 +
     1  (-152546419.092D0)*(KW_MC)**7 + (39701988.442D0)*(KW_MC)**6 + (-6387877.623D0)*(KW_MC)**5 +
     2  (611477.428D0)*(KW_MC)**4 + (-30718.082D0)*(KW_MC)**3 + (539.64D0)*(KW_MC)**2 + (-5.97D0)*(KW_MC) + (2.729D0) 
	      KW_G13 = (253459696.27D0)*(KW_MC)**10 + (-447376860.686D0)*(KW_MC)**9 + (339243867.971D0)*(KW_MC)**8 +
     1  (-143885156.461D0)*(KW_MC)**7 + (37158616.7D0)*(KW_MC)**6 + (-5937611.424D0)*(KW_MC)**5 +
     2  (564943.532D0)*(KW_MC)**4 + (-28227.578D0)*(KW_MC)**3 + (492.946D0)*(KW_MC)**2 + (-6.278D0)*(KW_MC) + (2.75D0) 
	      KW_G12 = (253459696.27D0)*(KW_MC)**10 + (-447376860.686D0)*(KW_MC)**9 + (339243867.971D0)*(KW_MC)**8 +
     1  (-143885156.461D0)*(KW_MC)**7 + (37158616.7D0)*(KW_MC)**6 + (-5937611.424D0)*(KW_MC)**5 +
     2  (564943.532D0)*(KW_MC)**4 + (-28227.578D0)*(KW_MC)**3 + (492.946D0)*(KW_MC)**2 + (-6.278D0)*(KW_MC) + (2.75D0)

      ELSE IF (MATNAME == 'TEST') THEN
	      KW_E1 = 1.D0
	      KW_E2 = 1.D0
	      KW_E3 = 1.D0
	      KW_nu23 = 0.3D0
	      KW_nu13 = 0.3D0
	      KW_nu12 = 0.3D0
	      KW_G23 = 0.3846D0
	      KW_G13 = 0.3846D0
	      KW_G12 = 0.3846D0
      
      ELSE
         PRINT *, 'Error: Material "', MATNAME, '" not supported.'
         RETURN
      END IF ! End material selection

      ! Print engineering constants to log file
      !PRINT *, "CSV Data for ", MATNAME, ":"
      !PRINT *, "w,E1,E2,E3,G23,G13,G12,v12,v13,v23"
      !PRINT *, KW_MC, KW_E1, KW_E2, KW_E3, KW_G23, KW_G13, KW_G12, KW_nu12, KW_nu13, KW_nu23

      ! Calculate elastic compliance from engineering constants
      KW_Compl = (0.D0)		 
      KW_Compl(1,1) = (One/KW_E1)
      KW_Compl(1,2) = -(KW_nu12/KW_E1)
      KW_Compl(1,3) = -(KW_nu13/KW_E1)
      KW_Compl(2,1) = -(KW_nu12/KW_E1)
      KW_Compl(2,2) = (One/KW_E2)
      KW_Compl(2,3) = -(KW_nu23/KW_E2)
      KW_Compl(3,1) = -(KW_nu13/KW_E1)
      KW_Compl(3,2) = -(KW_nu23/KW_E2)
      KW_Compl(3,3) = (One/KW_E3)
      KW_Compl(4,4) = (One/KW_G12)
      KW_Compl(5,5) = (One/KW_G13)
      KW_Compl(6,6) = (One/KW_G23)

      RETURN
      End
! ###########################################################################################

! ###########################################################################################
! ############################    Viscoelastic compliance tensor   ##########################
! ###########################################################################################
      SUBROUTINE KW_ViscoelastParam(KW_MC, MATNAME, KW_J0L, KW_Coeff, KW_Gamma_i, KW_Tau_i)
      Implicit none	
	  
      ! Input parameters
      DOUBLE PRECISION, INTENT(IN) :: KW_MC
      CHARACTER*(*), INTENT(IN) :: MATNAME
      DOUBLE PRECISION, INTENT(IN) :: KW_J0L
      DOUBLE PRECISION, INTENT(IN) :: KW_Coeff

      ! Output: Viscoelastic parameters
      Double precision, Dimension(4), Intent(OUT) :: KW_Gamma_i
      Double precision, Dimension(4), Intent(OUT) :: KW_Tau_i

      ! Complicane scaling factors 
      INCLUDE 'layers_gamma_values.inc'
      IF (KW_Coeff > 1) THEN
        KW_Gamma_i(1) = KW_Coeff
        KW_Gamma_i(2) = KW_Coeff
        KW_Gamma_i(3) = KW_Coeff
        KW_Gamma_i(4) = KW_Coeff
      END IF

      ! Retardation times in [h] 	 	 
      KW_Tau_i(1) = (0.1D0)
      KW_Tau_i(2) = (1.D0)
      KW_Tau_i(3) = (10.D0)
      KW_Tau_i(4) = (100.D0)

      RETURN
      End
! ###########################################################################################

! ###########################################################################################
! #########################    Viscoelastic time functions Tn & Tn+1   ######################
! ###########################################################################################
      SUBROUTINE KW_TimeFun(KW_T1Out, KW_T0Out, KW_DtimeIn, KW_TauIn) 
      Implicit None

      ! Input parameters
      Double precision :: KW_Zi
      ! Output: time functions
      Double precision,Intent(OUT) :: KW_T1Out,KW_T0Out 
      Double precision,Intent(IN) :: KW_DtimeIn,KW_TauIn
      
      KW_Zi = (KW_DtimeIn/KW_TauIn) ! xi = Dtime/tau
      KW_T1Out = (0.D0)
      KW_T0Out = (0.D0)
      KW_T1Out = (1.D0)-(((1.D0)/KW_Zi)*((1.D0)-Exp(-KW_Zi))) ! Tn+1
      KW_T0Out = (1.D0)-Exp(-KW_Zi)-KW_T1Out ! Tn

      RETURN
      END
! ###########################################################################################

! ###########################################################################################
! ##############################    Elastic residual vector (R1)   ##########################
! ###########################################################################################
      SUBROUTINE KW_ResR1(KW_R1Out, KW_Sigma_New, STRAN_In,
     1    KW_HygStr_In, KW_VisStr_In, KW_Compl_1In)
      Implicit None

      Integer i1
      Double precision, Dimension(6,1),Intent(OUT) :: KW_R1Out ! Elastic residuals = R1
      Double precision, Dimension(6,1),Intent(IN) :: KW_Sigma_New ! Stress
      Double precision, Dimension(6),Intent(IN) :: STRAN_In ! Total strain
      Double precision, Dimension(6),Intent(IN) :: KW_HygStr_In ! Hygroelastic strain
      Double precision, Dimension(4,6),Intent(IN) :: KW_VisStr_In ! Viscoelastic strain components
      Double precision, Dimension(6,6),Intent(IN) :: KW_Compl_1In ! Elastic compliance
      Double precision, Dimension(6,1) :: KW_C1S1 ! Elastic strain
      Double precision, Dimension(6) :: KW_Compl1S1 ! Elastic strain
      Double precision, Dimension(6) :: KW_SumVisStr ! Viscoelastic strain
      Double precision, Dimension(6) :: KW_R1_6 ! Residual vector

      KW_C1S1 = (0.D0)
      KW_C1S1 = Matmul(KW_Compl_1In, KW_Sigma_New)
      KW_Compl1S1(1:6) = KW_C1S1(1:6,1) ! Elastic strain

      Do i1=1, 6
	     KW_SumVisStr(i1) = Sum(KW_VisStr_In(1:4,i1)) ! Viscoelastic strain
      End Do

      KW_R1_6 = (0.D0)
      KW_R1_6 = KW_Compl1S1-(STRAN_In-KW_HygStr_In-KW_SumVisStr) 
      KW_R1Out = (0.D0)
      KW_R1Out(1:6,1) = KW_R1_6(1:6) ! Residual vector

      RETURN
      END
! ###########################################################################################

! ###########################################################################################
! ############################    Viscoelastic residual vector (R2)   #######################
! ###########################################################################################
      SUBROUTINE KW_ResR2(KW_R2iOut, KW_VEStr_OldIn, KW_VEStr_NewIn,
     1	  KW_Stress_OldIn, KW_Stress_NewIn, KW_Compl_0In, KW_Compl_1In, 
     2	  KW_Gamma_1In, KW_Gamma_0In, KW_T1In,KW_T0In, DTIME_In,
     3    KW_TauIn) 
  	  Implicit None

      Double precision,Intent(IN) :: KW_Gamma_1In ! KV parameters (gamma) at t(n+1)
      Double precision,Intent(IN) :: KW_Gamma_0In ! KV parameters (gamma) at t(n)
      Double precision,Intent(IN) :: KW_T1In, KW_T0In, DTIME_In ! Time functions at t(n) and t(n+1), and time increment
      Double precision,Intent(IN) :: KW_TauIn ! Retardations times (tau)
      Double precision,Dimension(6,1),Intent(OUT) :: KW_R2iOut ! Viscoelastic residuals = R2                   	  
      Double precision,Dimension(1,6),Intent(IN) :: KW_VEStr_OldIn ! Viscoelastic strain at t(n)
      Double precision,Dimension(1,6),Intent(IN) :: KW_VEStr_NewIn ! Viscoelastic strain at t(n+1)
      Double precision,Dimension(6),Intent(IN) :: KW_Stress_OldIn ! Stress at t(n)
      Double precision,Dimension(6,1),Intent(IN) :: KW_Stress_NewIn ! Stress at t(n+1)
      Double precision,Dimension(6,6),Intent(IN) :: KW_Compl_0In ! Elastic compliance at t(n)
      Double precision,Dimension(6,6),Intent(IN) :: KW_Compl_1In ! Elsatic compliance at t(n+1)
      Double precision,Dimension(6,1) :: KW_VEStr0,KW_VEStr1
      Double precision,Dimension(6,1) :: KW_TermV1,KW_TermV2
      Double precision,Dimension(6,1) :: KW_TermV3,KW_TermV4
      Double precision,Dimension(6,1) :: KW_S0
      Integer j1,j2
	  
      KW_VEStr0 = (0.D0)
      KW_VEStr1 = (0.D0)
		 
      Do j1=1, 6	  
         KW_VEStr0(j1,1) = KW_VEStr_OldIn(1,j1)
         KW_VEStr1(j1,1) = KW_VEStr_NewIn(1,j1)
      End Do
         KW_S0 = (0.D0)
      Do j2=1, 6		 
         KW_S0(j2,1) = KW_Stress_OldIn(j2)
      End Do		 
	  
      KW_TermV1 = (0.D0)
      KW_TermV1 = KW_VEStr0*(Exp(-DTIME_In/KW_TauIn))
      KW_TermV2 = (0.D0)
      KW_TermV2 = -KW_VEStr1
      KW_TermV3 = (0.D0)
      KW_TermV3 = (KW_T0In/KW_Gamma_0In)*Matmul(KW_Compl_0In,KW_S0) 
      KW_TermV4 = (0.D0)
      KW_TermV4 = (KW_T1In/KW_Gamma_1In)*
     1	  Matmul(KW_Compl_1In,KW_Stress_NewIn)
      KW_R2iOut = (0.D0)
      KW_R2iOut = (KW_TermV1+KW_TermV2+KW_TermV3+KW_TermV4)
	  
      RETURN
      END
! ###########################################################################################

! ###########################################################################################
! ################################    Viscoelastic strain change   ##########################
! ###########################################################################################
      SUBROUTINE KW_DeltaVE(KW_DVEiOut, KW_R2iIn, KW_Compl_1In,
     1	  KW_DStressIn, KW_Gamma_1In, KW_T1In)
      Implicit None

      Double precision,Intent(IN) :: KW_Gamma_1In ! KV parameters at t(n+1)
      Double precision,Intent(IN) :: KW_T1In ! Time function at t(n+1)
      Double precision,Dimension(1,6),Intent(OUT):: KW_DVEiOut ! Viscoelastic strain change
      Double precision,Dimension(6,1),Intent(IN) :: KW_R2iIn ! Viscoelastic residual vector
      Double precision,Dimension(6,6),Intent(IN) :: KW_Compl_1In ! Elastic compliance at t(n+1)
      Double precision,Dimension(6,1),Intent(IN) :: KW_DStressIn ! Stress change
      Double precision,Dimension(6,1) :: KW_Term1,KW_Term2
      Double precision,Dimension(6,1) :: KW_DVEi
      Integer i1
	  
      KW_Term1 = (0.D0)
      KW_Term1 = KW_R2iIn
	  
      KW_Term2 = (0.D0)
      KW_Term2 = (KW_T1In/KW_Gamma_1In)*Matmul(KW_Compl_1In,KW_DStressIn)
      KW_DVEi = (0.D0)
      KW_DVEi = (KW_Term1+KW_Term2)
      Do i1=1, 6
         KW_DVEiOut(1,i1) = KW_DVEi(i1,1)
      End Do		 

      RETURN
      END
! ###########################################################################################

! ===========================================================================================