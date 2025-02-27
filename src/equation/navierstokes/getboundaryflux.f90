!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2020 Andrés Rueda
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/flexi)
!
! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "defines.h"

!==================================================================================================================================
!> Routines to provide boundary conditions for the domain. Fills the boundary part of the fluxes list.
!==================================================================================================================================
MODULE MOD_GetBoundaryFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitBC
  MODULE PROCEDURE InitBC
END INTERFACE

INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE

#if PARABOLIC
INTERFACE Lifting_GetBoundaryFlux
  MODULE PROCEDURE Lifting_GetBoundaryFlux
END INTERFACE
#endif /*PARABOLIC*/

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE

PUBLIC::InitBC
PUBLIC::GetBoundaryFlux
#if PARABOLIC
PUBLIC::Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/
PUBLIC::FinalizeBC
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Initialize boundary conditions
!==================================================================================================================================
SUBROUTINE InitBC()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: nBCByType,BCSideID,BCdata
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,nBCs,BoundaryType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iBC,SideID,BCType
LOGICAL :: fillBC
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
   CALL abort(__STAMP__,&
     "InitBC not ready to be called or already called.")
END IF

! Count number of sides of each boundary
ALLOCATE(nBCByType(nBCs))
nBCByType=0
DO SideID=1,nBCSides
  DO iBC=1,nBCs
    IF(BC(SideID).EQ.iBC) nBCByType(iBC)=nBCByType(iBC)+1
  END DO
END DO

! Sort BCs by type, store SideIDs
ALLOCATE(BCSideID(nBCs,MAXVAL(nBCByType)))
nBCByType=0
DO SideID=1,nBCSides
  DO iBC=1,nBCs
    IF(BC(SideID).EQ.iBC)THEN
      nBCByType(iBC)=nBCByType(iBC)+1
      BCSideID(iBC,nBCByType(iBC))=SideID
    END IF
  END DO
END DO

fillBC=.FALSE.
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  SELECT CASE(BCType)
  CASE(20,220)
    fillBC=.TRUE.
  END SELECT
END DO

IF(fillBC)THEN
  ! Allocate buffer array to store temp data for all BC sides
  ALLOCATE(BCData(PP_nVar,0:PP_N,0:PP_N,nBCSides))
  BCData=0.
  
  ! Fill  BC data for steady BCs (BCtype = 20, 22)
  CALL FillBCdata(0.,BCdata)
END IF !fillBC

END SUBROUTINE InitBC


!==================================================================================================================================
!> precompute and store steadyState BC data for BCtype 20,220 (Dirichlet BC like 2,22) 
!==================================================================================================================================
SUBROUTINE FillBCdata(tIn,BCdata)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: Face_xGP
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: tIn                                  !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: BCdata(PP_nVar,0:PP_N,0:PP_N,nBCSides) !<boundary date filled
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
!==================================================================================================================================
! some BCdata is not changing over time, store into BCdata
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(20)
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    IF(BCState.EQ.0)THEN
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
            CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),BCdata(:,p,q,SideID))
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    ELSE !BCstate /= 0
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
            BCdata(:,p,q,SideID) = RefStateCons(BCState,:)
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    END IF !BCState=0
  CASE(220) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),BCdata(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  END SELECT ! BCType
END DO !iBC=1,nBCs
END SUBROUTINE FillBCdata


!==================================================================================================================================
!> Computes the boundary values for the navierstokes part
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Riemann      ,ONLY: Riemann,AdvRiemann
USE MOD_DG_Vars      ,ONLY: U_Master
USE MOD_Mesh_Vars    ,ONLY: nSides,nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: ConstoPrim_aux,ConstoPrim,PrimtoCons
USE MOD_Equation_Vars,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1,RefStatePrim,RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
USE MOD_Flux         ,ONLY: EvalEulerFlux1D
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: gradPx_Master,gradPy_Master,gradPz_Master
USE MOD_Flux         ,ONLY: EvalDiffFlux3D,EvalDiffFlux1D_Outflow
#endif /*PARABOLIC*/
#if FLUXO_HYPERSONIC
USE MOD_NFVSE_Vars   ,ONLY: alpha_vis
USE MOD_Mesh_Vars    ,ONLY: SideToElem
#endif
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseGetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: tIn                                 !<current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nSides)  !<Navierstokes boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,iVar,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N),N_loc(1:3,1:3)
REAL                                 :: Prim(1:8),ar,br
REAL                                 :: U_loc(PP_nVar)
real                                 :: pres, a, normalMachNo
#if PARABOLIC
REAL                                 :: BCGradMat(1:3,1:3)
REAL                                 :: Fd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: Gd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: Hd_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradPx_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradPy_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradPz_Face_loc(1:PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: gradPn_loc(PP_nVar),gradPt1_loc(PP_nVar),gradPt2_loc(PP_nVar)
REAL                                 :: gradVel(3,0:PP_N,0:PP_N)
#endif /*PARABOLIC*/
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(20,220) !steadyStateBCs
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),BCdata(:,:,:,SideID), &
#if PARABOLIC
                   gradPx_Master(:,:,:,SideID),gradPx_Master(:,:,:,SideID), &
                   gradPy_Master(:,:,:,SideID),gradPy_Master(:,:,:,SideID), &
                   gradPz_Master(:,:,:,SideID),gradPz_Master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID) &
#if FLUXO_HYPERSONIC
                   ! use alpha_vis directly because alpha_vis_Master is not prolonged to BC sides
                   , alpha_vis(SideToElem(S2E_ELEM_ID,SideID)) &
#endif
                   )
    END DO !iSide=1,nBCloc
  CASE(2)
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    IF(BCState.EQ.0)THEN
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
            CALL ExactFunc(IniExactFunc,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
          END DO ! p
        END DO ! q
        CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                     gradPx_Master(:,:,:,SideID),gradPx_Master(:,:,:,SideID), &
                     gradPy_Master(:,:,:,SideID),gradPy_Master(:,:,:,SideID), &
                     gradPz_Master(:,:,:,SideID),gradPz_Master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                     NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID) &
#if FLUXO_HYPERSONIC
                   , alpha_vis(SideToElem(S2E_ELEM_ID,SideID)) &
#endif
                   )
      END DO !iSide=1,nBCloc
    ELSE !BCstate /= 0
      DO q=0,PP_N
        DO p=0,PP_N
          U_Face_loc(:,p,q) = RefStateCons(BCState,:)
        END DO ! p
      END DO ! q
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                     gradPx_Master(:,:,:,SideID),gradPx_Master(:,:,:,SideID), &
                     gradPy_Master(:,:,:,SideID),gradPy_Master(:,:,:,SideID), &
                     gradPz_Master(:,:,:,SideID),gradPz_Master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                     NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID) &
#if FLUXO_HYPERSONIC
                   , alpha_vis(SideToElem(S2E_ELEM_ID,SideID)) &
#endif
                   )
      END DO !iSide=1,nBCloc
    END IF !BCState=0
  
  CASE(22) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
        END DO ! p
      END DO ! q
      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                   gradPx_Master(:,:,:,SideID),gradPx_Master(:,:,:,SideID), &
                   gradPy_Master(:,:,:,SideID),gradPy_Master(:,:,:,SideID), &
                   gradPz_Master(:,:,:,SideID),gradPz_Master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID) &
#if FLUXO_HYPERSONIC
                   , alpha_vis(SideToElem(S2E_ELEM_ID,SideID)) &
#endif
                   )
    END DO !iSide=1,nBCloc
  CASE(3) ! Adiabatic Wall
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ! ghost boundary state: reverse velocity, everything else same
          U_Face_loc(1,p,q) = U_Master(1,p,q,SideID)
          U_Face_loc(5,p,q) = U_Master(5,p,q,SideID)
          U_Face_loc(2:4,p,q) = -U_Master(2:4,p,q,SideID)
          CALL ConsToPrim_aux(Prim, U_Master(:,p,q,SideID))
#if PARABOLIC
          ! all gradients except temperature equal
          gradPx_Face_loc(2:5,p,q) = gradPx_Master(2:5,p,q,SideID)
          gradPx_Face_loc(1,p,q) = gradPx_Face_loc(5,p,q)/prim(7)*sKappaM1 ! (grad rho) = (grad p)/RT
          gradPy_Face_loc(2:5,p,q) = gradPy_Master(2:5,p,q,SideID)
          gradPy_Face_loc(1,p,q) = gradPy_Face_loc(5,p,q)/prim(7)*sKappaM1
          gradPz_Face_loc(2:5,p,q) = gradPz_Master(2:5,p,q,SideID)
          gradPz_Face_loc(1,p,q) = gradPz_Face_loc(5,p,q)/prim(7)*sKappaM1
#endif
        END DO ! p
      END DO ! q
      CALL AdvRiemann(Flux(:,:,:,SideId),U_Master(:,:,:,SideId),U_Face_loc,NormVec(:,:,:,SideID), &
        TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
#if PARABOLIC
      ! reset U_Face_loc to now have zero velocity (earlier it was ghost state)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ConsToPrim(Prim, U_Master(:,p,q,SideID))
          U_Face_loc(1,p,q) = U_Master(1,p,q,SideID)
          U_Face_loc(5,p,q) = Prim(5)*skappaM1
          U_Face_loc(2:4,p,q) = 0.0
        END DO
      END DO
      CALL EvalDiffFlux3D(Fd_Face_loc,Gd_Face_loc,Hd_Face_loc,U_Face_loc(:,:,:), &
                          gradPx_Face_loc,                                       &
                          gradPy_Face_loc,                                       &
                          gradPz_Face_loc)                           
      ! Sum up Euler and Diffusion Flux
      DO iVar=2,PP_nVar
#if FLUXO_HYPERSONIC
        Flux(iVar,:,:,SideID) = Flux(iVar,:,:,SideID) + (1-alpha_vis(SideToElem(S2E_ELEM_ID,SideID)))*( & 
                                NormVec(1,:,:,SideID)*Fd_Face_loc(iVar,:,:) + &
                                NormVec(2,:,:,SideID)*Gd_Face_loc(iVar,:,:) + &
                                NormVec(3,:,:,SideID)*Hd_Face_loc(iVar,:,:))
#else
        Flux(iVar,:,:,SideID) = Flux(iVar,:,:,SideID)               + & 
                                NormVec(1,:,:,SideID)*Fd_Face_loc(iVar,:,:) + &
                                NormVec(2,:,:,SideID)*Gd_Face_loc(iVar,:,:) + &
                                NormVec(3,:,:,SideID)*Hd_Face_loc(iVar,:,:)
#endif
      END DO ! ivar
#endif /*PARABOLIC*/
    END DO !iSide=1,nBCLoc
  CASE(4) ! Isothermal Wall
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ! set the ghost state: reverse velocity and everything else same
          U_Face_loc(1,p,q) = U_Master(1,p,q,SideID)
          U_Face_loc(5,p,q) = U_Master(5,p,q,SideID)
          U_Face_loc(2:4,p,q) = -U_Master(2:4,p,q,SideID)
        END DO ! p
      END DO ! q
      CALL Riemann(Flux(:,:,:,SideID),U_Master(:,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                   gradPx_Master(:,:,:,SideID),gradPx_Master(:,:,:,SideID), &
                   gradPy_Master(:,:,:,SideID),gradPy_Master(:,:,:,SideID), &
                   gradPz_Master(:,:,:,SideID),gradPz_Master(:,:,:,SideID), &
#endif /*PARABOLIC*/
                   NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID) &
#if FLUXO_HYPERSONIC
                   , alpha_vis(SideToElem(S2E_ELEM_ID,SideID)) &
#endif
                   )
    END DO !iSide=1,nBCLoc
  
  CASE(9) ! Euler Wall, slip wall, symmetry BC, v=0 strategy a la HALO (is very perfect)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          U_loc(1) = U_Master(1,p,q,SideID)
          U_loc(5) = U_Master(5,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! rotate momentum in normal direction
          U_loc(2)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,1))
          U_loc(3)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,2))
          U_loc(4)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,3))
          ! Compute the primitives
          CALL ConsToPrim_aux(Prim,U_loc)
          ! Compute the 1D wall Riemann problem pressure solution
          IF (Prim(2) .LE. 0.) THEN
            Prim(5)=Prim(5)*max((1.+0.5*KappaM1*Prim(2)/Prim(6)),0.0001)**(2.*Kappa*sKappaM1)
          ELSE
            ar=2.*sKappaP1/Prim(1)
            br=KappaM1*sKappaP1*Prim(5)
            Prim(5)=Prim(5)+Prim(2)/ar*0.5*(Prim(2)+sqrt(Prim(2)*Prim(2)+4.*ar*(Prim(5)+br)))
          END IF
          ! Now we compute the 1D Euler flux, but use the info that the normal component u=0
          ! we directly tranform the flux back into the Cartesian coords: F=(0,n1*p,n2*p,n3*p,0)^T
          Flux(  1,p,q,SideID) = 0.
          Flux(2:4,p,q,SideID) = Prim(5)*N_loc(:,1)
          Flux(  5,p,q,SideID) = 0.
  
#if PARABOLIC
          !! Diffusion
          ! We prepare the gradients and set the normal derivative to zero (symmetry condition!)
          ! BCGradMat = I - n * n^T = (gradient -normal component of gradient)
          BCGradMat(1,1) = 1. - N_loc(1,1)*N_loc(1,1)
          BCGradMat(2,2) = 1. - N_loc(2,1)*N_loc(2,1)
          BCGradMat(3,3) = 1. - N_loc(3,1)*N_loc(3,1)
          BCGradMat(1,2) = -N_loc(1,1)*N_loc(2,1)
          BCGradMat(1,3) = -N_loc(1,1)*N_loc(3,1)
          BCGradMat(3,2) = -N_loc(3,1)*N_loc(2,1)
          BCGradMat(2,1) = BCGradMat(1,2)
          BCGradMat(3,1) = BCGradMat(1,3)
          BCGradMat(2,3) = BCGradMat(3,2)
          gradPx_Face_loc(:,p,q) = BCGradMat(1,1) * gradPx_Master(:,p,q,SideID) &
                                 + BCGradMat(1,2) * gradPy_Master(:,p,q,SideID) &
                                 + BCGradMat(1,3) * gradPz_Master(:,p,q,SideID)
          gradPy_Face_loc(:,p,q) = BCGradMat(2,1) * gradPx_Master(:,p,q,SideID) &
                                 + BCGradMat(2,2) * gradPy_Master(:,p,q,SideID) &
                                 + BCGradMat(2,3) * gradPz_Master(:,p,q,SideID)
          gradPz_Face_loc(:,p,q) = BCGradMat(3,1) * gradPx_Master(:,p,q,SideID) &
                                 + BCGradMat(3,2) * gradPy_Master(:,p,q,SideID) &
                                 + BCGradMat(3,3) * gradPz_Master(:,p,q,SideID)
#endif /*PARABOLIC*/
        END DO ! p
      END DO ! q
#if PARABOLIC
      ! Evaluate 3D Diffusion Flux with interior state and symmetry gradients
      CALL EvalDiffFlux3D(Fd_Face_loc,Gd_Face_loc,Hd_Face_loc,U_Master(:,:,:,SideID), &
                          gradPx_Face_loc,gradPy_Face_loc,gradPz_Face_loc)
      ! Sum up Euler and Diffusion Flux
      DO iVar=2,PP_nVar
#if FLUXO_HYPERSONIC
        Flux(iVar,:,:,SideID) = Flux(iVar,:,:,SideID) + (1-alpha_vis(SideToElem(S2E_ELEM_ID,SideID)))*( & 
                                NormVec(1,:,:,SideID)*Fd_Face_loc(iVar,:,:) + &
                                NormVec(2,:,:,SideID)*Gd_Face_loc(iVar,:,:) + &
                                NormVec(3,:,:,SideID)*Hd_Face_loc(iVar,:,:))
#else
        Flux(iVar,:,:,SideID) = Flux(iVar,:,:,SideID)                       + & 
                                NormVec(1,:,:,SideID)*Fd_Face_loc(iVar,:,:) + &
                                NormVec(2,:,:,SideID)*Gd_Face_loc(iVar,:,:) + &
                                NormVec(3,:,:,SideID)*Hd_Face_loc(iVar,:,:)
#endif
      END DO ! ivar
#endif /*PARABOLIC*/
    END DO !iSide=1,nBCLoc
  CASE(10) ! outflow with gradient zero and solution intern
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ! Get inner state
          U_loc = U_Master(:,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! transform state into normal system
          U_Face_loc(1,p,q)= U_loc(1)
          U_Face_loc(5,p,q)= U_loc(5)
          U_Face_loc(2,p,q)= SUM(U_loc(2:4)*N_loc(:,1))
          U_Face_loc(3,p,q)= SUM(U_loc(2:4)*N_loc(:,2))
          U_Face_loc(4,p,q)= SUM(U_loc(2:4)*N_loc(:,3))
          
          ! get pressure from refstate if subsonic
          ! ignored if BCRefState is zero
          if (BCState .ne. 0) then
            pres = KappaM1*( U_loc(5) - 0.5 *(U_loc(2)**2 + U_loc(3)**2 + U_loc(4)**2)/U_loc(1) )
            a    = sqrt(Kappa*pres/U_loc(1))
            normalMachNo = ABS(U_Face_loc(2,p,q)/(a*U_loc(1)))
            if (normalMachNo<=1.0) then
              CALL ConsToPrim(Prim(1:5),U_Face_loc(:,p,q))
              Prim(5) = RefStatePrim(BCState,5)
              ! U_loc contains now the state with pressure from outside (refstate)
              CALL PrimToCons(Prim(1:5),U_Face_loc(:,p,q))
            end if
          end if
#if PARABOLIC
          ! for diffusion, we use the rotational invariance of the diffusion fluxes
          !   for this, we need to transform the gradients into the normal system (see GG Diss for details)
          !    two steps: tranform the derivatives (dx,dy,dz) into normal system, then compute normal and tangential components
          gradVel(1,p,q) = SUM((  N_loc(1,1)*gradPx_Master(2:4,p,q,SideID) &
                                + N_loc(2,1)*gradPy_Master(2:4,p,q,SideID) &
                                + N_loc(3,1)*gradPz_Master(2:4,p,q,SideID))*N_loc(:,1))
          gradVel(2,p,q) = SUM((  N_loc(1,2)*gradPx_Master(2:4,p,q,SideID) &
                                + N_loc(2,2)*gradPy_Master(2:4,p,q,SideID) &
                                + N_loc(3,2)*gradPz_Master(2:4,p,q,SideID))*N_loc(:,2))
          gradVel(3,p,q) = SUM((  N_loc(1,3)*gradPx_Master(2:4,p,q,SideID) &
                                + N_loc(2,3)*gradPy_Master(2:4,p,q,SideID) &
                                + N_loc(3,3)*gradPz_Master(2:4,p,q,SideID))*N_loc(:,3))
#endif /*PARABOLIC*/
        END DO ! p
      END DO ! q
      ! Compute 1D Euler Flux Fe_Face_loc
      CALL EvalEulerFlux1D(U_Face_loc,Flux(:,:,:,SideID))
#if PARABOLIC
      ! Compute 1D diffusion Flux Fd_Face_loc
      !   Use: tau_12 = 0, tau_13 = 0, q1 = 0 (Paper POINSOT and LELE, JCP 1992, page 113, Table IV)
      !   We use special evalflux routine
      CALL EvalDiffFlux1D_Outflow (Fd_Face_loc,U_Face_loc,gradVel)
      ! Sum up Euler and Diffusion Flux and tranform back into Cartesian system
#if FLUXO_HYPERSONIC
      Flux(:,:,:,SideID) = Flux(:,:,:,SideID) + (1-alpha_vis(SideToElem(S2E_ELEM_ID,SideID)))*Fd_Face_loc
#else
      Flux(:,:,:,SideID) = Flux(:,:,:,SideID) + Fd_Face_loc
#endif /*FLUXO_HYPERSONIC*/
#endif /*PARABOLIC*/
      DO q=0,PP_N
        DO p=0,PP_N
          Flux(2:4,p,q,SideID)=NormVec( :,p,q,SideID)*Flux(2,p,q,SideID) &
                              +TangVec1(:,p,q,SideID)*Flux(3,p,q,SideID) &
                              +TangVec2(:,p,q,SideID)*Flux(4,p,q,SideID)
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCLoc
  CASE DEFAULT !  check for BCtypes in Testcase
    CALL TestcaseGetBoundaryFlux(iBC,tIn,Flux)
  END SELECT ! BCType
END DO !iBC=1,nBCs
! Integrate over the surface
DO SideID=1,nBCSides
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO
  END DO
END DO! SideID
END SUBROUTINE GetBoundaryFlux


#if PARABOLIC
!==================================================================================================================================
!> Computes the boundary fluxes for the lifting procedure for all sides.
!==================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux(tIn,Flux)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_DG_Vars      ,ONLY: U_Master
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: ConsToPrim_aux,ConsToPrim,PrimToCons,ConsToEntropy
USE MOD_Equation_Vars,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1,RefStatePrim,RefStateCons
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID,BCData
USE MOD_Testcase_GetBoundaryFlux, ONLY: TestcaseLiftingGetBoundaryFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: tIn       !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nBCSides) !<lifting boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: ar,br,N_loc(1:3,1:3),Prim(8)
REAL                                 :: U_loc(PP_nVar)
REAL                                 :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: pres, a, normalMachNo
REAL                                 :: F(PP_nVar), F_m(PP_nVar)
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(20,220)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          Flux(:,p,q,SideID)=BCdata(:,p,q,SideID)
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(2)
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    ! Dirichlet means that we use the gradients from inside the grid cell
    ! BR1 uses arithmetic mean value of states for the Riemann flux
    IF(BCState.EQ.0)THEN
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
            CALL ExactFunc(IniExactFunc,tIn,face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    ELSE !BCstate /=0
      ! use BCState as refstate number
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        DO q=0,PP_N
          DO p=0,PP_N
              Flux(:,p,q,SideID)=RefStateCons(BCState,:)
          END DO ! p
        END DO ! q
      END DO !iSide=1,nBCloc
    END IF !BCstate=0
  
  CASE(22) ! exact BC = Dirichlet BC !!
    ! SPECIAL BC: BCState specifies exactfunc to be used!!
    ! Dirichlet means that we use the gradients from inside the grid cell
    ! BR1 uses arithmetic mean value of states for the Riemann flux
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCState,tIn,Face_xGP(:,p,q,SideID),Flux(:,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  
  CASE(3) ! Adiabatic Wall, Density and pressure from interior, velocity=0
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ! pressure and density from inside, velocity 0
          Flux(2:4,p,q,SideID) = 0
          CALL ConsToPrim(Prim, U_Master(:,p,q,SideID))
          Flux(1,p,q,SideID) = Prim(1)
          Flux(5,p,q,SideID) = Prim(5)*sKappaM1
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  
  CASE(4) ! Isothermal Wall
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          ! pressure from inside, temperature from refstate, velocity 0
          Flux(2:4,p,q,SideID) = 0
          CALL ConsToPrim(Prim, U_Master(:,p,q,SideID))
          Flux(1,p,q,SideID) = Prim(5)/RefStatePrim(BCState,5)*RefStatePrim(BCState,1)
          Flux(5,p,q,SideID) = Prim(5)*sKappaM1
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(9) ! Euler Wall, slip wall, symmetry BC, v=0 strategy a la HALO (is very perfect)
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          u_loc(1) = U_Master(1,p,q,SideID)
          u_loc(5) = U_Master(5,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! rotate momentum in normal direction
          U_loc(2)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,1))
          U_loc(3)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,2))
          U_loc(4)   = SUM(U_Master(2:4,p,q,SideID)*N_loc(:,3))
          ! Compute the primitives
          CALL ConsToPrim_aux(Prim,u_loc)
          ! Compute the 1D wall Riemann problem pressure solution
          IF (Prim(2) .LE. 0.) THEN
            Prim(5)=Prim(5)*max((1.+0.5*KappaM1*Prim(2)/Prim(6)),0.0001)**(2.*Kappa*sKappaM1)
          ELSE
            ar=2.*sKappaP1/Prim(1)
            br=KappaM1*sKappaP1*Prim(5)
            Prim(5)=Prim(5)+Prim(2)/ar*0.5*(Prim(2)+sqrt(Prim(2)*Prim(2)+4.*ar*(Prim(5)+br)))
          END IF
          Prim(2)=0.
          CALL PrimToCons(prim(1:5),u_loc)
          ! Compute Flux
          Flux(1,p,q,SideID)=u_loc(1)
          Flux(5,p,q,SideID)=0.5*(u_loc(5)+U_Master(5,p,q,SideID))
          ! Rotate momentum back in x,y,z coords
          Flux(2:4,p,q,SideID) = 0.5*(u_loc(2)*N_loc(:,1)+u_loc(3)*N_loc(:,2)+u_loc(4)*N_loc(:,3) &
                                 +U_Master(2:4,p,q,SideID))
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE(10) ! outflow with gradient zero and solution intern
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          U_Loc = U_Master(:,p,q,SideID)
          ! local normal system
          N_loc(:,1) = NormVec( :,p,q,SideID)
          N_loc(:,2) = TangVec1(:,p,q,SideID)
          N_loc(:,3) = TangVec2(:,p,q,SideID)
          ! transform state into normal system
          U_Face_loc(1,p,q)= U_loc(1)
          U_Face_loc(5,p,q)= U_loc(5)
          U_Face_loc(2,p,q)= SUM(U_loc(2:4)*N_loc(:,1))
          U_Face_loc(3,p,q)= SUM(U_loc(2:4)*N_loc(:,2))
          U_Face_loc(4,p,q)= SUM(U_loc(2:4)*N_loc(:,3))
          CALL ConsToPrim(Prim(1:5),U_Loc)
          ! get pressure from refstate if subsonic
          ! ignored if BCRefState is zero
          if (BCState .ne. 0) then
            pres = KappaM1*( U_Loc(5) - 0.5 *(U_Loc(2)**2 + U_Loc(3)**2 + U_Loc(4)**2)/U_Loc(1) )
            a    = sqrt(Kappa*pres/U_Loc(1))
            normalMachNo = ABS(U_Face_loc(2,p,q)/(a*U_Loc(1)))
            if (normalMachNo<=1.0) then
              CALL ConsToPrim(Prim(1:5),U_Face_loc(:,p,q))
              Prim(5) = RefStatePrim(BCState,5)
              ! U_Loc contains now the state with pressure from outside (refstate)
              CALL PrimToCons(Prim(1:5),U_Face_loc(:,p,q))
            end if
          end if
          CALL PrimToCons(Prim(1:5),U_Loc)
          Flux(:,p,q,SideID) = U_Loc
        END DO ! p
      END DO ! q
    END DO !iSide=1,nBCloc
  CASE DEFAULT !  check for BCtypes in Testcase
    CALL TestcaseLiftingGetBoundaryFlux(iBC,tIn,Flux)
  END SELECT ! BCType
END DO ! iBC
!for BR1 and BR2, lifting is in strong form, flux=U-Uminus...
DO SideID=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    ! convert the conservative variables to primitive/entropy variables based on PP_Lifting_Var
#if (PP_Lifting_Var==1)
    ! conservative variable lifting: do nothing
    F = Flux(:,p,q,SideID)
    F_m = U_Master(:,p,q,SideID)
#elif (PP_Lifting_Var==2)
    ! primitive variable lifting
    ConsToPrim(F, Flux(:,p,q,SideID))
    ConsToPrim(F_m, U_Master(:,p,q,SideID))
#elif (PP_Lifting_Var==3)
    ! entropy variable lifting
    F = ConsToEntropy(Flux(:,p,q,SideID))
    F_m = ConsToEntropy(U_Master(:,p,q,SideID))
#endif /*PP_Lifting_Var*/
    Flux(:,p,q,SideID)=(F-F_m)*SurfElem(p,q,SideID)
  END DO; END DO
END DO ! iSide
END SUBROUTINE Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/


!==================================================================================================================================
!> Finalize arrays used for boundary conditions.
!==================================================================================================================================
SUBROUTINE FinalizeBC()
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,nBCByType,BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(BCData)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC


END MODULE MOD_GetBoundaryFlux
