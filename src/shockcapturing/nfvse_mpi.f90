!==================================================================================================================================
! Copyright (c) 2020 - 2020 Andrés Rueda
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
!===================================================================================================================================
!> Module containing utilities for the NFVSE shock-capturing method
!> Attention: this module also contains the MPI routines for the communication of the blending coefficient in the NFVSE method:
!> -> Done here and not in DGTimeDerivative_weakForm for modularity!
!> -> Here the MPI approach is a bit different than in the main code: Both master and slave sides send and receive, then the
!>    blending coefficient correction is computed on both sides of the MPI interface. This is needed since we want to correct the DG
!>    volume integral BEFORE adding the diffusive contribution.
!===================================================================================================================================
module MOD_NFVSE_MPI
  implicit none
  
  private
  public ProlongBlendingCoeffToFaces, PropagateBlendingCoeff
#if FLUXO_HYPERSONIC
  public UpdateVisBlendingCoefficient
#endif
contains
!===================================================================================================================================
!> Prolong the blending coefficient to the faces
!===================================================================================================================================
  subroutine ProlongBlendingCoeffToFaces()
    use MOD_NFVSE_Vars, only: alpha, alpha_Master, alpha_Slave
    use MOD_Mesh_Vars , only: SideToElem, firstMortarInnerSide, nSides
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    integer :: sideID, ElemID, nbElemID
    !-------------------------------------------------------------------------------------------------------------------------------
    
!   First prolong alpha to all non-BC sides!
!   ----------------------------------------
    do sideID=firstMortarInnerSide, nSides
      ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side
      !master sides(ElemID,locSide and flip =-1 if not existing)
      if(ElemID.NE.-1) then ! element belonging to master side is on this processor
        alpha_Master(SideID) = alpha(ElemID) 
      end if
      
      nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
      !slave side (nbElemID,nblocSide and flip =-1 if not existing)
      if(nbElemID.NE.-1) then! element belonging to slave side is on this processor
        alpha_Slave (SideID) = alpha(nbElemID) 
      end if
    end do
    
!   Fill the small INNER non-conforming sides with the big sides' alpha
!   ATTENTION: This has to be done before calling Start_BlendCoeff_MPICommunication!!
!   ---------------------------------------------------------------------------------
    call Alpha_Mortar(alpha_Master,alpha_Slave,doMPISides=.FALSE.)
    
!   Send and receive alpha from master and slave sides
!   --------------------------------------------------
#if MPI
    call Start_BlendCoeff_MPICommunication()
#endif /*MPI*/
    
  end subroutine ProlongBlendingCoeffToFaces
!===================================================================================================================================
!> Propagate the blending coefficient to the neighbor elements
!> -> First the routine makes sure that the MPI communication is done
!> -> Then, the blending coefficient is set to alpha = max (alpha, alpha_neighbor)
!===================================================================================================================================
  subroutine PropagateBlendingCoeff()
    use MOD_NFVSE_Vars         , only: alpha, alpha_Master, alpha_Slave
    use MOD_Mesh_Vars          , only: firstSlaveSide, LastSlaveSide, SideToElem
    use MOD_Mesh_Vars          , only: firstMortarInnerSide, lastMortarInnerSide, firstMortarMPISide,lastMortarMPISide
    use MOD_NFVSE_Vars         , only: SpacePropFactor
    use MOD_Mesh_Vars          , only: MortarType,MortarInfo
#if MPI
    use MOD_MPI_Vars           , only: nNbProcs
    use MOD_MPI                , only: FinishExchangeMPIData
    use MOD_NFVSE_Vars         , only: MPIRequest_alpha
#endif /*MPI*/
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    integer :: sideID, ElemID, nbElemID
    real    :: minAlpha (firstSlaveSide:LastSlaveSide)
    integer :: MortarSideID, nMortars, locSide, iMortar
    !-------------------------------------------------------------------------------------------------------------------------------
    
!   Finish the MPI exchange
!   -----------------------
#if MPI
    call FinishExchangeMPIData(4*nNbProcs,MPIRequest_alpha)
#endif /*MPI*/
    
!   Compute the minimum alpha (for each pair master/slave) and enforce it on non-mortar sides
!   -----------------------------------------------------------------------------------------
    do sideID=firstSlaveSide, lastSlaveSide
      
      minAlpha(sideID) = max( alpha_Master(sideID), alpha_Slave(sideID) )
      minAlpha(sideID) = SpacePropFactor * minAlpha(sideID)
      
      ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side

      !master sides(ElemID,locSide and flip =-1 if not existing)
      if(ElemID.NE.-1) then ! element belonging to master side is on this processor
        if (alpha(ElemID) < minAlpha(sideID)) alpha(ElemID) = minAlpha(sideID)
      end if
      
      nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
      !slave side (nbElemID,nblocSide and flip =-1 if not existing)
      if(nbElemID.NE.-1) then! element belonging to slave side is on this processor
        if (alpha(nbElemID) < minAlpha(sideID)) alpha(nbElemID) = minAlpha(sideID)
      end if
    end do

!   Enforce minimum alpha on mortar sides
!   -------------------------------------
    
    ! Inner mortar sides
    DO MortarSideID=firstMortarInnerSide, lastMortarInnerSide

      !Save the small sides into master/slave arrays
      IF(MortarType(1,MortarSideID).EQ.1)THEN
        nMortars=4
      ELSE
        nMortars=2
      END IF !MortarType
      locSide=MortarType(2,MortarSideID)
      ElemID =SideToElem(S2E_ELEM_ID,MortarSideID) !element belonging to master side
      DO iMortar=1,nMortars
        SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
        if (alpha(ElemID) < minAlpha(sideID)) alpha(ElemID) = minAlpha(sideID)
      END DO !iMortar
    END DO !MortarSideID
    
    ! MPI mortar sides
    DO MortarSideID=firstMortarMPISide,lastMortarMPISide

      !Save the small sides into master/slave arrays
      IF(MortarType(1,MortarSideID).EQ.1)THEN
        nMortars=4
      ELSE
        nMortars=2
      END IF !MortarType
      locSide=MortarType(2,MortarSideID)
      ElemID =SideToElem(S2E_ELEM_ID,MortarSideID) !element belonging to master side
      DO iMortar=1,nMortars
        SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
        if (alpha(ElemID) < minAlpha(sideID)) alpha(ElemID) = minAlpha(sideID)
      END DO !iMortar
    END DO !MortarSideID

#if FLUXO_HYPERSONIC
    ! update viscous blending coefficient
    ! the viscous blending region is generally not a function of time, so this function can
    ! generally be called without any argument
    call UpdateVisBlendingCoefficient()
#endif
  end subroutine
#if MPI
!===================================================================================================================================
!> Start the MPI communication of the blending coefficient (send and receive)
!===================================================================================================================================
  subroutine Start_BlendCoeff_MPICommunication()
    use MOD_NFVSE_Vars         , only: alpha_Master, alpha_Slave
    use MOD_MPI                , only: StartReceiveMPIData,StartSendMPIData
    use MOD_Mesh_Vars          , only: firstSlaveSide, LastSlaveSide,firstMortarInnerSide,nSides
    use MOD_NFVSE_Vars         , only: MPIRequest_alpha
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    
!   Start off with the receive command
!   **********************************
    ! receive the slave
    call StartReceiveMPIData(alpha_Slave , 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,1), SendID=2) ! Receive MINE (sendID=2)
    
    ! receive the master
    call StartReceiveMPIData(alpha_Master, 1, firstMortarInnerSide,nSides, &
                             MPIRequest_alpha(:,2), SendID=1) ! Receive YOUR  (sendID=1) 
    
    ! transfer the mortar to the right MPI faces 
    CALL Alpha_Mortar(alpha_Master,alpha_Slave,doMPISides=.TRUE.)

    ! Send the slave
    call StartSendMPIData   (alpha_Slave , 1, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_alpha(:,3), SendID=2) ! SEND YOUR (sendID=2) 
    
    ! Send the master
    call StartSendMPIData   (alpha_Master, 1, firstMortarInnerSide, nSides, &
                             MPIRequest_alpha(:,4),SendID=1) 
    
  end subroutine Start_BlendCoeff_MPICommunication
#endif /*MPI*/
!==================================================================================================================================
!> Fills small non-conforming sides with data from the corresponding large side
!==================================================================================================================================
  subroutine Alpha_Mortar(alpha_Master,alpha_Slave,doMPISides)
    ! MODULES
    use MOD_Preproc
    use MOD_Mesh_Vars,   only: MortarType,MortarInfo
    use MOD_Mesh_Vars,   only: firstMortarInnerSide,lastMortarInnerSide
    use MOD_Mesh_Vars,   only: firstMortarMPISide,lastMortarMPISide
    use MOD_Mesh_Vars,   only: firstSlaveSide,lastSlaveSide
    use MOD_Mesh_Vars,   only: nSides 
    implicit none
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    real,INTENT(INOUT) :: alpha_Master(firstMortarInnerSide:nSides)  !< (INOUT) 
    real,INTENT(INOUT) :: alpha_Slave (firstSlaveSide:lastSlaveSide) !< (INOUT) 
    LOGICAL,INTENT(IN) :: doMPISides                                 !< flag whether MPI sides are processed
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER      :: iMortar,nMortars
    INTEGER      :: firstMortarSideID,lastMortarSideID
    INTEGER      :: MortarSideID,SideID,locSide,flip
    !==================================================================================================================================
    IF(doMPISides)THEN
      firstMortarSideID = firstMortarMPISide
      lastMortarSideID =  lastMortarMPISide
    ELSE
      firstMortarSideID = firstMortarInnerSide
      lastMortarSideID =  lastMortarInnerSide
    END IF !doMPISides

    DO MortarSideID=firstMortarSideID,lastMortarSideID

      !Save the small sides into master/slave arrays
      IF(MortarType(1,MortarSideID).EQ.1)THEN
        nMortars=4
      ELSE
        nMortars=2
      END IF !MortarType
      locSide=MortarType(2,MortarSideID)
      DO iMortar=1,nMortars
        SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
        flip  = MortarInfo(MI_FLIP,iMortar,locSide)
        SELECT CASE(flip)
          CASE(0) ! master side
            alpha_Master(SideID)=Alpha_Master(MortarSideID)
          CASE(1:4) ! slave side
            alpha_Slave (SideID)=Alpha_Master(MortarSideID)
        END SELECT !flip(iMortar)
      END DO !iMortar
    END DO !MortarSideID
  end subroutine Alpha_Mortar

! Updates the viscous blending coefficient based on current time and current values of alpha
#if FLUXO_HYPERSONIC
  subroutine UpdateVisBlendingCoefficient(t)
    use MOD_Mesh_Vars          , only: nElems, Elem_centers
    use MOD_NFVSE_Vars         , only: alpha, alpha_max, alpha_vis, viscous_blending_region_parser
    use iso_fortran_env        , only: output_unit

    implicit none
    real, intent(in), optional :: t
    real, dimension(4) :: parser_vals
    real :: parser_result
    integer :: eID

    do eID=1,nElems
      if(present(t)) then
        parser_vals = [Elem_centers(1,eID), Elem_centers(2,eID), Elem_centers(3,eID), t]
      else
        parser_vals = [Elem_centers(1,eID), Elem_centers(2,eID), Elem_centers(3,eID), 0.0]
      end if
      CALL viscous_blending_region_parser%evaluate(parser_vals, parser_result)
      if(viscous_blending_region_parser%error()) then
        call viscous_blending_region_parser%print_errors(output_unit)
      else
        if(parser_result > 0) then
          alpha_vis(eID) = alpha(eID)/alpha_max
        else
          alpha_vis(eID) = 0.0
        end if
      end if
    end do
  end subroutine UpdateVisBlendingCoefficient
#endif /*FLUXO_HYPERSONIC*/
end module MOD_NFVSE_MPI
