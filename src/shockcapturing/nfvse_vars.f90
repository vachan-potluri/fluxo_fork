!==================================================================================================================================
! Copyright (c) 2020 - 2021 Andrés Rueda
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
!> summary: Contains global variables for the NFVSE routines
!===================================================================================================================================
module MOD_NFVSE_Vars
#if FLUXO_HYPERSONIC
  use function_parser
#endif
  implicit none
  
  private
  public :: ComputeAlpha, alpha_max, alpha_min, ShockBlendCoef, sharpness, threshold, ModalThreshold
  public :: SubCellMetrics, sWGP, alpha, alpha_Master, alpha_Slave
#if FLUXO_HYPERSONIC
  public :: alpha_vis, alpha_vis_Master, alpha_vis_Slave, wall_blender_limit_parser, viscous_blending_region_parser
#endif
  public :: SubCellMetrics_t, InnerFaceMetrics_t
#if MPI
  public :: MPIRequest_alpha
#endif /*MPI*/
#if NFVSE_CORR
  public :: FFV_m_FDG, alpha_old, PositCorrFactor, PositMaxIter, maximum_alpha, amount_alpha, amount_alpha_steps
#endif /*NFVSE_CORR*/
  public :: sdxR, sdxL, rL, rR, U_ext 
  public :: Compute_FVFluxes, SubFVMethod
  public :: ReconsBoundaries, MPIRequest_Umaster
  public :: RECONS_CENTRAL, RECONS_NONE, RECONS_NEIGHBOR
  public :: SpacePropFactor, SpacePropSweeps, TimeRelFactor
  public :: TanDirs1, TanDirs2
!-----------------------------------------------------------------------------------------------------------------------------------
! New types
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Metric terms on an inner face
  ! -----------------------------
  type :: InnerFaceMetrics_t
    real, allocatable :: nv(:,:,:,:) !< Sub-cell normal vectors            : (1:3,0:N,0:N,-1:N)
    real, allocatable :: t1(:,:,:,:) !< Sub-cell tangent vectors (1)       : (1:3,0:N,0:N,-1:N)
    real, allocatable :: t2(:,:,:,:) !< Sub-cell tangent vectors (2)       : (1:3,0:N,0:N,-1:N)
    real, allocatable :: norm(:,:,:) !< Norm of the sub-cell normal vectors:     (0:N,0:N,-1:N)
!                             | | |_ Inner face index (FV BC idx)
!                             | |___ (FV subcell idx)
!                             |_____ (FV subcell idx)
  end type InnerFaceMetrics_t
  
  ! All the metric terms of an element
  ! ----------------------------------
  type :: SubCellMetrics_t
    type(InnerFaceMetrics_t) :: xi
    type(InnerFaceMetrics_t) :: eta
    type(InnerFaceMetrics_t) :: zeta
    contains
      procedure :: construct => SubCellMetrics_construct
      procedure :: TestMetricIdentities => SubCellMetrics_TestMetricIdentities
  end type SubCellMetrics_t
  
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables
!-----------------------------------------------------------------------------------------------------------------------------------
  
! General
! -------
  integer                             :: ComputeAlpha     !< Method to compute alpha
  real, target          , allocatable :: alpha(:)         !< Element-wise blending function (modified every time Ut is computed)
  real                  , allocatable :: alpha_Master(:)  !< Blending function on master sides
  real                  , allocatable :: alpha_Slave(:)   !< Blending function on slave sides
#if FLUXO_HYPERSONIC
  real, target, allocatable           :: alpha_vis(:)     !< Blending function for viscous residual (it is only scaled though)
  ! function parsers for wall blender limit and viscous blending region
  real, allocatable                   :: alpha_vis_Master(:), alpha_vis_Slave(:)
  type (fparser)                      :: wall_blender_limit_parser, viscous_blending_region_parser
#endif
  real                                :: threshold        !< Threshold for the shock indicator
  real, parameter                     :: sharpness = log((1.0-1.e-4)/1.e-4) !< Heuristically obtained sharpness for the shock indicator
  real                                :: alpha_max        !< Maximum blending factor
  real                                :: alpha_min        !< Minimum blending factor
  real                                :: ShockBlendCoef
  integer                             :: ModalThreshold
  type(SubCellMetrics_t), allocatable :: SubCellMetrics(:)      !< Metric terms for the native sub-cell finite volumes
  real                  , allocatable :: sWGP(:)                !< Inverse of the Gauss quadrature weights
#if MPI
  integer               , allocatable :: MPIRequest_alpha(:,:)  !< MPI request for the transfer of the blending coefficient
                                                                !  (nNbProcs,4)... 1: send slave, 2: send master, 3: receive slave, 4, receive master
#endif /*MPI*/
  
! Definition of the tangent directions for the "inner faces" (note that they differ from the definition of TangDirs in MOD_Mesh_Vars)
! ----------------------------------------------------------
  INTEGER,PARAMETER :: TanDirs1(6)  = (/ 1 , 1 , 2 , 1 , 2 , 1 /) !< first tangential vector direction for local "inner faces"
  INTEGER,PARAMETER :: TanDirs2(6)  = (/ 2 , 3 , 3 , 3 , 3 , 2 /) !< first tangential vector direction for local "inner faces"
  
! For the positivity limiter
! --------------------------
#if NFVSE_CORR
  real, allocatable :: FFV_m_FDG(:,:,:,:,:)
  real, allocatable :: alpha_old(:)                         !< Element-wise blending function (before correction)
  real              :: PositCorrFactor  ! Limiting factor for NFVSE correction
  integer           :: PositMaxIter
  real              :: maximum_alpha     ! Maximum alpha for the analyze routines
  real              :: amount_alpha
  integer           :: amount_alpha_steps
#endif /*NFVSE_CORR*/

! For the FV method
! -----------------
  integer                                    :: SubFVMethod
  procedure(i_sub_Compute_FVFluxes), pointer :: Compute_FVFluxes => null()

! Space propagation and time relaxation
  real              :: SpacePropFactor  ! Space propagation factor 
  real              :: TimeRelFactor    ! Time relaxation factor 
  integer           :: SpacePropSweeps  ! Number of space propagation sweeps

! For the reconstruction procedure
! --------------------------------
  real                      , allocatable :: sdxR(:), sdxL(:)       !< Inverse of subgrid sizes for reconstruction
  real                      , allocatable :: rR(:), rL(:)
  real                      , allocatable :: U_ext(:,:,:,:,:)       !< External solution for reconstruction on boundaries (PP_nVar,0:N,0:N,locside,iElem)
  integer                   , allocatable :: MPIRequest_Umaster(:,:)
  
  integer                                 :: ReconsBoundaries = 1
  
  integer, parameter :: RECONS_NONE     = 1
  integer, parameter :: RECONS_CENTRAL  = 2
  integer, parameter :: RECONS_NEIGHBOR = 3
!-----------------------------------------------------------------------------------------------------------------------------------
! Interfaces
!-----------------------------------------------------------------------------------------------------------------------------------
  abstract interface
    subroutine i_sub_Compute_FVFluxes (U, F , G , H , &
#if NONCONS
                                          FR, GR, HR, &
#endif /*NONCONS*/
                                          sCM, iElem )
      use MOD_PreProc
      import SubCellMetrics_t
      real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N), intent(in)    :: U   !< The element solution
      real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: F   !< Left flux in xi
      real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: G   !< Left flux in eta
      real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: H   !< Left flux in zeta
#if NONCONS
      real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: FR  !< Right flux in xi
      real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: GR  !< Right flux in eta
      real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: HR  !< Right flux in zeta
#endif /*NONCONS*/
      type(SubCellMetrics_t)                         , intent(in)    :: sCM       !< Sub-cell metric terms
      integer                                        , intent(in)    :: iElem
    end subroutine i_sub_Compute_FVFluxes
  end interface
!===================================================================================================================================
  contains
    elemental subroutine SubCellMetrics_construct(this,N)
      implicit none
      class(SubCellMetrics_t), intent(inout) :: this
      integer                , intent(in)    :: N
      
      allocate ( this % xi   % nv (1:3,0:N,0:N,-1:N) )
      allocate ( this % xi   % t1 (1:3,0:N,0:N,-1:N) )
      allocate ( this % xi   % t2 (1:3,0:N,0:N,-1:N) )
      allocate ( this % xi   % norm   (0:N,0:N,-1:N) )
      
      allocate ( this % eta  % nv (1:3,0:N,0:N,-1:N) )
      allocate ( this % eta  % t1 (1:3,0:N,0:N,-1:N) )
      allocate ( this % eta  % t2 (1:3,0:N,0:N,-1:N) )
      allocate ( this % eta   % norm   (0:N,0:N,-1:N) )
      
      allocate ( this % zeta % nv (1:3,0:N,0:N,-1:N) )
      allocate ( this % zeta % t1 (1:3,0:N,0:N,-1:N) )
      allocate ( this % zeta % t2 (1:3,0:N,0:N,-1:N) )
      allocate ( this % zeta % norm   (0:N,0:N,-1:N) )
      
    end subroutine SubCellMetrics_construct
    
    subroutine SubCellMetrics_TestMetricIdentities(this,iSC,jSC,kSC,div)
      use MOD_Interpolation_Vars , only: wGP
      implicit none
      class(SubCellMetrics_t), intent(inout) :: this
      integer                , intent(in)    :: iSC, jSC, kSC !< Subcell index: 1...N-1
      real, dimension(3)     , intent(out)   :: div
      !----------------------------------------------
      real, dimension(3) :: a, b, c, d, e, f
      real               :: anorm, bnorm, cnorm, dnorm, enorm, fnorm
      !----------------------------------------------
        
      a     = this % xi   % nv(:,jSC,kSC,iSC)
      anorm = this % xi   % norm(jSC,kSC,iSC)
      b     = this % xi   % nv(:,jSC,kSC,iSC-1)
      bnorm = this % xi   % norm(jSC,kSC,iSC-1)
      c     = this % eta  % nv(:,iSC,kSC,jSC)
      cnorm = this % eta  % norm(iSC,kSC,jSC)
      d     = this % eta  % nv(:,iSC,kSC,jSC-1)
      dnorm = this % eta  % norm(iSC,kSC,jSC-1)
      e     = this % zeta % nv(:,iSC,jSC,kSC)
      enorm = this % zeta % norm(iSC,jSC,kSC)
      f     = this % zeta % nv(:,iSC,jSC,kSC-1)
      fnorm = this % zeta % norm(iSC,jSC,kSC-1)
      
      div = (a * anorm - b * bnorm) / wGP(iSC) + (c * cnorm - d * dnorm) / wGP(jSC) + (e * enorm - f * fnorm) / wGP(kSC)
    end subroutine SubCellMetrics_TestMetricIdentities
end module MOD_NFVSE_Vars
