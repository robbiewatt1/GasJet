!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Mach 3 wind tunnel
!!  problem.
!!
!!  References:  Emery, A. E., 1968, JCP, 2, 306
!!               Woodward, P. and Colella, P., 1984, JCP, 54, 115
!!
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_pAmbient, sim_tAmbient, sim_gamma, &
     &  sim_smallP, sim_smallX, eos_singlespeciesa, sim_nozIn, sim_nozExit, &
     sim_nozHeight, sim_pBacking
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkPtr, Grid_releaseBlkPtr


  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer,intent(IN) :: blockID
  real,pointer :: solnData(:,:,:,:)
  integer :: axis(MDIM)  
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: xMax, yMax, zMax
  real :: wall_thickness
  integer:: i, j, k

  real :: rho_zone, velx_zone, vely_zone, velz_zone, pres_zone, &
       temp_zone, ener_zone, ekin_zone, eint_zone, rho_chamber

  integer :: Wall_spec = 1, Cham_spec = 2, Vac_spec = 3 
  integer :: species


!===============================================================================


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))
  allocate(zcent(blkLimitsGC(HIGH, KAXIS)))
  call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., &
       zcent, blkLimitsGC(HIGH, KAXIS))


  xMax = xcent(size(xcent))
  yMax = xcent(size(ycent))
  zMax = xcent(size(zcent))
  wall_thickness = (xMax - sim_nozIn) / 2.0


! In this problem the initial conditions are spatially uniform.
  
  pres_zone = sim_pAmbient
  temp_zone = sim_tAmbient
  rho_zone = (pres_zone / temp_zone) / (8.3107e7 / eos_singlespeciesA)
  rho_chamber = (sim_pBacking / temp_zone) / (8.3107e7 / eos_singlespeciesA) 

  velx_zone = 0.0
  vely_zone = 0.0
  velz_zone = 0.0


  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
  ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)

  eint_zone = pres_zone / (sim_gamma-1.)
  eint_zone = eint_zone / rho_zone
  ener_zone = eint_zone + ekin_zone
  ener_zone = max(ener_zone, sim_smallP)

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

          if (ycent(j) <= sim_nozHeight ) then
            if (xcent(i) <= wall_thickness ) then !.and. xcent(i) <= wall_thickness + sim_nozIn) then
              species = Cham_spec
            else 
              species = Wall_spec
            endif
          else
            species = Vac_spec
          endif

          if (species == Cham_spec) then
            call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, -1.0)
            call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho_chamber)
          elseif (species == Wall_spec) then
            call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho_zone*2000)
            call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, 1.0)
          elseif (species == Vac_spec) then
            call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, -1.0)
            call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho_zone)
          endif

           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, ener_zone)    
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velx_zone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vely_zone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velz_zone)

        enddo
      enddo
  enddo




 ! call Grid_releaseBlkPtr(blockID, solnData, CENTER)
  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

  return
end subroutine Simulation_initBlock



