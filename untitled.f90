  call Grid_getBlkPtr(blockID, solnData, CENTER)
#if NSPECIES > 0
!  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallX
!  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallX
#endif

  ! store the variables in the block's unk data
  solnData(DENS_VAR,:,:,:) = rho_zone
  solnData(PRES_VAR,:,:,:) = pres_zone
  solnData(ENER_VAR,:,:,:) = ener_zone
#ifdef EINT_VAR
  solnData(EINT_VAR,:,:,:) = eint_zone
#endif
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(GAME_VAR,:,:,:) = sim_gamma

  solnData(VELX_VAR,:,:,:) = velx_zone
  solnData(VELY_VAR,:,:,:) = vely_zone
  solnData(VELZ_VAR,:,:,:) = velz_zone

#ifdef BDRY_VAR
  if (blockID == 10) then
    solnData(BDRY_VAR,:,:,:) = 1.0
    solnData(DENS_VAR,:,:,:) = 10
  else
    solnData(BDRY_VAR,:,:,:) = -1.0
  endif
#endif