!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

program GenNASATileFile
  use ESMF
  !  use NUOPC
  use NASACatchFile
  
  implicit none
  integer :: localrc
  type(ESMF_Field) :: lndField
  type(ESMF_Mesh) :: lndMesh
  type(ESMF_Grid) :: atmGrid
  type(ESMF_Grid) :: ocnGrid
  type(ESMF_XGrid) :: xgrid
  type(ESMF_Mesh) :: xgridMesh

  
  ! Initialize ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
    defaultCalkind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
#if 0
  ! Create Field from file
  lndField=NCF_CreateField("catch_rast_small.nc", rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif


  ! Create C12 Atm Grid
  atmGrid= ESMF_GridCreateCubedSphere(tilesize=12, name="ATM-Grid", &
       staggerLocList=[ESMF_STAGGERLOC_CENTER,ESMF_STAGGERLOC_CORNER], &
       coordSys=ESMF_COORDSYS_SPH_DEG, rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  
  ! Create lnd catchment mesh from file
  lndMesh=NCF_CreateCatchmentMesh("catch_rast_small.nc", rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

#if 0  
  ! Debug output of catchment mesh
  call ESMF_MeshWrite(lndMesh,"lndMesh", rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif  

  ! Create ocean grid
  ocnGrid = ESMF_GridCreate1PeriDimUfrm(maxIndex=(/72,36/), &
       minCornerCoord=(/-180._ESMF_KIND_R8, -90._ESMF_KIND_R8/), &
       maxCornerCoord=(/180._ESMF_KIND_R8, 90._ESMF_KIND_R8/), &
       staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), name="LND-Grid", &
       rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Create XGrid  
  xgrid=ESMF_XGridCreate(&
       sideAGrid=(/atmGrid/), &
       sideBGrid=(/ocnGrid/), &
       sideBGridPriority=(/1/), &
       sideBMesh=(/lndMesh/), &
       sideBMeshPriority=(/2/), &
       storeOverlay=.true., &
       rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Get Mesh (for debugging purposes)
  call ESMF_XGridGet(xgrid, mesh=xgridMesh, rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Debug output of xgrid mesh
  call ESMF_MeshWrite(xgridMesh,"xgridMesh", rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)




  
  ! Get rid of XGrid
  call ESMF_XGridDestroy(xgrid, rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Get rid of lnd Mesh
  call ESMF_MeshDestroy(lndMesh, rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  ! Get rid of ocn Grid
  call ESMF_GridDestroy(ocnGrid, rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Get rid of atm Grid
  call ESMF_GridDestroy(atmGrid, rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  

  
  
  
  ! Finalize ESMF
  call ESMF_Finalize()
  
end program  
