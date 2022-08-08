!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

  ! This is a module that provides some basic functionality for doing some 
  ! things with ESMF and the NASA Catchment file

module NASACatchFile
  use ESMF
  use NetCDF
  
  implicit none

#define NCF_SRCLINE __FILE__,__LINE__
#define NCF_CONTEXT file=__FILE__,line=__LINE__
#define NCF_PASSTHRU msg="Internal subroutine call returned Error"

  ! Public interfaces available from this module
  public NCF_CreateCatchmentMesh
  public NCF_CreateField
  public NCF_CreateGrid

  ! Make things private, except that which is explicitly make public
  private
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

    ! Create an ESMF Mesh to represent the catchments in a file
    function NCF_CreateCatchmentMesh(fileName, rc)
#undef  NCF_METHOD
#define NCF_METHOD "NCF_CreateCatchmentMesh()"

      ! Return value
      type(ESMF_Mesh) :: NCF_CreateCatchmentMesh
      !
      ! !ARGUMENTS:
      character (len=*),         intent(in),  optional :: fileName
      integer,                   intent(out), optional :: rc

      ! Local variables
      integer :: localrc
      type(ESMF_Grid) :: grid
      type(ESMF_Field) :: field
      type(ESMF_Array) :: array

      ! Create Field from file
      field=NCF_CreateField(fileName, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return


      ! Get Grid and Array
      call ESMF_FieldGet(field, &
           grid=grid, array=array, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Create Mesh from raster info
      NCF_CreateCatchmentMesh=ESMF_MeshCreate(&
           rasterGrid=grid, rasterArray=array, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return
      
      
      ! Get rid of Field
      call ESMF_FieldDestroy(field, rc=localrc) 
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end function NCF_CreateCatchmentMesh
    

    ! Creates an ESMF Field to hold the raster data,
    ! index space, and coordinates of a NASA catchment file
    function NCF_CreateField(fileName, rc)
#undef  NCF_METHOD
#define NCF_METHOD "NCF_CreateField()"

      ! Return value
      type(ESMF_Field) :: NCF_CreateField
      !
      ! !ARGUMENTS:
      character (len=*),         intent(in),  optional :: fileName
      integer,                   intent(out), optional :: rc

      ! Local variables
      integer :: localrc
      type(ESMF_Grid) :: grid

      ! Debug
      integer :: lbnd(2),ubnd(2)
      integer(ESMF_KIND_I4), pointer :: fieldPtr(:,:)

      
      ! Create Grid from file
      grid=NCF_CreateGrid(fileName, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Create Field on Grid
      NCF_CreateField=ESMF_FieldCreate(grid, &
           typekind=ESMF_TYPEKIND_I4, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Read raster info into Field
      call ESMF_FieldRead(NCF_CreateField, &
           fileName, variableName="CatchIndex", &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

#if 0
      ! Debug output
      call ESMF_FieldGet(NCF_CreateField, localDE=0, &
           farrayPtr=fieldPtr, rc=localrc) 
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return         

      lbnd(1)=lbound(fieldPtr,1)
      lbnd(2)=lbound(fieldPtr,2)
      write(*,*) "Field(",lbnd(1),lbnd(2),")=",fieldPtr(lbnd(1),lbnd(2))

      ubnd(1)=ubound(fieldPtr,1)
      ubnd(2)=ubound(fieldPtr,2)
      write(*,*) "Field(",ubnd(1),ubnd(2),")=",fieldPtr(ubnd(1),ubnd(2))
#endif
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end function NCF_CreateField

    
    ! Creates a Grid to represent the index space and coordinates of
    ! a NASA catchment file
    function NCF_CreateGrid(fileName, rc)
#undef  NCF_METHOD
#define NCF_METHOD "NCF_CreateGrid()"

      ! Return value
      type(ESMF_Grid) :: NCF_CreateGrid
      !
      ! !ARGUMENTS:
      character (len=*),         intent(in),  optional :: fileName
      integer,                   intent(out), optional :: rc

      ! Local variables
      integer :: ncStatus, localrc
      integer :: ncid, dimid
      integer :: lonSize, latSize
      
      write(*,*) "In cGFNCF filename=",fileName

      ! Open file using netcdf
      ncStatus = nf90_open (path=trim(filename), mode=nf90_nowrite, ncid=ncid)
      if (NetCDFCheckError (ncStatus, &
           NCF_METHOD,  &
           NCF_SRCLINE,  &
           "Opening file "//trim(filename), &
           rc)) return
      
      
      ! Get longitude (x) size
      ncStatus = nf90_inq_dimid(ncid, "N_lon", dimid)
      if (NetCDFCheckError (ncStatus, &
           NCF_METHOD, &
           NCF_SRCLINE,&
           "Opening N_lon dimension in "//trim(fileName),&
           rc)) return
      
      ncStatus = nf90_inquire_dimension(ncid, dimId, len=lonSize)
      if (NetCDFCheckError (ncStatus, &
           NCF_METHOD, &
           NCF_SRCLINE,&
           "Reading N_lon dimension from "//trim(fileName),&
           rc)) return

      ! Get latitude (y) size
      ncStatus = nf90_inq_dimid(ncid, "N_lat", dimid)
      if (NetCDFCheckError (ncStatus, &
           NCF_METHOD, &
           NCF_SRCLINE,&
           "Opening N_lon dimension in "//trim(fileName),&
           rc)) return
      
      ncStatus = nf90_inquire_dimension(ncid, dimId, len=latSize)
      if (NetCDFCheckError (ncStatus, &
           NCF_METHOD, &
           NCF_SRCLINE,&
           "Reading N_lon dimension from "//trim(fileName),&
           rc)) return

      ! Close file
      ncStatus = nf90_close(ncid)
      if (NetCDFCheckError (ncStatus, &
           NCF_METHOD,  &
           NCF_SRCLINE, &
           "Closing file "//trim(filename), &
           rc)) return
      
      
      ! Debug output
      write(*,*) "Grid size=",lonSize," x ",latSize

      ! Create Grid
      NCF_CreateGrid=ESMF_GridCreate1PeriDim( &
           maxIndex=(/lonSize,latSize/), &
           coordSys=ESMF_COORDSYS_SPH_DEG, &
           coordDep1=(/1/),  coordDep2=(/2/), &
           indexflag=ESMF_INDEX_GLOBAL,  rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Add Center coords
      call ESMF_GridAddCoord(NCF_CreateGrid, &
           staggerloc=ESMF_STAGGERLOC_CENTER, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Fill longitude center coordinates from file
      call readCenterCoordsFromNASACatchFile(NCF_CreateGrid, &
           filename, coordDim=1, varname="longitude", rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Fill latitude center coordinates from file
      call readCenterCoordsFromNASACatchFile(NCF_CreateGrid, &
           filename, coordDim=2, varname="latitude", rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      
      
      ! Add Corner coords
      call ESMF_GridAddCoord(NCF_CreateGrid, &
           staggerloc=ESMF_STAGGERLOC_CORNER, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Calculate corner coordinates from center
      call calcCornerLongCoordsFromCenter(NCF_CreateGrid, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      call calcCornerLatCoordsFromCenter(NCF_CreateGrid, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end function NCF_CreateGrid


    ! Read center coordinates and put into Grid 
    ! This only supports reading 1D coordinates
    ! (The coordinates in NASA files are 1D)
    subroutine readCenterCoordsFromNASACatchFile(grid, filename, &
         coordDim, varname, rc)
      
      ! Arguments
      type(ESMF_Grid) :: grid
      character(len=*), intent(in)  :: filename
      integer :: coordDim
      character(len=*), intent(in)  :: varname
      integer, intent(out),optional :: rc      

      ! Local variables
      integer :: localrc 
      type(ESMF_Array) :: coordArray, tmpCoordArray
      type(ESMF_Distgrid) :: gridArrayDistGrid
      type(ESMF_Distgrid) :: factorDistGrid
      integer :: localDECount,rank
      type(ESMF_LocalArray), allocatable :: localArrayList(:)
      real(ESMF_KIND_R8), pointer :: farrayPtr(:)
 
     
      ! Get lon. coordinate Array from Grid
      call ESMF_GridGetCoord(grid, coordDim=coordDim, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           array=coordArray, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      


      ! Get DistGrid and localDECount from Array
      call ESMF_ArrayGet(coordArray,  &
           rank=rank, &
           localDECount=localDECount, &
           distgrid=gridArrayDistGrid, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      
      ! Only support 1D coordinates right now
      if (rank .ne. 1) then
         call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
              msg="Only Grids with 1D coordinate arrays are currently supported.", &
              NCF_CONTEXT, rcToReturn=rc)
         return
      endif

      ! Allocate space for localArrayList
      allocate(localArrayList(localDECount))

      ! Get localArrayList
      call ESMF_ArrayGet(coordArray,  &
           localArrayList=localArrayList, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

    
      ! Get index Info from DistGrid
      call createFactorDistGrid(gridArrayDistgrid, coordDim, factorDistgrid, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Create factored Array
      tmpCoordArray=ESMF_ArrayCreate(factorDistgrid, &
           localArrayList=localArrayList, &
           dataCopyFlag=ESMF_DATACOPY_REFERENCE, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return
      
      ! Read lon. coordinates into temporary Array
      call ESMF_ArrayRead(tmpCoordArray, &
           fileName=filename, variableName=varname, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      


      ! DEBUG OUTPUT
      !! call ESMF_ArrayGet(coordArray, localDE=0, farrayPtr=farrayPtr, rc=localrc) 
      !! if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
      !!     NCF_CONTEXT, rcToReturn=rc)) return
      !! 
      !! write(*,*) trim(varname)," = ",farrayPtr
      
      ! Get rid of Array
      call ESMF_ArrayDestroy(tmpCoordArray, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Get rid of Distgrid
      call ESMF_DistGridDestroy(factorDistgrid, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine readCenterCoordsFromNASACatchFile

   
    
    ! Create a distgrid that's just 1D along the factor dim
    subroutine createFactorDistGrid(distgrid, factorDim, factorDistgrid, rc)
      ! Arguments
      type(ESMF_DistGrid) :: distgrid
      integer :: factorDim
      type(ESMF_DistGrid) :: factorDistgrid
      integer, intent(out),optional :: rc      

      ! Local variables
      type(ESMF_Index_Flag) :: indexflag
      type(ESMF_DELayout) :: DELayout
      integer :: localrc 
      integer :: minIndexPTile(2,1), maxIndexPTile(2,1)
      integer :: tileCount,dimCount,deCount
      integer, allocatable :: deBlockList(:,:,:)
      integer, allocatable :: minIndexPDE(:,:)
      integer, allocatable :: maxIndexPDE(:,:)
      integer :: d
      
      ! Get Info from Distgrid
      call ESMF_DistGridGet(distgrid, &
           tileCount=tileCount, &
           dimCount=dimCount, &
           deCount=deCount, &
           indexFlag=indexflag, &
           delayout=delayout, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      
      
      ! Only support 2D 1 tile grids right now
      if (tileCount .ne. 1) then
         call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
              msg="Only 1 tile Grids currently supported.", &
              NCF_CONTEXT, rcToReturn=rc)
         return
      endif

      if (dimCount .ne. 2) then
         call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
              msg="Only 2D Grids currently supported.", &
              NCF_CONTEXT, rcToReturn=rc)
         return
      endif

      ! Allocate min and max de lists
      allocate(minIndexPDe(2,deCount))
      allocate(maxIndexPDe(2,deCount))
      
      
      ! Get Info from Distgrid
      call ESMF_DistGridGet(distgrid, &
           minIndexPTile=minIndexPTile, &
           maxIndexPTile=maxIndexPTile, &
           minIndexPDe=minIndexPDE, &
           maxIndexPDe=maxIndexPDE, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      

      ! Create deBlockList from min and max DE lists
      allocate(deBlockList(1,2,deCount))
      do d=1,deCount
         deBlockList(1,1,d)=minIndexPDE(factorDim,d)
         deBlockList(1,2,d)=maxIndexPDE(factorDim,d)
      enddo
      
      ! Create new distgrid
      factorDistgrid=ESMF_DistGridCreate(minIndex=minIndexPTile(factorDim:factorDim,1), &
           maxIndex=maxIndexPTile(factorDim:factorDim,1), &
           deBlockList=deBlockList, &
           delayout=delayout, &
           indexFlag=indexflag, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      

      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine createFactorDistGrid


    ! Calculate corner longitude coordinates from center longitude coordinates
    ! Only supports 1D R8 coordinate arrays
    subroutine calcCornerLongCoordsFromCenter(grid, rc)
      
      ! Arguments
      type(ESMF_Grid) :: grid
      integer, intent(out),optional :: rc      

      ! Local variables
      integer :: localrc 
      type(ESMF_Array) :: centerCoordArray
      type(ESMF_Array) :: centerCoordArrayWHalo
      integer :: lDE,localDECount
      real(ESMF_KIND_R8), pointer :: centerwHaloPtr(:)
      type(ESMF_Array) :: cornerCoordArray
      real(ESMF_KIND_R8), pointer :: cornerPtr(:)
      integer :: i,lbnd,ubnd
      
     
      ! Get lon. center coordinate Array from Grid
      call ESMF_GridGetCoord(grid, coordDim=1, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           array=centerCoordArray, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      

      ! Get lon. corner coordinate Array from Grid
      call ESMF_GridGetCoord(grid, coordDim=1, &
           staggerloc=ESMF_STAGGERLOC_CORNER, &
           array=cornerCoordArray, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      
      
      ! Get localDECount from center Array
      ! center and corner localDECount should match
      ! since they are from the same Grid
      call ESMF_ArrayGet(centerCoordArray,  &
           localDECount=localDECount, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Create halo array
      call createAndFillHaloArray(centerCoordArray, centerCoordArraywHalo, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Loop over DEs calculating corners
      do lDE=0,localDECount-1
      
         ! Get center coord array
         call ESMF_ArrayGet(centerCoordArrayWHalo, localDE=lDE, &
              farrayPtr=centerWHaloPtr, rc=localrc) 
         if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
              NCF_CONTEXT, rcToReturn=rc)) return

         ! Get corner coord array
         call ESMF_ArrayGet(cornerCoordArray, localDE=lDE, &
              farrayPtr=cornerPtr, rc=localrc) 
         if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
              NCF_CONTEXT, rcToReturn=rc)) return         

         ! Debug Output
#if 0
         write(*,*) "center bnds=",lbound(centerWHaloPtr,1),ubound(centerWHaloPtr,1)
         write(*,*) "corner bnds=",lbound(cornerPtr,1),ubound(cornerPtr,1)
#endif

         ! Get bounds of DE
         lbnd=lbound(cornerPtr,1)
         ubnd=ubound(cornerPtr,1)

         ! If it's the bottom, then use -180.0 and adjust lower bound
         if (lbnd == 1) then
            cornerPtr(1)=-180.0_ESMF_KIND_R8
            lbnd = 2
         endif
         
         ! For the rest, Loop calculating corner coords as half way between centers
         do i=lbnd, ubnd
            cornerPtr(i)=(centerWHaloPtr(i-1)+centerWHaloPtr(i))/2.0_ESMF_KIND_R8
         enddo

         ! Debug output
#if 0
         lbnd=lbound(cornerPtr,1)
         ubnd=ubound(cornerPtr,1)
         write(*,*) "corner[",lbnd,"]=",cornerPtr(lbnd)
         write(*,*) "corner[",lbnd+1,"]=",cornerPtr(lbnd+1)
         write(*,*) "corner[",ubnd-1,"]=",cornerPtr(ubnd-1)
         write(*,*) "corner[",ubnd,"]=",cornerPtr(ubnd)
#endif

#if 0
         lbnd=lbound(centerWHaloPtr,1)
         ubnd=ubound(centerWHaloPtr,1)
         write(*,*) "center[",lbnd,"]=",centerWHaloPtr(lbnd)
         write(*,*) "center[",lbnd+1,"]=",centerWHaloPtr(lbnd+1)
         write(*,*) "center[",ubnd-1,"]=",centerWHaloPtr(ubnd-1)
         write(*,*) "center[",ubnd,"]=",centerWHaloPtr(ubnd)
#endif
         
      enddo

! XMRKX !

      ! Get rid of Halo Array
      call ESMF_ArrayDestroy(centerCoordArrayWHalo, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return         
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine calcCornerLongCoordsFromCenter

    
    ! Calculate corner latitude coordinates from center latitude coordinates
    ! Only supports 1D R8 coordinate arrays
    subroutine calcCornerLatCoordsFromCenter(grid, rc)
      
      ! Arguments
      type(ESMF_Grid) :: grid
      integer, intent(out),optional :: rc      

      ! Local variables
      integer :: localrc 
      type(ESMF_Array) :: centerCoordArray
      type(ESMF_Array) :: centerCoordArrayWHalo
      integer :: lDE,localDECount
      real(ESMF_KIND_R8), pointer :: centerwHaloPtr(:)
      type(ESMF_Array) :: cornerCoordArray
      real(ESMF_KIND_R8), pointer :: cornerPtr(:)
      integer :: i,lbnd,ubnd
      integer :: maxIndex(2)
      integer :: max_lat_ubnd
      

      ! Get Grid max upper bound for lat corners
      call ESMF_GridGet(grid, &
           tile=1, & ! Only support 1 Tile
           staggerloc=ESMF_STAGGERLOC_CORNER, &
           maxIndex=maxIndex, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Latitude is 2nd index
      max_lat_ubnd=maxIndex(2)
      
           
      ! Get lat. center coordinate Array from Grid
      call ESMF_GridGetCoord(grid, coordDim=2, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           array=centerCoordArray, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      

      ! Get lat. corner coordinate Array from Grid
      call ESMF_GridGetCoord(grid, coordDim=2, &
           staggerloc=ESMF_STAGGERLOC_CORNER, &
           array=cornerCoordArray, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      
      
      ! Get localDECount from center Array
      ! center and corner localDECount should match
      ! since they are from the same Grid
      call ESMF_ArrayGet(centerCoordArray,  &
           localDECount=localDECount, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Create halo array
      call createAndFillHaloArray(centerCoordArray, centerCoordArraywHalo, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Loop over DEs calculating corners
      do lDE=0,localDECount-1
      
         ! Get center coord array
         call ESMF_ArrayGet(centerCoordArrayWHalo, localDE=lDE, &
              farrayPtr=centerWHaloPtr, rc=localrc) 
         if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
              NCF_CONTEXT, rcToReturn=rc)) return

         ! Get corner coord array
         call ESMF_ArrayGet(cornerCoordArray, localDE=lDE, &
              farrayPtr=cornerPtr, rc=localrc) 
         if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
              NCF_CONTEXT, rcToReturn=rc)) return         

#if 0
         ! DEBUG OUTPUT
         write(*,*) "center bnds=",lbound(centerWHaloPtr,1),ubound(centerWHaloPtr,1)
         write(*,*) "corner bnds=",lbound(cornerPtr,1),ubound(cornerPtr,1)
#endif
         
         ! Get bounds of DE
         lbnd=lbound(cornerPtr,1)
         ubnd=ubound(cornerPtr,1)

         ! If it's the bottom, then use -90.0 and adjust lower bound
         if (lbnd == 1) then
            cornerPtr(1)=-90.0_ESMF_KIND_R8
            lbnd = lbnd + 1
         endif

         ! If it's the top, then use 90.0 and adjust upper bound down
         if (ubnd == max_lat_ubnd) then
            cornerPtr(ubnd)=90.0_ESMF_KIND_R8
            ubnd = ubnd - 1
         endif
         
         ! For the rest, Loop calculating corner coords as half way between centers
         do i=lbnd, ubnd
            cornerPtr(i)=(centerWHaloPtr(i-1)+centerWHaloPtr(i))/2.0_ESMF_KIND_R8
         enddo

         ! DEbug output
#if 0
         lbnd=lbound(cornerPtr,1)
         ubnd=ubound(cornerPtr,1)
         write(*,*) "corner[",lbnd,"]=",cornerPtr(lbnd)
         write(*,*) "corner[",lbnd+1,"]=",cornerPtr(lbnd+1)
         write(*,*) "corner[",ubnd-1,"]=",cornerPtr(ubnd-1)
         write(*,*) "corner[",ubnd,"]=",cornerPtr(ubnd)
#endif

#if 0
         lbnd=lbound(centerWHaloPtr,1)
         ubnd=ubound(centerWHaloPtr,1)
         write(*,*) "center[",lbnd,"]=",centerWHaloPtr(lbnd)
         write(*,*) "center[",lbnd+1,"]=",centerWHaloPtr(lbnd+1)
         write(*,*) "center[",ubnd-1,"]=",centerWHaloPtr(ubnd-1)
         write(*,*) "center[",ubnd,"]=",centerWHaloPtr(ubnd)
#endif
         
      enddo

! XMRKX !

      ! Get rid of Halo Array
      call ESMF_ArrayDestroy(centerCoordArrayWHalo, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return         
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine calcCornerLatCoordsFromCenter

    ! Create a copy of an Array that has a +1 halo on the lower side of each DE
    ! NOTE: THIS WILL CURRENTLY NOT WORK SINCE HALO IS SMM AND SMM DOESN'T WORK ON
    !       FACTORIZED ARRAYS. FIX IT BY USING THE CREATEFACTORIZEDDISTGRID METHOD
    !       TO MAKE A NEW DISTGRID FOR wHALOARRAY
    subroutine createAndFillHaloArray(origArray, wHaloArray, rc)
      ! Arguments
      type(ESMF_Array) :: origArray
      type(ESMF_Array) :: wHaloArray
      integer, intent(out),optional :: rc      

      ! Local variables
      integer :: localrc
      type(ESMF_Index_Flag) :: indexflag
      type(ESMF_Distgrid) :: distgrid
      integer :: distgridToArrayMap(2)
      integer :: totalLWidth(1)
      integer :: dimCount
      integer :: lDE, localDECount
      real(ESMF_KIND_R8), pointer :: origPtr(:),wHaloPtr(:)
      integer :: lbnd,ubnd
      type(ESMF_RouteHandle) :: haloHandle
      
      ! Get info from original Array
      call ESMF_ArrayGet(origArray, &
           distgrid=distgrid, &
           indexflag=indexflag, &
           dimCount=dimCount, &
           distGridToArrayMap=distgridToArrayMap, &
           localDECount=localDECount, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return         

      ! Only 2D grids supported right now
      if (dimCount .ne. 2) then
         call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
              msg="Only 2D Grids currently supported.", &
              NCF_CONTEXT, rcToReturn=rc)
         return
      endif

      ! Set lower halo
      totalLWidth(:)=1
      
      ! Create new Array
      wHaloArray=ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, &
           indexflag=indexflag, distgridToArrayMap=distgridToArrayMap, &
           totalLWidth=totalLWidth, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return         

      ! Copy data
      do lDE=0,localDECount-1

         ! Get orig array data
         call ESMF_ArrayGet(origArray, localDE=lDE, &
              farrayPtr=origPtr, rc=localrc) 
         if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
              NCF_CONTEXT, rcToReturn=rc)) return

         ! Get corner coord array
         call ESMF_ArrayGet(wHaloArray, localDE=lDE, &
              farrayPtr=wHaloPtr, rc=localrc) 
         if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
              NCF_CONTEXT, rcToReturn=rc)) return         

         ! Get bounds of DE
         lbnd=lbound(origPtr,1)
         ubnd=ubound(origPtr,1)

         ! Init
         wHaloPtr(:)=0.0
         
         ! Copy
         wHaloPtr(lbnd:ubnd)=origPtr(lbnd:ubnd)
      enddo


      ! Parallel Halo of data
      call ESMF_ArrayHaloStore(wHaloArray, routeHandle=haloHandle, rc=localrc) 
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      call ESMF_ArrayHalo(wHaloArray, routeHandle=haloHandle, rc=localrc) 
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      call ESMF_ArrayHaloRelease(routeHandle=haloHandle, rc=localrc) 
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return         

      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine createAndFillHaloArray


    
        
!
!  check NetCDF file error code
!
function NetCDFCheckError (ncStatus, methodName, fileName, lineNo, errmsg, rc)
#undef  NCF_METHOD
#define NCF_METHOD "NetCDFCheckError()"
  ! Return value
  logical                       :: NetCDFCheckError

  ! Arguments
  integer,          intent(in)  :: ncStatus
  character(len=*), intent(in)  :: methodName
  character(len=*), intent(in)  :: fileName
  integer,          intent(in)  :: lineNo
  character(len=*), intent(in)  :: errmsg
  integer, intent(out),optional :: rc

  ! Internal variables
  integer, parameter :: nf90_noerror = 0

  ! Init return
  NetCDFCheckError = .FALSE.

#ifdef ESMF_NETCDF
  if ( ncStatus .ne. nf90_noerror) then
     call ESMF_LogWrite (  &
          msg="netCDF Error: " // trim (errmsg) // ": " // trim (nf90_strerror(ncStatus)),  &
          logmsgFlag=ESMF_LOGMSG_ERROR, &
          line=lineNo, file=fileName, method=methodName)
     print '("NetCDF Error: ", A, " : ", A)', &
          trim(errmsg),trim(nf90_strerror(ncStatus))
     call ESMF_LogFlush()
     if (present(rc)) rc = ESMF_FAILURE
     NetCDFCheckError = .TRUE.
  else
     if (present(rc)) rc = ESMF_SUCCESS
     return
  end if
#else
  if (ESMF_LogFoundError(ESMF_RC_LIB_NOT_PRESENT, &
       msg="- ESMF_NETCDF not defined when lib was compiled", &
       ESMF_CONTEXT, rcToReturn=rc)) return
#endif
  
end function NetCDFCheckError


    
end module
