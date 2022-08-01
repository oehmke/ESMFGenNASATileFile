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

  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

    ! TODO: CONSIDER ADDING NCF_ TO THE FRONT OF ALL PUBLIC METHODS IN THIS MODULE,
    !       AND THEN SHORTEN THEIR NAMES TO GET RID OF NASACatchFile

    ! TODO: MAKE ALL NON-PUBLIC METHODS PRIVATE!

    
    ! Creates a Grid to represent the index space and coordinates of
    ! a NASA catchment file
    function createGridFromNASACatchFile(fileName, rc)
#undef  NCF_METHOD
#define NCF_METHOD "createGridFromNASACatchFile()"

      ! Return value
      type(ESMF_Grid) :: createGridFromNASACatchFile
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

      
      ! Debug output
      write(*,*) "Grid size=",lonSize," x ",latSize

      ! Create Grid
      createGridFromNASACatchFile=ESMF_GridCreate1PeriDim( &
           maxIndex=(/lonSize,latSize/), &
           coordSys=ESMF_COORDSYS_SPH_DEG, &
           coordDep1=(/1/),  coordDep2=(/2/), &
           indexflag=ESMF_INDEX_GLOBAL,  rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Add Center coords
      call ESMF_GridAddCoord(createGridFromNASACatchFile, &
           staggerloc=ESMF_STAGGERLOC_CENTER, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Fill longitude center coordinates from file
      call readCenterCoordsFromNASACatchFile(createGridFromNASACatchFile, &
           filename, coordDim=1, varname="longitude", rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return

      ! Fill latitude center coordinates from file
      call readCenterCoordsFromNASACatchFile(createGridFromNASACatchFile, &
           filename, coordDim=2, varname="latitude", rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      
      
      ! Add Corner coords
      call ESMF_GridAddCoord(createGridFromNASACatchFile, &
           staggerloc=ESMF_STAGGERLOC_CORNER, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return
      
      ! Close file
      ncStatus = nf90_close(ncid)
      if (NetCDFCheckError (ncStatus, &
           NCF_METHOD,  &
           NCF_SRCLINE, &
           "Closing file "//trim(filename), &
           rc)) return
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end function createGridFromNASACatchFile


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
