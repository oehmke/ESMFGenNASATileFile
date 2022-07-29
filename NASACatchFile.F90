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

      ! Fill coordinates from file
      call readCenterCoordsFromNASACatchFile(createGridFromNASACatchFile, &
           filename, localrc)
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
    subroutine readCenterCoordsFromNASACatchFile(grid, filename, rc)
      ! Arguments
      type(ESMF_Grid) :: grid
      character(len=*), intent(in)  :: filename
      integer, intent(out),optional :: rc      

      ! Local variables
      integer :: localrc 
      type(ESMF_Array) :: lonCoordArray


      ! Get lon. coordinate Array from Grid
      call ESMF_GridGetCoord(grid, coordDim=1, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           array=lonCoordArray, rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      
      
      ! Read lon. coordinates into Array
      call ESMF_ArrayRead(lonCoordArray, &
           fileName=filename, variableName="longitude", &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NCF_PASSTHRU, &
           NCF_CONTEXT, rcToReturn=rc)) return      
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine readCenterCoordsFromNASACatchFile
  
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
