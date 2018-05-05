c=======================================================================
c
c  GenWFL:  Generate WPhotpmc Fole Lists by locating all epoch
c           directorie under a specified tile directory, read each
c           binary FITS tablem to determine whether ascending or
c           descending or both, and create the files for WPhotpmc
c           to access unWISE coadd images
c      TBD: what to do with "mixed" frames lists
c           SKIP THEM
c           whether to rescale "std" images for W1 3-band-cryo epochs
c           DO UPON REQUEST
c           whether to rename W1 cryo epochs as W3
c           YES when "-1b" is specified
c  vsn 1.0  B71214: initial version
c      1.1  B80111: made temporary-file pathname specifiable
c      1.2  B80112: added percentage requirement to define Asce/Desc;
c                   added PA range to epoch/scandir xref table;
c                   included full-depth coadd & "n" image for W1 -> W3
c           B80117: installed PA-range subroutine ChkPAhist
c           B80119: installed randomized work-file name
c           B80223: made FBT file names 500 characters long
c           B80226: changed "end of cryo" to "end of cryo PSF"
c           B80314: put error stop if temp file can't be created
c           B80327: fixed bug that wrote to unit 12 despite GotOX = F
c           B80504: fixed formatting bug in the epochs table
c
c=======================================================================
c
      Character*500 TileDirNam, FLnamA, FLnamD, XrefNam, cpLine,
     +              Std3bcNam, W3Nam, TmpNam, WorkNam, TempDir,
     +              FBTnam1, FBTnam2, TileNam
      Character*62  XRefLine, TmpStr
      Character*11  Vsn, Flag, EpNam
      Character*8   cDate, cTime
      Character*1   W1char, W2char, OKchar
      Integer*4     IArgC, NArg, NArgs, LNBlnk, FileID, system, status,
     +              nEpochs, nEpA, nEpD, nEpMx, readwrite, blocksize,
     +              nRows, nCols, hdutype, k, felem, nelems, nPA, nAsce,
     +              nMJD, nIncd, n, n0, n1, nDec, nRA, ScanDir,
     +              nUsed, n4bc, n3bc, npc, nMisMch, Access, nEpAM,
     +              nEpDM, PAhist(3600), nPct, nOutA, nOutD
      Real*4        d2r, RA, Dec, long, lat, dot, fac, fwdlong, sunlong,
     +              pa, pa3, RS3bcFac, rs3bc, fsd, fscan, PAmin, PAmax
      Real*8        dRA, dDec, MJD, MJD0, MJD3, MJD4
      Logical*4     GotTileDir, dbg, NEO3, IzCryo, Iz3band, GotW1,
     +              GotOA, GotOD, GotOX, opt1b, anynull, DoRS3b, GotW2,
     +              Wrt20, Wrt21, OK
      byte          incd, asce
c
      data Vsn/'1.2  B80504'/, GotTileDir,dbg/2*.false./, nMisMch/0/,
     +     nEpochs,nEpA,nEpD,nEpMx,nEpAM,nEpDM/6*0/, TempDir/'.'/,
     +     GotOA,GotOD,GotOX/3*.false./, nOutA,nOutD/2*0/
     +     opt1b/.false./, d2r/1.745329252e-2/, fsd/0.95/,
     +     rs3bc/1.2/, DoRS3b/.false./, WorkNam/'genwfltmp'/
      data MJD0/57467.6875/,     ! 3/20/2016 04 30 00; used for sunlong
     +     fac/0.9856101/,       ! 360/365.256     
     +     MJD4/55414.932/,      ! August 6, 2010, end of 4-band cryo
     +     MJD3/55480.0/         ! JD = 2455480.0, October 11, 2010
c                                ! Peter's Preference, end of "cryo PSF"
c    +     MJD3/55469.277509259/ ! real end of 3-band cryo
c
      common /vdt/ cdate,ctime,vsn
c
c=======================================================================
c
      NARgs = IArgC()
1     if (NARgs .lt. 6) then
        print *,'GenWFL vsn ', Vsn
        print *
        print *,'Usage:   genwfl <flags specifications>'
        print *
        print *,'Where the REQUIRED flags and specifications are:'
        print *,' -t   pathname for the tile directory (RaDecID)'
        print *,' -oa  output file name for ascending frames list'
        print *,' -od  output file name for descending frames list'
        print *
        print *,'The OPTIONAL flags and specifications are:'
        print *,' -1a  support option 1a (default)'
        print *,' -1b  support option 1b'
        print *,' -rs  #.# rescaling factor for W1 3-band-cryo std'
        print *,' -mf  minimum fraction of scans in one direction to'
        print *,'      define an epoch as Asce or Desc (0.95)'
        print *,' -ox  output file name for epoch/scandir xref table'
        print *,' -td  temporary work file directory name (.)'
        print *,' -d   debug mode'
        print *
        stop
      end if
c
c=======================================================================
c
      NArg = 0
2     NArg = NArg + 1
      call GetArg(NArg,Flag)
      call UpCase(Flag)
c                                      ! Tile Directory Name
      If (Flag .eq. '-T') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TileDirNam)
        GotTileDir = .true.
        if (dbg) print *,'Tile Directory: ',
     +                    TileDirNam(1:lnblnk(TileDirNam))
      Else if (Flag .eq. '-OA') then    ! Ascending Frames List Name
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,FLNamA)
        GotOA = .true.
        if (dbg) print *,'Ascending Frames List name: ',
     +                    FLNamA(1:lnblnk(FLNamA))
      Else if (Flag .eq. '-OD') then    ! Descending Frames List Name
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,FLNamD)
        GotOD = .true.
        if (dbg) print *,'Descending Frames List name: ',
     +                    FLNamD(1:lnblnk(FLNamD))
      Else if (Flag .eq. '-OX') then    ! Xref Name
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,XrefNam)
        GotOX = .true.
        if (dbg) print *,'Epoch/ScanDir xref table name: ',
     +                    XrefNam(1:lnblnk(XrefNam))
      Else if (Flag .eq. '-RS') then    ! Rescale factor for W1 3-band cryo std
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3007) rs3bc
        DoRS3b = .true.
        if (dbg) print *,'Rescale factor for W1 3-band cryo std: ',
     +                    rs3bc
      Else if (Flag .eq. '-MF') then    ! Min fraction of given scan dir
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3007) fsd
        if (dbg) print *,'Min fraction of given scan dir: ', fsd
      Else if (Flag .eq. '-1A') then
        opt1b = .false.
      Else if (Flag .eq. '-1B') then
        opt1b = .true.
      Else if (Flag .eq. '-D') then
        dbg = .true.
      Else if (Flag .eq. '-TD') then    ! Temporary Work Directory Name
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TempDir)
        if (dbg) print *,'Work file directory name: ',
     +                    TempDir(1:lnblnk(TempDir))
        k = lnblnk(TempDir)
        if (TempDir(k:k) .ne. '/') TempDir(k+1:k+1) = '/'
      end if
      If (NArg .lt. NArgs) Go to 2
c
c=======================================================================
c
      if (.not.GotTileDir) then
        print *,'ERROR: Full-depth directory name not specified'
        NArgs = 0
        go to 1
      end if
      if (.not.GotOA) then
        print *,'ERROR: Ascending Frames List name not specified'
        NArgs = 0
        go to 1
      end if
      if (.not.GotOD) then
        print *,'ERROR: Descending Frames List name not specified'
        NArgs = 0
        go to 1
      end if
c
c-----------------------------------------------------------------------
c      
      call signon('genwfl')
c
      k = lnblnk(TileDirNam)
      If ((TileDirNam(k:k) .eq. '/') .or. (TileDirNam(k:k) .eq. '/'))
     +     TileDirNam(k:k) = ' '
      TileNam = TileDirNam
5     k = index(TileNam,'/')
      if (k .eq. 0) k = index(TileNam,'\')
      if (k .gt. 0) then
        TileNam = TileNam(k+1:lnblnk(TileNam))
        go to 5
      end if
      if (dbg) print *,'Tile name: ', TileNam(1:lnblnk(TileNam))
c
      call GetTmpNam(WorkNam,TempDir,OK,dbg)
      if (.not.OK) then
	    print *,'Unable to create temporary work file'
        call exit(64)
      end if
      status = system('ls -1d '//TileDirNam(1:lnblnk(TileDirNam))
     +      //'/e* > '//WorkNam(1:lnblnk(WorkNam)))
      if (status .ne. 0) then
        print *,'ERROR: attempt to get list of epochs failed'
        print *,'       system status =', status
        print *,'Execution terminated'
        call exit(64)
      end if
      open (10, file = WorkNam)
      felem  = 1
      nelems = 1
      if (GotOX) then
        open(12, file = XRefNam)
        write(12,'(a)') '| Ep|   ScanDir   |   CryoType    |'
     +              //'PAmin | PAmax| dPA  |Pct|OK|'
      end if
      open  (20, file = FLNamD)
      write (20,'(a)')
     +   '|  path              |      base     | b1 | b2 | b3 | b4 |'
      write (20,'(a)')
     +   '|   c                |       c       |  i |  i |  i |  i |'
      write (20,'(a)') '|                    |'
      write (20,'(a)') '|                    |'
      open  (21, file = FLNamA)
      write (21,'(a)')
     +   '|  path              |      base     | b1 | b2 | b3 | b4 |'
      write (21,'(a)')
     +  '|   c                |       c       |  i |  i |  i |  i |'
      write (21,'(a)') '|                    |'
      write (21,'(a)') '|                    |'
c
c-----------------------------------------------------------------------
c                                      ! Loop through all epochs
10    read (10, '(a)', end = 500) TmpNam
12    k = index(TmpNam,'/')
      if (k .eq. 0) k = index(TmpNam,'\')
      if (k .gt. 0) then
        TmpNam = TmpNam(k+1:lnblnk(TmpNam))
        go to 12
      end if
      EpNam = TmpNam
      FBTnam1 = TileDirNam(1:lnblnk(TileDirNam))//'/'
     +       //EpNam(1:lnblnk(EpNam))//'/unwise-'
     +       //TileNam(1:lnblnk(TileNam))//'-w1-frames.fits'
      FBTnam2 = TileDirNam(1:lnblnk(TileDirNam))//'/'
     +       //EpNam(1:lnblnk(EpNam))//'/unwise-'
     +       //TileNam(1:lnblnk(TileNam))//'-w2-frames.fits'
      if (dbg) print *,'W1 FBT name: ',FBTnam1(1:lnblnk(FBTnam1))
      if (dbg) print *,'W2 FBT name: ',FBTnam2(1:lnblnk(FBTnam2))
      GotW1 = (Access(FBTnam1(1:lnblnk(FBTnam1)),' ') .eq. 0)
      GotW2 = (Access(FBTnam2(1:lnblnk(FBTnam2)),' ') .eq. 0)
      if (.not.(GotW1 .or. GotW2)) then
        print *,'WARNING: FBT Files not found:'
        print *,FBTnam1(1:lnblnk(FBTnam1))
        print *,FBTnam2(1:lnblnk(FBTnam2))
        go to 3333
      end if
      nUsed = 0
      n0    = 0                       ! no. of desc-scan frames
      n1    = 0                       ! no. of asce-scan frames
      n4bc  = 0                       ! no. of 4-band-cryo frames
      n3bc  = 0                       ! no. of 3-band-cryo frames
      npc   = 0                       ! no. of post-cryo frames
      IzCryo  = .false.
      Iz3band = .false.
      PAmin   =  9999.99
      PAmax   = -9999.99
      PAhist  = 0
c      
      status = 0
      call ftgiou(FileID,status)         ! get unit number
      readwrite = 0                      ! open FITS file
      if (GotW1) then
        call ftopen(FileID,FBTnam1,readwrite,blocksize,status)
      else
        call ftopen(FileID,FBTnam2,readwrite,blocksize,status)
      end if
      if (status .ne. 0) then
        if (GotW1) then
          write(6,'(a)') 'ERROR: Could not open '//trim(FBTnam1)
        else
          write(6,'(a)') 'ERROR: Could not open '//trim(FBTnam2)
        end if
        go to 3333
      endif
      nEpochs = nEpochs + 1
      call ftmahd(FileID,2,hdutype,status) ! move to header #2
      if (dbg) print *,'hdutype for nhdu-2, status:     ',
     +                  hdutype, status
      status = 0
      NEO3   = .false.
      call ftgnrw(FileID,nrows,status)     ! get #rows
      if (status .ne. 0) go to 3001
      if (dbg) print *,'No. of rows in table and status:', nrows, status
      call ftgncl(FileID,ncols,status)     ! get #cols
      if (status .ne. 0) go to 3002
      if (dbg) print *,'No. of cols in table and status:', ncols, status
      call ftgcno(FileID, .false., 'pa      ', nPA,  status)    ! row# for PA
      NEO3 = (status .eq. 0)
      status = 0
      call ftgcno(FileID, .false., 'ascending', nAsce, status)  ! row# for ascending
      NEO3 = NEO3 .and. (status .eq. 0)
      if (dbg .or.(nEpochs .eq. 1)) then
        If (NEO3) then
          print *,'FBT is NEO3 format'
        else
          print *,'FBT is NEO2 format'
        end if
      end if
      if (dbg .and. NEO3) then
        print *,'Col. no and status for pa:       ', nPA, status
        print *,'Col. no and status for ascending:', nAsce, status
      end if
      status = 0
      call ftgcno(FileID, .false., 'mjd     ', nMJD, status)  ! row# for mjd
      if (dbg) print *,'Col. no and status for mjd:     ', nMJD, status
      if (status .ne. 0) go to 3005
      call ftgcno(FileID, .false., 'included', nIncd, status) ! row# for included
      if (status .ne. 0) go to 3006
      if (dbg) print *,'Col. no and status for included:', nIncd, status
      call ftgcno(FileID, .false., 'ra      ', nRA,  status)  ! row# for RA
      if (status .ne. 0) go to 3003
      if (dbg) print *,'Col. no and status for ra:      ', nRA, status
      call ftgcno(FileID, .false., 'dec     ', nDec, status)  ! row# for Dec
      if (status .ne. 0) go to 3004
      if (dbg) print *,'Col. no and status for dec:     ', nDec, status
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                 ! Loop through all records for this epoch    
      do 100 n = 1, nrows
c      
        k = 1
        call ftgcvd(FileID,nRA,n,felem,nelems,0.0d0,dRA,anynull,status)
        if (status .ne. 0) go to 3000
        ra = dRA
        k = 2
        call ftgcvd(FileID,nDec,n,felem,nelems,0.0d0,dDec,anynull,status)
        if (status .ne. 0) go to 3000
        dec = dDec
        k = 3
        call ftgcvd(FileID,nMJD,n,felem,nelems,0.0d0,MJD,anynull,status)
        if (status .ne. 0) go to 3000
        k = 4
        call ftgcvb(FileID,nIncd,n,felem,nelems,0,incd,anynull,status)
        if (status .ne. 0) go to 3000
c        
        if (dbg .and. (n .eq. 1)) then
          print *
          print *,'        RA              Dec              Long'
     +         //'            Lat             SunLong          FwdLong'
     +         //'            ScanDir    PA               MJD'
        end if
c
        if (incd .eq. 0) go to 100       
        call cel2ec(ra, dec, long, lat)
        call cel2pa(ra, dec, pa)
        sunlong = (mjd - mjd0)*fac
20      if (sunlong .gt. 360.0) then
          sunlong = sunlong - 360.0
          go to 20
        end if
30      if (sunlong .lt. 0.0) then
          sunlong = sunlong + 360.0
          go to 30
        end if
        fwdlong = sunlong - 90.0
        if (cos(d2r*(long-fwdlong)) .gt. 0.0) then   ! looking forward
          ScanDir = 0                                ! desc
          n0 = n0 + 1
        else                                         ! looking backward
          ScanDir = 1                                ! asce
          pa = pa + 180.0
          n1 = n1 + 1
        end if
        if (NEO3) then
          k = 5
          call ftgcve(FileID,nPA,n,felem,nelems,0.0d0,pa3,anynull,status)
          if (status .ne. 0) go to 3000
          k = 6
          call ftgcvb(FileID,nAsce,n,felem,nelems,0.0d0,asce,anynull,status)
          if (status .ne. 0) go to 3000
          if (asce .ne. ScanDir) then
            print *,'WARNING: Scan direction differs for NEO2 and NEO3'
            print *,'         NEO3:', asce,'; NEO2:', ScanDir
            nMisMch = nMisMch + 1
          end if
          dot = cos(d2r*pa)*cos(d2r*pa3)
     +        + sin(d2r*pa)*sin(d2r*pa3)
          if (dot .lt. 0.9961947) then
            print *,'WARNING: NEO2/NEO3 PAs differ by > 5 deg'
            print *,'         NEO2:', pa
            print *,'         NEO3:', pa3
          end if
          pa = pa3
        end if
        if (pa .gt. 360.0) pa = pa - 360.0
        if (pa .lt.   0.0) pa = pa + 360.0
        k = NInt(10.0*pa+0.5)
        if (k .lt. 1)   k = 1
        if (k .gt. 3600) k = 3600
        PAhist(k) = PAhist(k) + 1
c
        if (dbg) print *, ra, dec, long, lat,
     +                    sunlong, fwdlong, ScanDir, pa, mjd
        nUsed = nUsed + 1
        if (pa .lt. PAmin) PAmin = pa
        if (pa .gt. PAmax) PAmax = pa
c
        if (mjd .le. MJD4) then
          n4bc   = n4bc + 1
          IzCryo = .true.
        else if (mjd .le. MJD3) then
          n3bc   = n3bc + 1
          IzCryo  = .true.
          Iz3band = .true.
        else
          npc = npc + 1
        end if
c
100   continue      
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Finish processing this epoch    
      if (DoRS3b .and. Iz3Band .and. GotW1) then
        RS3bcFac = rs3bc*float(n3bc)/float(nUsed)
        Std3bcNam = TileDirNam(1:lnblnk(TileDirNam))//'/'
     +        //EpNam(1:lnblnk(EpNam))//'/unwise-'
     +        //TileNam(1:lnblnk(TileNam))
     +        //'-w1-std-m.fits'
        call Rescale(Std3bcNam,RS3bcFac)
      end if
c
      if (GotW1) then
        W1char = '1'
      else
        W1char = '0'
      end if
      if (GotW2) then
        W2char = '1'
      else
        W2char = '0'
      end if
c
      Wrt20 = (n0 .gt. 0)
      Wrt21 = (n1 .gt. 0)
      XRefLine = EpNam
      XRefline(62:62) = 'N'
      if ((n0 .gt. 0) .and. (n1 .gt. 0)) then
        if (n1 .gt. n0) then
          XRefLine(9:18) = 'Mixed/Asce'
          Wrt20 = .false.
          fscan = float(n1)/float(n0+n1)
          if (fscan .lt. fsd) Wrt21 = .false.
          nEpAM = nEpAM + 1
        end if
        if (n1 .lt. n0) then
          XRefLine(9:18) = 'Mixed/Desc'
          Wrt21 = .false.
          fscan = float(n0)/float(n0+n1)
          if (fscan .lt. fsd) Wrt20 = .false.
          nEpDM = nEpDM + 1
        end if
        if (n1 .eq. n0) then
          XRefLine(9:18) = 'Mixed     '
          nEpMx = nEpMx + 1
          Wrt20 = .false.
          Wrt21 = .false.
        end if
        print *,'WARNING: ',EpNam(1:lnblnk(EpNam)),
     +          ' has both Asce and Desc frames'
      end if
c      
      if (Wrt20) then
        if (n1 .eq. 0) then
          XRefLine(9:12) = 'Desc'
          nEpD = nEpD + 1
        end if
        if (opt1b .and. IzCryo) then
          write (20,'(a)') TileDirNam(1:lnblnk(TileDirNam))//'/'
     +        //EpNam(1:lnblnk(EpNam))//'/ unwise-'
     +        //TileNam(1:lnblnk(TileNam))//' 0 '
     +        //W2char//' '//W1char//' 0'
          nOutD = nOutD + 1
          XRefline(62:62) = 'Y'
          if (GotW1) then
            W3Nam = TileDirNam(1:lnblnk(TileDirNam))
     +          //'/'//EpNam(1:lnblnk(EpNam))//'/unwise-'
     +          //TileNam(1:lnblnk(TileNam))
     +          //'-w3-img-m.fits'
            if (Access(W3Nam(1:LNBlnk(W3Nam)),' ') .ne. 0) then
              cpLine = 'cp '//TileDirNam(1:lnblnk(TileDirNam))
     +            //'/'//EpNam(1:lnblnk(EpNam))//'/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w1-img-m.fits '//W3Nam
              k = system(cpLine)
              print *, cpLine(1:lnblnk(cpLine))
              cpLine = 'cp '//TileDirNam(1:lnblnk(TileDirNam))//'/'
     +            //EpNam(1:lnblnk(EpNam))//'/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w1-std-m.fits '//TileDirNam(1:lnblnk(TileDirNam))
     +            //'/'//EpNam(1:lnblnk(EpNam))//'/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w3-std-m.fits'
              k = system(cpLine)
              print *, cpLine(1:lnblnk(cpLine))
            end if
            W3Nam = TileDirNam(1:lnblnk(TileDirNam))
     +          //'/Full/unwise-'
     +          //TileNam(1:lnblnk(TileNam))
     +          //'-w3-img-m.fits'
            if (Access(W3Nam(1:LNBlnk(W3Nam)),' ') .ne. 0) then
              cpLine = 'cp '//TileDirNam(1:lnblnk(TileDirNam))
     +            //'/full/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w1-img-m.fits '//W3Nam
              k = system(cpLine)
              print *, cpLine(1:lnblnk(cpLine))
              W3Nam = TileDirNam(1:lnblnk(TileDirNam))
     +            //'/Full/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w3-n-m.fits'
              cpLine = 'cp '//TileDirNam(1:lnblnk(TileDirNam))
     +            //'/full/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w1-n-m.fits '//W3Nam
              k = system(cpLine)
              print *, cpLine(1:lnblnk(cpLine))
            end if
          end if
        else
          write (20,'(a)') TileDirNam(1:lnblnk(TileDirNam))//'/'
     +        //EpNam(1:lnblnk(EpNam))//'/ unwise-'
     +        //TileNam(1:lnblnk(TileNam))//' '
     +        //W1char//' '//W2char//' 0 0'
          nOutD = nOutD + 1
          XRefline(62:62) = 'Y'
        end if
      end if
c      
      if (Wrt21) then
        if (n0 .eq. 0) then
          XRefLine(9:12) = 'Asce'
          nEpA = nEpA + 1
        end if
        if (opt1b .and. IzCryo) then
          write (21,'(a)') TileDirNam(1:lnblnk(TileDirNam))//'/'
     +        //EpNam(1:lnblnk(EpNam))//'/ unwise-'
     +        //TileNam(1:lnblnk(TileNam))//' 0 '
     +        //W2char//' '//W1char//' 0'
          nOutA = nOutA + 1
          XRefline(62:62) = 'Y'
          if (GotW1) then
            W3Nam = TileDirNam(1:lnblnk(TileDirNam))
     +          //'/'//EpNam(1:lnblnk(EpNam))//'/unwise-'
     +          //TileNam(1:lnblnk(TileNam))
     +          //'-w3-img-m.fits'
            if (Access(W3Nam(1:LNBlnk(W3Nam)),' ') .ne. 0) then
              cpLine = 'cp '//TileDirNam(1:lnblnk(TileDirNam))
     +            //'/'//EpNam(1:lnblnk(EpNam))//'/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w1-img-m.fits '//W3Nam
              k = system(cpLine)
              print *, cpLine(1:lnblnk(cpLine))
              cpLine = 'cp '//TileDirNam(1:lnblnk(TileDirNam))//'/'
     +            //EpNam(1:lnblnk(EpNam))//'/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w1-std-m.fits '//TileDirNam(1:lnblnk(TileDirNam))
     +            //'/'//EpNam(1:lnblnk(EpNam))//'/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w3-std-m.fits'
              k = system(cpLine)
              print *, cpLine(1:lnblnk(cpLine))
            end if
            W3Nam = TileDirNam(1:lnblnk(TileDirNam))
     +          //'/Full/unwise-'
     +          //TileNam(1:lnblnk(TileNam))
     +          //'-w3-img-m.fits'
            if (Access(W3Nam(1:LNBlnk(W3Nam)),' ') .ne. 0) then
              cpLine = 'cp '//TileDirNam(1:lnblnk(TileDirNam))
     +            //'/full/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w1-img-m.fits '//W3Nam
              k = system(cpLine)
              print *, cpLine(1:lnblnk(cpLine))
              W3Nam = TileDirNam(1:lnblnk(TileDirNam))
     +            //'/Full/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w3-n-m.fits'
              cpLine = 'cp '//TileDirNam(1:lnblnk(TileDirNam))
     +            //'/full/unwise-'
     +            //TileNam(1:lnblnk(TileNam))
     +            //'-w1-n-m.fits '//W3Nam
              k = system(cpLine)
              print *, cpLine(1:lnblnk(cpLine))
            end if
          end if
        else
          write (21,'(a)') TileDirNam(1:lnblnk(TileDirNam))//'/'
     +        //EpNam(1:lnblnk(EpNam))//'/ unwise-'
     +        //TileNam(1:lnblnk(TileNam))//' '
     +        //W1char//' '//W2char//' 0 0'
          nOutA = nOutA + 1
          XRefline(62:62) = 'Y'
        end if
      end if
c
      OKchar = XrefLine(62:62)     
      if ((n4bc .gt. 0) .and. (n3bc .gt. 0)) then
c                                  '   CryoType    '
        XRefLine = XrefLine(1:20)//' 4&3-band cryo'
      else if ((npc .gt. 0) .and. (n3bc .gt. 0)) then
        XRefLine = XrefLine(1:19)//'3-band&post-cryo'
      else if (n4bc .gt. 0) then
        XRefLine = XrefLine(1:20)//' 4-band cryo'
      else if (n3bc .gt. 0) then
        XRefLine = XrefLine(1:20)//' 3-band cryo'
      else if (npc .gt. 0) then
        XRefLine = XrefLine(1:20)//' post-cryo'
      end if
      XrefLine(62:62) = OKchar
      if (dbg) print *,'PAmin, PAmax before ChkPAhist:', PAmin, PAmax
      call ChkPAhist(PAhist,PAmin,PAmax,dbg)
      if (dbg) print *,'PAmin, PAmax after ChkPAhist: ', PAmin, PAmax
      write (XrefLine(36:41),'(F6.1)') PAmin
      write (XrefLine(43:48),'(F6.1)') PAmax
      write (XrefLine(50:55),'(F6.1)') PAmax-PAmin
      nPct = 100*Max0(n0,n1)/(n0+n1)
      write (XrefLine(57:59),'(I3)') nPct
      if (GotOX) write (12,'(a)') XRefLine(1: lnblnk(XRefLine))
      if (dbg) print *, XRefLine
c
      if (dbg) then
        print *,'No. of rows read: ', nrows
        print *,'No. of rows used: ', nUsed
        print *,'No. ascending:',n1,'; no. descending:',n0
      end if
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
      go to 10
c
c-----------------------------------------------------------------------
c                                      ! Termination processing
500   close(10)
      if (.not.dbg) status = system('rm '//WorkNam(1:lnblnk(WorkNam)))
      print *,'No. epochs processed:        ', nEpochs
      print *,'No. pure  ascending epochs:  ', nEpA
      print *,'No. pure  descending epochs: ', nEpD
      print *,'No. mixed ascending epochs:  ', nEpAM
      print *,'No. mixed descending epochs: ', nEpDM
      print *,'No. evenly mixed epochs:     ', nEpMx
      print *,'No. ascending epochs output: ', nOutA
      print *,'No. descending epochs output:', nOutD
      if (nMisMch .gt. 0)
     +    print *,'No. NEO2/NEO3 ASCE/DESC mismatches:', nMisMch
c
      call signoff('genwfl')
      stop
c
c-----------------------------------------------------------------------
c                                      ! Error/warning messages
3000  print *,'ERROR reading table row no.', n,
     +        '; parameter', k,', epoch '//EpNam(1:lnblnk(EpNam)),
     +        '; status =',status
      go to 3333
3001  print *,'ERROR: unable to get @rows; status =',status
      go to 3333
3002  print *,'ERROR: unable to get #cols; status =',status
      go to 3333
3003  print *,'ERROR: ra not found; status =',status
      go to 3333
3004  print *,'ERROR: dec not found; status =',status
      go to 3333
3005  print *,'ERROR: mjd not found; status =',status
      go to 3333
3006  print *,'ERROR: included not found; status =',status
      go to 3333
3007  print *,'ERROR: bad value for "-rs":', TmpStr
      call exit(64)
c
3333  print *,'skipping this epoch'
      go to 10      
c      
      end
c      
c=======================================================================
c
      subroutine SignOn(pgmnam)
c
c *** signon- routine which provides sign-on and sign-off messages
c             (orig by John Fowler- mod by Howard McCallon-041214-SIRTF)
c
c     inputs:  pgmnam = program name                                 [call arg]
c
c     outputs: message to stdout
c
      character*(*) pgmnam
      character vsn*11,cdate*8,ctime*8,Fmt*11,FLen*4
      integer*4 onoff,jdate(3),jtime(3),lnblnk
      real*4    dummyt,second(2),etime
c
      common /vdt/ cdate,ctime,vsn
c##
      onoff = 1
c
c         i. obtain date
c
100   cdate = '00-00-00'
      call idate(jdate)    ! Linux call
c
      jdate(3) = mod(jdate(3), 100)
      write(cdate(1:2), '(i2)') jdate(2)
      write(cdate(4:5), '(i2)') jdate(1)
      write(cdate(7:8), '(i2)') jdate(3)
c
      if(cdate(4:4) .eq. ' ') cdate(4:4) = '0'
      if(cdate(7:7) .eq. ' ') cdate(7:7) = '0'
c
c         ii. obtain time
c
      ctime = '00:00:00'
      call itime(jtime)
      write(ctime(1:2), '(i2)') jtime(1)
      write(ctime(4:5), '(i2)') jtime(2)
      write(ctime(7:8), '(i2)') jtime(3)
c
      if(ctime(4:4) .eq. ' ') ctime(4:4) = '0'
      if(ctime(7:7) .eq. ' ') ctime(7:7) = '0'
c
c         iii. set up format for pgmnam
c
      write(Flen,'(I4)') lnblnk(pgmnam)
      Fmt = '(A'//Flen//'$)'
c
c         iv. write out results
c
      write(*,Fmt) pgmnam
      if(onoff .eq. 1) then                      ! sign on
        write(*,301) vsn,cdate,ctime
      else                                       ! sign off
        dummyt = etime(second)
        write(*,302) vsn,cdate,ctime,second
      endif
  301 format(' version: ',a11,' - execution begun on ',a8,' at ',a8)
  302 format(' version: ',a11,' - execution ended on ',a8,' at ',a8
     *    /1x,f9.2,' cpu seconds used;',f8.2,' system seconds used.')
c
      return
c
      entry SignOff(pgmnam)
      OnOff = 2
      go to 100
c
      end
c
c=======================================================================
c
      Subroutine NextNarg(NArg,NArgs)
c
      integer NArg, NArgs
c
c-----------------------------------------------------------------------
c
      if (NArg .lt. NArgs) then
        NArg = NArg + 1
        return
      else
        print *,'ERROR: expected another argument but none found'
        call exit(64)
      end if
      return
      end
c
c=======================================================================
c
      subroutine upcase(string)
      character*(*) string
      integer*4 j, lnblnk
c
      do 10 j = 1,lnblnk(string)
         if(string(j:j) .ge. "a" .and. string(j:j) .le. "z") then
            string(j:j) = achar(iachar(string(j:j)) - 32)
         end if
10    continue
      return
      end
c      
c=======================================================================
c
      subroutine Cel2Ec(RA, Dec, Long, Lat)
c
      real*4 RA, Dec, Long, Lat, SOb, Cob, X, Y, Z, d2r,
     +       cRA, cDec, sRA, sDec, X2, Y2
c
c   Obliquity(2015) in J2000: 23.43734105     
c
      data d2r/1.745329252e-2/, cOb, sOb/0.9174956, 0.39777459/
c
c-----------------------------------------------------------------------
c
      cRA   = cos(d2r*RA)
      cDec  = cos(d2r*Dec)
      sRA   = sin(d2r*RA)
      sDec  = sin(d2r*Dec)
c
      X =  sDec
      Y = -cDec*sRA
      Z =  cDec*cRA
c      
      X2 =  X*Cob + Y*Sob
      Y2 = -X*Sob + Y*Cob
c     Z2 =  Z
c
      if (X2 .gt.  1.0) X2 =  1.0
      if (X2 .lt. -1.0) X2 = -1.0
c
      Lat  = asin(X2)/d2r
      Long = atan2(-Y2,Z)/d2r
      if (Long .lt. 0.0) Long = Long + 360.0
c
      return
      end
c      
c=======================================================================
c
      subroutine Cel2pa(RA, Dec, pa)
c
      real*4 RA, Dec, pa, SOb, Cob, d2r, r2d
c
c   Obliquity(2015) in J2000: 23.43734105     
c
      data d2r/1.745329252e-2/, cOb, sOb/0.9174956, -0.39777459/,
     +     r2d/57.29578/ 
c
c-----------------------------------------------------------------------
c
      pa = 90.0 - r2d*atan2(-sOb*sin(d2r*Dec)*sin(d2r*RA)
     +                      +cOb*cos(d2r*Dec), sOb*cos(d2r*RA))
      if (pa .gt. 360.0) pa = pa - 360.0
c     
      return
      end
c      
c=======================================================================
c          
      subroutine Rescale(FNam,Fac)
c
      Character*500 FNam, OutFNam
      Character*80  hdrline(100), record
      Real*4        Fac, Std(2048,2048), nullval
      Integer*4     FileID, NCols, NRows, NPlanes, nkey, nkeys, nsp,
     +              status, i, NPix, naxes(2), Access
      Logical*4     anynull
c
c-----------------------------------------------------------------------
c
      OutFNam = FNam
      FNam = FNam(1:lnblnk(FNam)-4)//'old.fits'
      if (Access(FNam(1:LNBlnk(FNam)),' ') .eq. 0) then
        print *
        print *,'WARNING: File already exists: ', FNam(1:LNBlnk(FNam))
        print *,'         STD rescaling skipped'
        print *
        return
        end if
c
      Call GetNAX(FNam,NCols,NRows,NPlanes,1,FileID)
      if ((NCols .ne. 2048) .or. (NRows .ne. 2048)
     +  .or. (NPlanes .ne. 1)) then
        print *,'ERROR: wrong NAXIS info for std:'
        print *,'       ',FNam(1:lnblnk(FNam))
        print *,'        file has NCols, NRows, NPlanes:',
     +                   NCols, NRows, NPlanes
        print *,'        rescaling will be aborted'
        return
      end if
c
      call ftghsp(FileID,nkeys,nsp,status)
      nkey = 0
      do 10 i=1,nkeys
          status = 0
          call ftgrec(FileID,i,record,status)
          if (record(1:6) .eq. 'SIMPLE')         go to 10
          if (record(1:6) .eq. 'BITPIX')         go to 10
          if (record(1:5) .eq. 'NAXIS')          go to 10
          if (record(17:24) .eq. 'Flexible')     go to 10
          if (record(15:26) .eq. 'Astrophysics') go to 10
          nkey = nkey + 1
          hdrline(nkey) = record
10    continue
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,Std,anynull,status)
      if (status .ne. 0) then
        print *,'ERROR reading  std file; status = ',status
        print *,'      rescaling will be aborted'
        return
      end if
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
      i = system('mv '//OutFNam(1:lnblnk(OutFNam))//' '//FNam)
      print *,'mv '//OutFNam(1:lnblnk(OutFNam))//' '//FNam
      Std = Fac*Std
      naxes(1) = 2048
      naxes(2) = 2048
      call wrfitsr(2,naxes,NPix,Std,OutFNam,status,nkey,hdrline)
c     
      return
      end
c
c=======================================================================
c
      Subroutine GetNAX(FilNam,NAXIS1,NAXIS2,NAXIS3,IDOp,FileID)
c-----------------------------------------------------------------------
c
c    Gets the dimensions of the file named FilNam; if IDOp = 0, the
c    file is closed, otherwise it's left open with handle FileID
c
c-----------------------------------------------------------------------
c
      Character*200 FilNam
      Integer*4     NAXIS, NAXIS1, NAXIS2, NAXIS3, IStat, ImOpen,
     +              ImRKeyI, FileID, IDOp, LNBlnk, ImClose
      integer status,readwrite,blocksize,naxes(3)
c
c-----------------------------------------------------------------------
c
C  The STATUS parameter must always be initialized.
      status=0
C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(FileID,status)
c     print *,'GetNAX: FileID, status:',fileid,status ! dbg
C  Open the FITS file 
      readwrite=0
      call ftopen(FileID,FilNam,readwrite,blocksize,status)
      if (status /= 0) then
          write(6,'(a)') 'GetNAX: Could not read '//trim(FilNam)
          istat = 3
          return
      endif
c
C  Determine the size of the image.
      call ftgknj(FileID,'NAXIS',1,2,naxes,NAXIS,status)
c     print *,'GetNAX: naxes, naxis, status:',naxes, naxis, status ! dbg
c
C  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (NAXIS .lt. 2)then
          print *,'GetNAX: Failed to read the NAXISn keywords.'
        istat = 4
          return
       end if
c
      NAXIS1 = naxes(1)
      NAXIS2 = naxes(2)
c
      If (NAXIS .gt. 2) then
        NAXIS3 = naxes(3)
      Else
        NAXIS3 = 1
      End If
c
      If (IDOp .eq. 0) then
        call ftclos(FileID, status)
        call ftfiou(FileID, status)
      end if
c
      return
c
      end
c      
c=======================================================================
c
      subroutine wrfitsr(naxis,naxes,lsize,array,fout,status,
     +                   nhdrline,hdrline)
c-----------------------------------------------------------------------
c
c  General-purpose real*4 FITS output subroutine
c
c    NOTE (JWF B60826): this version of FITSIO *MUST* have the first few
c                       header parameters in a certain order or else it 
c                       will fail on a read attempt; the following has worked:
c
c SIMPLE  =                    T / Written by IDL:  Fri Mar 25 12:38:25 2011      
c BITPIX  =                  -32 / Number of bits per data pixel                  
c NAXIS   =                    2 / Number of data axes                            
c NAXIS1  =                  641 /Number of positions along axis 1                
c NAXIS2  =                  641 /Number of positions along axis 2                
c  ...
c END                                                                             
c
c-----------------------------------------------------------------------
c
      integer*4     status,blocksize,outunit,naxis,naxes(naxis),bitpix,
     +              fpixel, group, nelements, nhdrline, n, lsize 
      real(4)       array(lsize)
      character*(*) fout
      character(80) hdrline(nhdrline)
      logical*4     simple, extend
c
      data          blocksize/1/, bitpix/-32/, simple/.true./,
     +              extend/.false./, fpixel/1/, group/1/
c
c-----------------------------------------------------------------------
c
c  The STATUS parameter must be initialized before using FITSIO.  A
c  positive value of STATUS is returned whenever a serious error occurs.
c  FITSIO uses an `inherited status' convention, which means that if a
c  subroutine is called with a positive input value of STATUS, then the
c  subroutine will exit immediately, preserving the status value. For 
c  simplicity, this program only checks the status value at the end of 
c  the program, but it is usually better practice to check the status 
c  value more frequently.
c     print *,'wrfitsr: file name = ',fout(1:lnblnk(fout))    ! dbg
c     do 10 n = 1, nhdrline                                   ! dbg
c       print *,'wrfitsr: header line ',n                     ! dbg
c       print *,hdrline(n)(1:lnblnk(hdrline(n)))              ! dbg
c10   continue                                                ! dbg
c
C  Get  unused Logical Unit Numbers to use to open the FITS files.
      status = 0 
      call ftgiou(outunit,status)
      if (status .ne. 0) return
c
      call ftinit(outunit,fout,blocksize,status)
      if (status .ne. 0) return
c
C  Initialize parameters about the FITS image.
C  BITPIX = 16 means that the image pixels will consist of 16-bit
C  integers.  The size of the image is given by the NAXES values. 
C  The EXTEND = TRUE parameter indicates that the FITS file
C  may contain extensions following the primary array.
c
C  Write the required header keywords to the file
      call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)
      if (status .ne. 0) return
c
      if (nhdrline .gt. 0) then
        do 100 n = 1, nhdrline
          status = 0
          call ftprec(outunit,hdrline(n),status)
          if (status .ne. 0) then
            print *, 'wrfitsr WARNING - failed to place header line:'
            print *,hdrline(n)(1:lnblnk(hdrline(n)))
          end if
100     continue
      end if      
c
C  Write the array to the FITS file.
C  The last letter of the subroutine name defines the datatype of the
C  array argument; in this case the 'J' indicates that the array has an
C  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
C  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
C  total number of pixels.  GROUP is seldom used parameter that should
C  almost always be set = 1.
      nelements = naxes(1)*naxes(2)
      if (naxis .eq. 3) nelements = nelements*naxes(3)
      status = 0
      call ftppre(outunit,group,fpixel,nelements,array,status)
      if (status .ne. 0) return
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(outunit, status)
      call ftfiou(outunit, status)
c
      return
      end
c
c=======================================================================
c
      subroutine ChkPAhist(PAhist,PAmin,PAmax,dbg)
c
      real*4    PAmin, PAmax
      integer*4 PAhist(3600), k, kLo, kHi, nEmptyOrig, nEmptyOcc, kk,
     +          nOcc, kDelt, kEmpty1, kDeltMax, kDelt1
      logical*4 ZeroX, dbg
c
c-----------------------------------------------------------------------
c
      nOcc  = 0                        ! count no. of occupied cells
      do 10 k = 1, 3600
        if (PAhist(k) .gt. 0) nOcc = nOcc + 1
10    continue
c
      if (nOcc .eq. 0) then            ! histogram is empty!
        PAmin = 0.0
        PAmax = 0.0
        return
      end if
c
      if (nOcc .gt. 3500) then         ! histogram is (almost) full!
        PAmin = 0.0
        PAmax = 360.0
        return
      end if
c                                      ! chk whether ~0 is populated
      kLo = 1
      do 20 k = 1, 3600
        if (PAhist(k) .gt. 0) then
          kLo = k
          go to 30
        end if
20    continue
c
30    kHi = 3600
      do 40 kk = 1, 3600               ! chk whether ~360 is populated
        k = 3601 - kk
        if (PAhist(k) .gt. 0) then
          kHi = k
          go to 50
        end if
40    continue
c
50    ZeroX = (kLo .lt. 25) .and. (kHi .gt. 3575)
      if (dbg) print *,'kLo, kHi:', kLo, kHi
      if (ZeroX) go to 500 
      nEmptyOrig = 3600 - kHi + kLo    ! no zero crossing; if the empty
      nEmptyOcc  = kHi - kLo - nOcc    ! region around zero is bigger
      if (dbg) then                    ! than in the occupied fraction,
        print *,                       ! return limits of latter
     +   'no zero crossing; nEmptyOrig, nEmptyOcc:',
     +    nEmptyOrig, nEmptyOcc
      end if
c
      if (nEmptyOrig .ge. nEmptyOcc) then
        PAmin = (float(kLo)-0.5)/10.0  ! return limits of latter
        PAmax = float(kHi)/10.0        
        return                         
      end if
c
60    kDelt    = 0                     ! find stretch of zeros not 
      kEmpty1  = -9                    ! around origin and bigger yet
      kDeltMax = -9
      kDelt1   = 0
      do 70 k = kLo, kHi
        if (PAhist(k) .eq. 0) then
          kDelt = kDelt + 1
          if (kDelt .eq. 1) kEmpty1 = k
        else
          if (kDelt .gt. kDeltMax) then
            kDeltMax = kDelt
            kDelt1  = kEmpty1
          end if
          kDelt = 0
        end if
70    continue
c
      if (dbg) print *,'kDelt1, kDeltMax:', kDelt1, kDeltMax
      if (kDeltMax .gt. nEmptyOrig) then  ! embedded empty region 
        PAmax = (kDelt1 - 1)/10.0         ! is bigger
        PAmin = (kDelt1 + kDeltMax)/10.0 - 360.0
        return
      else
        PAmin = (float(kLo)-0.5)/10.0  ! biggest embedded empty region
        PAmax = float(kHi)/10.0        ! is smaller than around origin
        return                              
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c                                      ! we have a zero crossing; this
500     kLo = 1                        ! boils down to previous case with
        kHi = 3600                     ! zeros around origin but no zeros
        nEmptyOrig = 0
        if (dbg) print *,'Zero crossing detected'
        go to 60
c
      return
      end
c
c=======================================================================
c
      subroutine GetTmpNam(TmpNam, TempDir, OK, dbg)
c      
      character*500 TmpNam0, TmpNam, WrkStrg, TempDir
      Character*11  Vsn
      Character*8   CDate, CTime
      integer*4     nTries, idum, k, lnblnk, Access, nCall
      real*4        ran1
      logical*4     OK, dbg
      data          nCall/0/
c      
      Common / VDT / CDate, CTime, Vsn
c
      if (nCall .eq. 0) then
        idum = -1
        nTries  = ran1(idum)        ! initialize random number generator
        nCall = nCall + 1
      end if
      TmpNam0 = TmpNam
      OK      = .true.
      nTries  = 0
c
10    nTries = nTries + 1
      if (nTries .gt. 999) then
        OK = .false.
        return
      end if
      k = 1.e6*ran1(idum)
      write (WrkStrg,'(i8)') k
20    k = index(WrkStrg,' ')
      if (k .lt. lnblnk(WrkStrg)) then 
        WrkStrg(k:k) = '%'
        go to 20
      end if      
      WrkStrg = WrkStrg(1:lnblnk(WrkStrg))//CDate//CTime
      k = index(WrkStrg,':')
      WrkStrg(k:k) = '%'
      k = index(WrkStrg,':')
      WrkStrg(k:k) = '%'
30    k = index(WrkStrg,' ')
      if (k .lt. lnblnk(WrkStrg)) then 
        WrkStrg(k:k) = '%'
        go to 30
      end if      
      TmpNam = TempDir(1:lnblnk(TempDir))//TmpNam0(1:lnblnk(TmpNam0))
     +         //WrkStrg(1:lnblnk(WrkStrg))
      if (Access(TmpNam(1:lnblnk(TmpNam)),' ') .eq. 0) go to 10
      if (dbg) print *,'returning file name: |'
     +                //TmpNam(1:lnblnk(TmpNam))//'|'
      return
      end
c
c=======================================================================
c
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*4 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
