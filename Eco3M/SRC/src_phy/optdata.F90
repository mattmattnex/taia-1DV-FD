!
! ======================================================================
!
       subroutine OPTDATA
       Use comdynmod	
       Use comrunmod
       Use mod_varphy_coupl
!
       integer :: err
!
!       open(93,file=trim(cdir)//'/'//coptique,status='old')
       open(93,file=trim(cdir)//'/inputs/'//coptique,iostat=err)
!
      
       do jm=1,19
         read(93,*) (SRDIR(jl,jm),jl=1,jpzt)
       enddo
       do jm=1,19
         read(93,*) (SRDIF(jl,jm),jl=1,jpzt)
       enddo
       do jw=1,4
         read(93,*) (RDIR(jm,jw),jm=1,19)
       enddo
       read(93,*) rdif
       read(93,*) (XKW(jl),jl=1,jpzt)
       read(93,*) (CHI(jl),jl=1,jpzt)
       read(93,*) (EMO(jl),jl=1,jpzt)
       read(93,*) (APHY(jl),jl=1,jpzt)
       read(93,*) aphymax
       read(93,*) (AW(jl),jl=1,jpzt)
       read(93,*) (ASTAR(jl),jl=1,jpzt)
       read(93,*) (BW(jl),jl=1,jpzt)
       read(93,*) bstar  
       read(93,*) cstar  
       xdwl=5.*nfrog

       close(93)

       return
       end
