!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr; 
!
!This software (Eco3M) is a computer program whose purpose is to perform 
!biogeochemical or coupled physical-biogeochemical modelling.
!
!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software. You can  use, 
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info". 
!
!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software''s author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability. 
!
!In this respect, the user''s attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software''s suitability as regards their
!requirements in conditions enabling the security of their systems and/or 
!data to be ensured and,  more generally, to use and operate it in the 
!same conditions as regards security. 
!
!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.
!***************************************************************************
!***************************************************************************
   MODULE eco3m_string
!
!> Module that contains the functions for string processing
!! It contains also the definition of the length of variable strings and
!! long strings.
!! \author Melika Baklouti, Vincent Faure
!! \date 2007-06-27
!----------------------------------------------------------------------------
 Implicit None
 Integer, Parameter :: L_CHAIN = 2500  !< Length of long string 
 Integer, Parameter :: L_VAR = 25      !< Length of strings for variables name
 Integer, Parameter :: L_VAR_LG = 40   !< Length of strings for long variables name 

CONTAINS

!-----------------------------------------------------------------------------------------------------------------
    character(L_CHAIN) function f_chain(chain,ipos,separateur)
    
!> This function returns the sub-string located between the ipos-1-th and
!! the ipos-th separator of a string
!-----------------------------------------------------------------------------------------------------------------

        ! arguments
        Character(len=L_CHAIN), intent(in) :: chain !< Input string
        Character(1), intent(in)   :: separateur !< Separator
        Integer, intent(in)        :: ipos !< Position of the separator

        character(len=L_CHAIN)  :: ch_tempo

        ! local variables:
        Integer::ii, ideb, ifin, ip

        ch_tempo = adjustl(chain)
        ch_tempo = trim(ch_tempo) // separateur

        ideb = 1
        ifin = len_trim(ch_tempo)
        ip = 0
        do
            ii=index(ch_tempo(ideb:ifin),separateur)+ideb-1
            ip = ip+1
            if (ip == ipos) then
                f_chain=ch_tempo(ideb:ii-1)
                f_chain=adjustl(f_chain)
                f_chain=trim(f_chain)
                exit
            endif
            ideb=ii+1
        enddo
        Return
    End function f_chain
!-----------------------------------------------------------------------------------------------------------------
    Integer function f_nschain(chain,separateur)

!> Returns the number of sub-string within an input string, separated by the input separator 
!-----------------------------------------------------------------------------------------------------------------

        Implicit None

        ! arguments
        Character(len=L_CHAIN), intent(in)  :: chain !< Input string
        Character(len=1), intent(in) :: separateur !< Separator

        ! local variables
        Integer::ii,ideb,ifin,ip
        character(len=L_CHAIN)  :: ch_tempo2

        if(len_trim(chain)/=0) then
            ch_tempo2 = adjustl(chain)
            ch_tempo2 = trim(ch_tempo2) // separateur

            ideb = 1
            ifin = len_trim(ch_tempo2)
            ip = 0 ! nb of occurences of the "separateur" string
            do
                ii=index(ch_tempo2(ideb:ifin),separateur)+ideb-1
                ip = ip+1
                ideb=ii+1
                if (ideb >= ifin) exit
            enddo

        else
            ip=0
        endif
        f_nschain = ip

        Return
    End function f_nschain

!-----------------------------------------------------------------------------------------------------------------
  
    Integer function f_nschain_seq(chain,sequence,lseq)
    
!> Returns the number of sub-string within an input string, separated by a string sequence of length lseq
!! (i.e same as f_nschain, but instead of using a char separator, it uses a string sequence) 

!-----------------------------------------------------------------------------------------------------------------
        Implicit None
        integer, intent(in) :: lseq !> Length of the separating string sequence
        Character(len=L_CHAIN), intent(in) :: chain !> String to process
        Character(len=lseq), intent(in) :: sequence !> Separating string sequence

        ! local variables
        Integer::ii,ideb,ifin,ip
        character(len=L_CHAIN)  :: ch_tempo2

        if(len_trim(chain)/=0) then
            ch_tempo2 = adjustl(chain)
            ch_tempo2 = trim(ch_tempo2) // sequence

            ideb = 1
            ifin = len_trim(ch_tempo2)
            ip = 0 ! nb of occurrence of the ":" string
            do
                ii=index(ch_tempo2(ideb:ifin),sequence)+ideb-1
                ip = ip+1
                ideb=ii+lseq
                if (ideb >= ifin) exit
            enddo

        else
            ip=0
        endif
        f_nschain_seq = ip

        Return
    End function f_nschain_seq

!-----------------------------------------------------------------------------------------------------------------
!
    character (L_CHAIN)  Function f_Int2chain (int,longchain)
!
!> Converts an integer into a left-ragged string, whose length is also returned
! \author M Baklouti, V Faure
!-----------------------------------------------------------------------------------------------------------------
    Implicit None
    Integer, intent(in) :: int !< Input integer that will be converted into strings
    Integer, intent(inout) :: longchain !< Input integer that will contain the length of the trimmed output string

  ! local variables
    character(L_CHAIN) :: chain

    write(chain,*) int
    chain = trim(adjustl(chain))
    longchain = len_trim(chain)
    f_Int2chain = chain(1:longchain)

    End Function f_Int2chain
!-----------------------------------------------------------------------------------------------------------------
    character (L_CHAIN)  Function f_Int2chain_file(int,longchain)

    !> Converts an integer into a left-ragged string, whose length is also returned. It is zero-padded.
    !! For instance, f_Int2chain_file(4, lenout) will return the string "00004", with "lenout=5" 
    !! \author Nicolas Barrier
!-----------------------------------------------------------------------------------------------------------------

    Implicit None
    Integer, intent(in) :: int !< Input integer that will be converted into strings
    Integer, intent(inout) :: longchain  !< Input integer that will contain the length of the trimmed output string

    character(L_CHAIN) :: chain

    write(chain,"(I5.5)") int
    chain = trim(adjustl(chain))
    longchain = len_trim(chain)
    f_Int2chain_file = chain(1:longchain)

    End Function f_Int2chain_file
!-----------------------------------------------------------------------------------------------------------------

END MODULE eco3m_string
