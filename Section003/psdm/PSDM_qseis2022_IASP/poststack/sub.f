c       ***************************
c       get string's actually length
c       ***************************
        Integer Function Lengths(str)
        Character*(*) str
**        write(istderr,*) len(str)
*        write(*,*) len(str)
        do 15 i = len(str),1,-1
                if(str(i:i) .NE. ' ') go to 20
15      continue
20      Lengths = i                 ! modified
        End


