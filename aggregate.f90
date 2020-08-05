!aggregate.f90

!Algoritmo para deteccao de comunidades em rede, baseado no artigo:
!------------------------------------------------------------------
!Paper: Fast unfolding of communities in large networks
!Authors: Vincent D Blondel, Jean-Loup Guillaume, Renaud Lambiotte and Etienne Lefebvre
!doi:10.1088/1742-5468/2008/10/P10008
!------------------------------------------------------------------
!Este programa eh composto de tres partes: louvain.f90; move_nodes.f90; aggregate.f90

      subroutine aggregate()
      use variaveis
      integer*8 i, j, n
      
      n=numb_comun    !Este n nao eh o mesmo do louvain.f90

!----Renumerar comunidades-----------------
      C1 = -1
      numb_comun = 1    
      do i=1, n
        where( (C2.eq.C2(i)).and.(C1.lt.0) )
          C1 = numb_comun  
        end where
        numb_comun = maxval(C1)+1    
      end do      

      
      numb_comun = maxval(C1)
      write(*,*)'Comunidade nova do no C(i):'
      do i=1,n
        write(*,*) 'C1(',i,')=', C1(i)
      end do

!---------Agregar comunidades-----------     
      allocate (G2(numb_comun, numb_comun))
      G2 = 0
      do i=1,n
        do j=1,n
              G2(C1(i),C1(j)) = G2(C1(i),C1(j)) + G1(i,j)
        end do
      end do
!------------------------------------------
      deallocate (G1)
      allocate (G1(numb_comun, numb_comun))
!------------Novo Grafo (G1)-----------
      G1 = G2 
      deallocate (G2)
      end subroutine aggregate
