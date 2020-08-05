!move_nodes.f90

!Algoritmo para deteccao de comunidades em rede, baseado no artigo:
!------------------------------------------------------------------
!Paper: Fast unfolding of communities in large networks
!Authors: Vincent D Blondel, Jean-Loup Guillaume, Renaud Lambiotte and Etienne Lefebvre
!doi:10.1088/1742-5468/2008/10/P10008
!------------------------------------------------------------------
!Este programa eh composto de tres partes: louvain.f90; move_nodes.f90; aggregate.f90

!A variavel ki(i) nesta subrotina eh diferente da usada em louvain.f90
      !variaveis comuns somente a move_nodes.f90 e suas funcoes
      module variaveis_move   
      integer*8, allocatable ::  ki(:)
      end module
 
      subroutine move_nodes() 
      use variaveis      
      use variaveis_move

      real*8 gain_q, q, best_q, best_c
      integer*8 kc, e_ic
      integer*8 i, j, move, n


!-------grau dos nos (ki) do grafo G1---
      allocate (ki(numb_comun))
      do i=1,numb_comun
            ki(i) = 0
         do j =1,numb_comun
            ki(i) = ki(i)+G1(i,j)
         end do
      end do 
!---------------------------------------
      
      deallocate(C1, C2)
      allocate (C1(numb_comun),C2(numb_comun)) 
      do i=1,numb_comun !Aloca o no i na comunidade Ci
         C2(i)=i 
         C1(i)=C2(i) 
      end do 

      move = numb_comun
      do while (move.ne.0)!mova enquanto houver ganho de modularidade
      move = numb_comun 
        do i=1,numb_comun
!          i= int(rand(0)*(numb_comun+1-1))+1
           best_q = -1000
           best_c = C2(i)
           do j=1,numb_comun
             if(G1(i,j).gt.0) then
             gain_q = (1.d0/m)*(e_ic(C2(j),i)- &
                      (ki(i)*kc(C2(j),i)/(2.d0*m)))
               if(best_q.lt.gain_q) then
                  best_q = gain_q
                  best_c = C2(j)
               end if
              end if
           end do
           C2(i) = best_c
           if(C1(i).eq.C2(i)) then
             move = move - 1
           end if
           C1(i) = C2(i)
        end do
      end do
      
      deallocate (ki)
      end subroutine move_nodes


!-------grau da comunidade (kc)--------------
         function kc(cmt,node)
         use variaveis 
         use variaveis_move
         integer*8 i, kc,node, cmt
         kc = 0
         do i = 1,numb_comun
              if (C2(i).eq.cmt) then 
                kc = kc+ki(i)
              end if
         end do
        if (C2(node).eq.cmt) then
            kc = kc - ki(node)
         end if
         return
         end function kc 

!-soma dos pesos dos links entre o no i e a comunidade C(cmt)-
         function e_ic(cmt,node)
         use variaveis         
         use variaveis_move
         integer*8 j, e_ic, cmt, node
         e_ic = 0
         do j = 1,numb_comun
         if(node.ne.j) then
           if ( (C2(j).eq.cmt).and.(G1(node,j).gt.0) ) then 
             e_ic = e_ic+G1(node,j)
           end if
         end if
         end do
         return
         end function e_ic 

