! louvain.f90

!Algoritmo para deteccao de comunidades em rede, baseado no artigo:
!------------------------------------------------------------------
!Paper: Fast unfolding of communities in large networks
!Authors: Vincent D Blondel, Jean-Loup Guillaume, Renaud Lambiotte and Etienne Lefebvre
!doi:10.1088/1742-5468/2008/10/P10008
!------------------------------------------------------------------
!Codigo escrito por: Daniel Prado (pradofisica@gmail.com) e Derick
!Este programa eh composto de tres partes: louvain.f90; move_nodes.f90; aggregate.f90

      !variaveis comuns as 3 partes do programa
      module variaveis     
      integer*8, allocatable :: C1(:), C2(:)
      integer*8, allocatable :: G1(:,:), G2(:,:)
      integer*8  m, numb_comun 
      end module variaveis

      !variaveis comuns somente a louvain.f90 e suas funcoes
      module variaveis_louvain 
      integer*8 n      
      integer*8, allocatable :: v(:,:), C0(:)
      integer*8, allocatable :: ki(:)
      end module variaveis_louvain

      program louvain
      use variaveis
      use variaveis_louvain
      implicit none
      real*8 Q1, Q2, modularity
      integer*8 n1
      integer*8 i, j
      integer*8, allocatable ::  lar(:)
      character(len=80) arq_graph, file_out 

      open(50, file = 'entra.dat') 
       read(50,*) arq_graph   !Nome do arquivo com a matriz de adjacencia 
       read(50,*) n           !numero de nos do grafo
      close(50)

      allocate (v(n,n), lar(n))

!--------leitura do grafo----------
!      open(100, file = "grafo5space.dat")
       open(100, file = arq_graph)
       
        do i=1,n
           read(100,500) (lar(j),j=1,n)
           do j=1,n      
           v(i,j) = lar(j) 
           end do         
        end do 
         
!---------escrita da matriz de adjacencia na tela----------
        write(*,*) 'Grafo:', arq_graph
        do i=1,n
           write(*,500) (v(i,j), j=1,n)
        end do 

      close(100)
!---------verifica a simetria na matriz de adjacencia---------------------
      do i=1,n
        do j=1,n
         if(v(i,j).ne.v(j,i)) then
           write(*,*) 'sem simetria','i=', i,'j=',j
         end if
        end do
      end do
!---------------------------------------           
500   format (10000I2)
!-------grau dos nos (ki) do grafo v--------
      allocate (ki(n))
      write(*,*)'grau do no i:'
      do i=1,n
            ki(i) = 0
         do j =1,n
            ki(i) = ki(i)+v(i,j) 
         end do
      end do 
!------- soma das arestas--------------------
      m = sum(v)/2
      write(*,*)'Numero de arestas(m) =', m  

!************Louvain*************************
      allocate (C0(n), C1(n), C2(n), G1(n,n))

      do i=1,n !comunidade inicial do no i
         C0(i)=i  
      end do 
  
      numb_comun = n
      G1 = v          
      Q1 = -1000.0
      Q2 = modularity()

      do while(Q1.lt.Q2)
         Q1=Q2
         call move_nodes()       !primero passo do metodo Louvain
         call aggregate()        !segundo passo do metodo Louvain
      
          !coloca o no do grafo original na nova comunidade
          n1 = maxval(C0)   
          do i=1, n1
             where(C0.eq.i)
               C0 = C1(i)
             end where
          end do
      
          !calcula modularidade 
          Q2 = modularity()

      end do

!****************Escrever dados de saida****************      
      file_out =  'modularity_'//arq_graph
      open(51,file= file_out)    
      !Escreve os dados de saida 
      write(51,*) arq_graph
      
      write(51,*) 'Matriz de adjacencia entre as comunidades:'
      do i=1,numb_comun
        write(51,*) (G1(i,j), j=1,numb_comun)
      end do 
      write(51,fmt='(/A)') 'Comunidades:'
      do i=1, numb_comun        
        write(51,fmt='(A,I4,A)', advance='no') 'C', i, '['
        do j=1,n
           if(C0(j).eq.i) then
             write(51,fmt='(I5,A1)', advance='no') j,',' 
           end if
        end do
        write(51,fmt='(A)')']'
      end do
      write(51,fmt='(/A,I6)')'numero de comunidades = ', numb_comun
      write(51,*)'modularidade = ', Q2

      close(51)
      write(*,*)'O arquivo de saida ja esta na pasta ', file_out
      end program louvain 

!********************************************************
!--------Calcula modularidade-----------------
      function modularity()
      use variaveis
      use variaveis_louvain
      real*8 modularity
      integer*8 i, j

      modularity = 0.d0
      do i=1,n
         do j=1,n
           if(C0(i).eq.C0(j)) then
             modularity = modularity+v(i,j)-ki(i)*ki(j)/(2.0*m)
           end if 
         end do
      end do

      modularity = modularity/(2*m)
      return
      end function modularity

!------------sub-rotinas--------------------
      include 'move_nodes.f90'
      include 'aggregate.f90'
