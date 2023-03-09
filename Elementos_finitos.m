%%Programa que resuelve desplazamiento axial en vigas por el método de
%%elementos finitos
%%Este codigo debe ser capaz de realizar el metodo por interpolacion lineal
%%o cuadratica, solicitar la longitud de los elementos, si son uniformes o
%%no y generar las matrices de inercia y rigidez, y encontrar la matriz de
%%rigidez global.

%%Propiedades de la barra en cuestion
L_1 = 4/12; %%Longitud de primera barra
    d_0 = 28; %%diametro inicial de L1
    d_1 = 20; %%diametro final de L1
L_2 = 6/12; %%longitud de segunda barra
L_3 = 2/12; %%longitud de tercera barra
    d_3 = 20; %%diametro inicial de L3
    d_4 = 16; %%diametro final de L3
L_4 = 3/12; %%longitud de cuarta barra
E = 10152642; %Modulo de elasticidad del aluminio en psi
rho = 0.0975437; %%densidad del aluminio en lbm/in^3
syms x %Define a x como una variable
fprintf('\n');
disp('METODO DE ELEMENTOS FINITOS PARA BARRAS')
fprintf('\n');
disp('Seleccione el tipo de interpolacion, presione 1 para lineal, 2 para cuadratica')
fprintf('\n');
n=input('Mi tipo de interpolacion='); %%opcion de tipo de interpolacion

if n == 1
    fprintf('\n');
    disp('Seleccione 1 para elementos de longitud constante, 2 para longitud variable')
    o = input('tipo de elemento finito'); %%seleción de longitud de elemento finito
    if o == 1
        disp('Cuantos elementos finitos desea por seccion')
        fprintf('\n');
        k=input('cantidad de elementos: '); %%cantidad de elementos finitos por barra
        dis_elem_1 = zeros(1,k); %vectores que guardan la longitud de los vectores de elementos finitos
        dis_elem_2 = zeros(1,k); % estos solo crean los vectores
        dis_elem_3 = zeros(1,k);
        dis_elem_4 = zeros(1,k);
        for i=1:k
            dis_elem_1(i) = L_1/k; %son llenados con el tamaño de cada elemento
            dis_elem_2(i) = L_2/k;
            dis_elem_3(i) = L_3/k;
            dis_elem_4(i) = L_4/k;
        end
            lim_elem_1 = zeros(k,2);
            lim_elem_2 = zeros(k,2);
            lim_elem_3 = zeros(k,2);
            lim_elem_4 = zeros(k,2);
        for i=0:(k-1)
            for j=1:2
                y=j-1;
            lim_elem_1(i+1,j) = [((i+y)*dis_elem_1(1))]; %discretización por elemento finito sec-1
            lim_elem_2(i+1,j) = [((i+y)*dis_elem_2(1)+(4/12))]; %discretización por elemento finito sec-2 
            lim_elem_3(i+1,j) = [((i+y)*dis_elem_3(1)+(10/12))]; %discretización por elemento finito sec-3
            lim_elem_4(i+1,j) = [((i+y)*dis_elem_4(1)+(1))]; %discretización por elemento finito sec-4
            end
        end
            %Funciones de interpolación
            syms x %x1 limite inferior x2 limite superior
            x1=1;
            x2=2;
            N = [((x-x2)/(x1-x2)),((x-x1)/(x2-x1))];
            
            %calculo de las matrices de Rigidez por elemento finito
            A1 = pi*(-2*x+14)^(2);
            %Primera seccion
            K_1 = cell(k,1);
            kk = zeros(2,2);
            for v=1:k
                x1=lim_elem_1(v,1);
                x2=lim_elem_1(v,2);
                for i=1:2
                    for j=1:2
                        G =A1*diff(N(i))*diff(N(j));
                        kk(i,j) = E*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                    end
                end
                K_1{v,1}=kk;
            end
            %Segunda Seccion
            A2 = pi*(d_1/2)^2; %Area constante de la seccion 2
            K_2 = cell(k,1); %Matriz que guarda las matrices de rigidez de los elementos de seccion2
            for v=1:k
                x1=lim_elem_2(v,1); %Limites de integracion por elemento
                x2=lim_elem_2(v,2); %limites de integracion por elemento
                for i=1:2
                    for j=1:2
                        G =A2*diff(N(i))*diff(N(j));
                        kk(i,j) = E*int(G,[lim_elem_2(v,1) lim_elem_2(v,2)]); %asignacion de elementos mat. rigidez
                    end
                end
                K_2{v,1}=kk;
            end
            %Tercera seccion
            A3 = pi*(-x+20); %funcion de area para seccion 3
            K_3 = cell(k,1); %Matriz que guarda las matrices de rigidez de los elementos de seccion4
            for v=1:k
                x1=lim_elem_3(v,1); %Limites de integracion por elemento
                x2=lim_elem_3(v,2); %limites de integracion por elemento
                for i=1:2
                    for j=1:2
                        G =A3*diff(N(i))*diff(N(j));
                        kk(i,j) = E*int(G,[lim_elem_3(v,1) lim_elem_3(v,2)]); %asignacion de elementos mat. rigidez
                    end
                end
                K_3{v,1}=kk;
            end
            
            %Cuarta seccion
            A4 = pi*(8)^2; %Area de seecion 4
            K_4 = cell(k,1); %Matriz que guarda las matrices de rigidez de los elementos de seccion4
            for v=1:k
                x1=lim_elem_4(v,1); %Limites de integracion por elemento
                x2=lim_elem_4(v,2); %limites de integracion por elemento
                for i=1:2
                    for j=1:2
                        G =A4*diff(N(i))*diff(N(j));
                        kk(i,j) = E*int(G,[lim_elem_4(v,1) lim_elem_4(v,2)]); %asignacion de elementos mat. rigidez
                    end
                end
                K_4{v,1}=kk;
            end
            
                %ENSAMBLAJE DE LAS MATRICES DE RIGIDEZ, RIGIDEZ GLOBAL
                
                %Matriz Global de rigidez seccion 1
                Kglobal_1 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Kglobal_1(i+(c-1),j+(c-1))=K_1{c}(i,j)+Kglobal_1(i+(c-1),j+(c-1));
                      end
                    end
                end
                %Matriz Global de rigidez seccion 2
                Kglobal_2 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Kglobal_2(i+(c-1),j+(c-1))=K_2{c}(i,j)+Kglobal_2(i+(c-1),j+(c-1));
                      end
                    end
                end
                %Matriz Global de rigidez seccion 3
                Kglobal_3 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Kglobal_3(i+(c-1),j+(c-1))=K_3{c}(i,j)+Kglobal_3(i+(c-1),j+(c-1));
                      end
                    end
                end
                %MatrizGlobal de rigidez seccion 4
                Kglobal_4 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                        for j=1:2
                            Kglobal_4(i+(c-1),j+(c-1))=K_4{c}(i,j)+Kglobal_4(i+(c-1),j+(c-1));
                        end
                    end
                end
               %Matriz Global de rigidez total
               
               
               %ponemos todas las matrices gobales en un arreglo tipo
               %celula
               CelK = cell(4,1);
               CelK{1} = Kglobal_1;
               CelK{2} = Kglobal_2;
               CelK{3} = Kglobal_3;
               CelK{4} = Kglobal_4;
               Kglobaltotal = zeros( (4*(k-1)+5),(4*(k-1)+5) );
               for c=1:4
                    for i= (1:(k+1))
                        for j= (1:(k+1))
                            Kglobaltotal(i+(k)*(c-1),j+(k)*(c-1))= CelK{c}(i,j)+Kglobaltotal(i+(k)*(c-1),j+(k)*(c-1));
                        end
                    end
               end
               
            %CALCULO DE LAS MATRICES DE INERCIA POR ELEMENTO FINITO
             A1 = pi*(-2*x+14)^(2);
             A2 = pi*(d_1/2)^2;
             A3 = pi*(-x+10);
             A4 = pi*(8)^2;
            %Matriz de inercia seccion 1
            M_1 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(2,2); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_1(v,1);
                x2=lim_elem_1(v,2);
                for i=1:2
                    for j=1:2
                        G =A1*(N(i))*(N(j));
                        MM(i,j) = rho*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                    end
                end
                M_1{v,1}=MM;
            end
            
            %Matriz de inercia seccion 2
            M_2 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(2,2); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_2(v,1);
                x2=lim_elem_2(v,2);
                for i=1:2
                    for j=1:2
                        G =A2*(N(i))*(N(j));
                        MM(i,j) = rho*int(G,[lim_elem_2(v,1) lim_elem_2(v,2)]);
                    end
                end
                M_2{v,1}=MM;
            end
            
            %Matriz de inercia seccion 3
            M_3 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(2,2); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_3(v,1);
                x2=lim_elem_3(v,2);
                for i=1:2
                    for j=1:2
                        G =A3*(N(i))*(N(j));
                        MM(i,j) = rho*int(G,[lim_elem_3(v,1) lim_elem_3(v,2)]);
                    end
                end
                M_3{v,1}=MM;
            end
               
            %Matriz de inercia seccion 4  
            M_4 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(2,2); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_4(v,1);
                x2=lim_elem_4(v,2);
                for i=1:2
                    for j=1:2
                        G =A4*(N(i))*(N(j));
                        MM(i,j) = rho*int(G,[lim_elem_4(v,1) lim_elem_4(v,2)]);
                    end
                end
                M_4{v,1}=MM;
            end
           
            %Matriz de emsamblaje seccion 1
            Mglobal_1 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Mglobal_1(i+(c-1),j+(c-1))=M_1{c}(i,j)+Mglobal_1(i+(c-1),j+(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 2
            Mglobal_2 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Mglobal_2(i+(c-1),j+(c-1))=M_2{c}(i,j)+Mglobal_2(i+(c-1),j+(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 3
            Mglobal_3 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Mglobal_3(i+(c-1),j+(c-1))=M_3{c}(i,j)+Mglobal_3(i+(c-1),j+(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 4
            Mglobal_4 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Mglobal_4(i+(c-1),j+(c-1))=M_4{c}(i,j)+Mglobal_4(i+(c-1),j+(c-1));
                      end
                    end
                end
                
           %Ensamblaje de matriz de inercia de toda la barra
           %ponemos todas las matrices gobales en un arreglo tipo
               %celula
               elK = cell(4,1);
               elK{1} = Mglobal_1;
               elK{2} = Mglobal_2;
               elK{3} = Mglobal_3;
               elK{4} = Mglobal_4;
               Mglobaltotal = zeros( (4*(k-1)+5),(4*(k-1)+5) );
               for c=1:4
                    for i= (1:(k+1))
                        for j= (1:(k+1))
                            Mglobaltotal(i+(k)*(c-1),j+(k)*(c-1))= elK{c}(i,j)+Mglobaltotal(i+(k)*(c-1),j+(k)*(c-1));
                        end
                    end
               end
               
    %Imposición de condiciones de frontera de la barra en cuestion
    Kglobaltotal_t = zeros(4*(k-1)+4,4*(k-1)+4);
    Mglobaltotal_t = zeros(4*(k-1)+4,4*(k-1)+4);
    for i=2:(4*(k-1)+5)
        for j=2:(4*(k-1)+5)
            Kglobaltotal_t(i-1,j-1)=Kglobaltotal(i,j);
            Mglobaltotal_t(i-1,j-1)=Mglobaltotal(i,j);
        end
    end
  
    disp('Matriz global de rigidez con condiciones de frontera')
    Kglobaltotal_t
    disp('Matriz global de inercia con condiciones de frontera')
    Mglobaltotal_t
    disp('Modos y Frecuencias naturales')
    
               %Calculo de modos y frecuencias naturales
               MI = inv (Mglobaltotal_t);
               KM = MI*Kglobaltotal_t;
               [L,V]= eig(KM) %Vector de modo L, Frecuencias naturales V
              
    elseif o == 2
        disp('Cuantos elementos finitos desea por seccion')
        fprintf('\n');
        k=input('cantidad de elementos'); %%cantidad de elementos finitos por barra
        dis_elem_1 = zeros(1,k); %vectores que guardan la longitud de los vectores de elementos finitos
        dis_elem_2 = zeros(1,k); % estos solo crean los vectores
        dis_elem_3 = zeros(1,k);
        dis_elem_4 = zeros(1,k);
        
        disp('Introducir longitudes de la primera seccion de la barra')
        fprintf('\n');
        for u=1:(k-1)
            y=input('Longitud de subsegmento n');
            dis_elem_1(u)=y %se deja ver para que el usuario pueda apreciar que valores va introduciendo
        end
        sum1 = 0; %suma para encontrar el valor del ultimo subsegmento
        for t=1:(k-1)
            sum1=sum1+dis_elem_1(t);
        end
        if sum1 >= L_1
            disp('Las medidas de sus elementos en la seccion 1 excede la longitud total del mismo, reinicie')
        end
        
        dis_elem_1(k)= L_1 - sum1;
        
        disp('Introducir longitudes de la segunda barra seccion de la barra')
        fprintf('\n');
        for u=1:(k-1)
            y=input('Longitud de subsegmento n');
            dis_elem_2(u)=y %se deja ver para que el usuario pueda apreciar que valores va introduciendo
        end
        sum2 = 0; %suma para encontrar el valor del ultimo subsegmento
        for t=1:(k-1)
            sum2=sum2+dis_elem_2(t);
        end
        if sum1 >= L_2
            disp('Las medidas de sus elementos en la seccion 2 excede la longitud total del mismo, reinicie')
        end
        
        dis_elem_1(k)= L_2 - sum2;
        
        disp('Introducir longitudes de la tercera seccion de la barra')
        fprintf('\n');
        for u=1:(k-1)
            y=input('Longitud de subsegmento n');
            dis_elem_3(u)=y %se deja ver para que el usuario pueda apreciar que valores va introduciendo
        end
        sum3 = 0; %suma para encontrar el valor del ultimo subsegmento
        for t=1:(k-1)
            sum3=sum3+dis_elem_3(t);
        end
        if sum3 >= L_3
            disp('Las medidas de sus elementos en la seccion 3 excede la longitud total del mismo, reinicie')
        end
        
        dis_elem_3(k)= L_3 - sum3;
        
        disp('Introducir longitudes de la cuarta seccion de la barra')
        fprintf('\n');
        for u=1:(k-1)
            y=input('Longitud de subsegmento n');
            dis_elem_4(u)=y %se deja ver para que el usuario pueda apreciar que valores va introduciendo
        end
        sum4 = 0; %suma para encontrar el valor del ultimo subsegmento
        for t=1:(k-1)
            sum1=sum1+dis_elem_4(t);
        end
        if sum4 >= L_4
            disp('Las medidas de sus elementos en la seccion 1 excede la longitud total del mismo, reinicie')
        end
        
        dis_elem_4(k)= L_4 - sum4;
      
        lim_elem_1 = zeros(k,2);
            lim_elem_2 = zeros(k,2);
            lim_elem_3 = zeros(k,2);
            lim_elem_4 = zeros(k,2);
        for i=0:(k-1)
            for j=1:2
                y=j-1;
            lim_elem_1(i+1,j) = [((i+y)*dis_elem_1(1))]; %discretización por elemento finito sec-1
            lim_elem_2(i+1,j) = [((i+y)*dis_elem_2(1)+(4/12))]; %discretización por elemento finito sec-2 
            lim_elem_3(i+1,j) = [((i+y)*dis_elem_3(1)+(10/12))]; %discretización por elemento finito sec-3
            lim_elem_4(i+1,j) = [((i+y)*dis_elem_4(1)+(1))]; %discretización por elemento finito sec-4
            end
        end
            %Funciones de interpolación
            syms x %x1 limite inferior x2 limite superior
            x1=1;
            x2=2;
            N = [((x-x2)/(x1-x2)),((x-x1)/(x2-x1))];
            
            %calculo de las matrices de Rigidez por elemento finito
            A1 = pi*(-2*x+14)^(2);
            %Primera seccion
            K_1 = cell(k,1);
            kk = zeros(2,2);
            for v=1:k
                x1=lim_elem_1(v,1);
                x2=lim_elem_1(v,2);
                for i=1:2
                    for j=1:2
                        G =A1*diff(N(i))*diff(N(j));
                        kk(i,j) = E*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                    end
                end
                K_1{v,1}=kk;
            end
            %Segunda Seccion
            A2 = pi*(d_1/2)^2; %Area constante de la seccion 2
            K_2 = cell(k,1); %Matriz que guarda las matrices de rigidez de los elementos de seccion2
            for v=1:k
                x1=lim_elem_2(v,1); %Limites de integracion por elemento
                x2=lim_elem_2(v,2); %limites de integracion por elemento
                for i=1:2
                    for j=1:2
                        G =A2*diff(N(i))*diff(N(j));
                        kk(i,j) = E*int(G,[lim_elem_2(v,1) lim_elem_2(v,2)]); %asignacion de elementos mat. rigidez
                    end
                end
                K_2{v,1}=kk;
            end
            %Tercera seccion
            A3 = pi*(-x+20); %funcion de area para seccion 3
            K_3 = cell(k,1); %Matriz que guarda las matrices de rigidez de los elementos de seccion4
            for v=1:k
                x1=lim_elem_3(v,1); %Limites de integracion por elemento
                x2=lim_elem_3(v,2); %limites de integracion por elemento
                for i=1:2
                    for j=1:2
                        G =A3*diff(N(i))*diff(N(j));
                        kk(i,j) = E*int(G,[lim_elem_3(v,1) lim_elem_3(v,2)]); %asignacion de elementos mat. rigidez
                    end
                end
                K_3{v,1}=kk;
            end
            
            %Cuarta seccion
            A4 = pi*(8)^2; %Area de seecion 4
            K_4 = cell(k,1); %Matriz que guarda las matrices de rigidez de los elementos de seccion4
            for v=1:k
                x1=lim_elem_4(v,1); %Limites de integracion por elemento
                x2=lim_elem_4(v,2); %limites de integracion por elemento
                for i=1:2
                    for j=1:2
                        G =A4*diff(N(i))*diff(N(j));
                        kk(i,j) = E*int(G,[lim_elem_4(v,1) lim_elem_4(v,2)]); %asignacion de elementos mat. rigidez
                    end
                end
                K_4{v,1}=kk;
            end
            
                %ENSAMBLAJE DE LAS MATRICES DE RIGIDEZ, RIGIDEZ GLOBAL
                
                %Matriz Global de rigidez seccion 1
                Kglobal_1 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Kglobal_1(i+(c-1),j+(c-1))=K_1{c}(i,j)+Kglobal_1(i+(c-1),j+(c-1));
                      end
                    end
                end
                %Matriz Global de rigidez seccion 2
                Kglobal_2 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Kglobal_2(i+(c-1),j+(c-1))=K_2{c}(i,j)+Kglobal_2(i+(c-1),j+(c-1));
                      end
                    end
                end
                %Matriz Global de rigidez seccion 3
                Kglobal_3 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Kglobal_3(i+(c-1),j+(c-1))=K_3{c}(i,j)+Kglobal_3(i+(c-1),j+(c-1));
                      end
                    end
                end
                %MatrizGlobal de rigidez seccion 4
                Kglobal_4 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                        for j=1:2
                            Kglobal_4(i+(c-1),j+(c-1))=K_4{c}(i,j)+Kglobal_4(i+(c-1),j+(c-1));
                        end
                    end
                end
               %Matriz Global de rigidez total
               
               
               %ponemos todas las matrices gobales en un arreglo tipo
               %celula
               CelK = cell(4,1);
               CelK{1} = Kglobal_1;
               CelK{2} = Kglobal_2;
               CelK{3} = Kglobal_3;
               CelK{4} = Kglobal_4;
               Kglobaltotal = zeros( (4*(k-1)+5),(4*(k-1)+5) );
               for c=1:4
                    for i= (1:(k+1))
                        for j= (1:(k+1))
                            Kglobaltotal(i+(k)*(c-1),j+(k)*(c-1))= CelK{c}(i,j)+Kglobaltotal(i+(k)*(c-1),j+(k)*(c-1));
                        end
                    end
               end
               
            %CALCULO DE LAS MATRICES DE INERCIA POR ELEMENTO FINITO
             A1 = pi*(-2*x+14)^(2);
             A2 = pi*(d_1/2)^2;
             A3 = pi*(-x+10);
             A4 = pi*(8)^2;
            %Matriz de inercia seccion 1
            M_1 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(2,2); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_1(v,1);
                x2=lim_elem_1(v,2);
                for i=1:2
                    for j=1:2
                        G =A1*(N(i))*(N(j));
                        MM(i,j) = rho*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                    end
                end
                M_1{v,1}=MM;
            end
            
            %Matriz de inercia seccion 2
            M_2 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(2,2); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_2(v,1);
                x2=lim_elem_2(v,2);
                for i=1:2
                    for j=1:2
                        G =A2*(N(i))*(N(j));
                        MM(i,j) = rho*int(G,[lim_elem_2(v,1) lim_elem_2(v,2)]);
                    end
                end
                M_2{v,1}=MM;
            end
            
            %Matriz de inercia seccion 3
            M_3 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(2,2); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_3(v,1);
                x2=lim_elem_3(v,2);
                for i=1:2
                    for j=1:2
                        G =A3*(N(i))*(N(j));
                        MM(i,j) = rho*int(G,[lim_elem_3(v,1) lim_elem_3(v,2)]);
                    end
                end
                M_3{v,1}=MM;
            end
               
            %Matriz de inercia seccion 4  
            M_4 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(2,2); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_4(v,1);
                x2=lim_elem_4(v,2);
                for i=1:2
                    for j=1:2
                        G =A4*(N(i))*(N(j));
                        MM(i,j) = rho*int(G,[lim_elem_4(v,1) lim_elem_4(v,2)]);
                    end
                end
                M_4{v,1}=MM;
            end
           
            %Matriz de emsamblaje seccion 1
            Mglobal_1 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Mglobal_1(i+(c-1),j+(c-1))=M_1{c}(i,j)+Mglobal_1(i+(c-1),j+(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 2
            Mglobal_2 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Mglobal_2(i+(c-1),j+(c-1))=M_2{c}(i,j)+Mglobal_2(i+(c-1),j+(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 3
            Mglobal_3 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Mglobal_3(i+(c-1),j+(c-1))=M_3{c}(i,j)+Mglobal_3(i+(c-1),j+(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 4
            Mglobal_4 = zeros((k+1),(k+1));
                for c=1:k
                    for i=1:2
                      for j=1:2
                        Mglobal_4(i+(c-1),j+(c-1))=M_4{c}(i,j)+Mglobal_4(i+(c-1),j+(c-1));
                      end
                    end
                end
                
           %Ensamblaje de matriz de inercia de toda la barra
           %ponemos todas las matrices gobales en un arreglo tipo
               %celula
               elK = cell(4,1);
               elK{1} = Mglobal_1;
               elK{2} = Mglobal_2;
               elK{3} = Mglobal_3;
               elK{4} = Mglobal_4;
               Mglobaltotal = zeros( (4*(k-1)+5),(4*(k-1)+5) );
               for c=1:4
                    for i= (1:(k+1))
                        for j= (1:(k+1))
                            Mglobaltotal(i+(k)*(c-1),j+(k)*(c-1))= elK{c}(i,j)+Mglobaltotal(i+(k)*(c-1),j+(k)*(c-1));
                        end
                    end
               end
    
          %Imposición de condiciones de frontera de la barra en cuestion
    Kglobaltotal_t = zeros(4*(k-1)+4,4*(k-1)+4);
    Mglobaltotal_t = zeros(4*(k-1)+4,4*(k-1)+4);
    for i=2:(4*(k-1)+5)
        for j=2:(4*(k-1)+5)
            Kglobaltotal_t(i-1,j-1)=Kglobaltotal(i,j);
            Mglobaltotal_t(i-1,j-1)=Mglobaltotal(i,j);
        end
    end
 
    disp('Matriz global de rigidez con condiciones de frontera')
    Kglobaltotal_t
    disp('Matriz global de inercia con condiciones de frontera')
    Mglobaltotal_t
    disp('Modos y Frecuencias naturales')
    
               %Calculo de modos y frecuencias naturales
               MI = inv (Mglobaltotal);
               KM = MI*Kglobaltotal;
               [L,V]= eig(KM) %Vector de modo L, Frecuencias naturales V
              
    else 
         fprintf('\n');
         disp('**El valor introducido no es valido**')
    end
    
elseif n == 2
   fprintf('\n');
    disp('Seleccione 1 para elementos de longitud constante, 2 para longitud variable')
    w = input('tipo de elemento finito'); %%seleción de longitud de elemento finito
     if w == 1
        disp('Cuantos elementos finitos desea por seccion')
        fprintf('\n');
        k=input('cantidad de elementos'); %%cantidad de elementos finitos por barra
        
        Dis_elem_1 = zeros(1,k); %vectores que guardan la longitud de los vectores de elementos finitos
        Dis_elem_2 = zeros(1,k); % estos solo crean los vectores
        Dis_elem_3 = zeros(1,k);
        Dis_elem_4 = zeros(1,k);
        for i=1:k
            dis_elem_1(i) = L_1/k; %son llenados con el tamaño de cada elemento
            dis_elem_2(i) = L_2/k;
            dis_elem_3(i) = L_3/k;
            dis_elem_4(i) = L_4/k;
        end
            lim_elem_1 = zeros(k,2);
            lim_elem_2 = zeros(k,2);
            lim_elem_3 = zeros(k,2);
            lim_elem_4 = zeros(k,2);
        for i=0:(k-1)
            for j=1:2
                y=j-1;
            lim_elem_1(i+1,j) = [((i+y)*dis_elem_1(1))]; %discretización por elemento finito sec-1
            lim_elem_2(i+1,j) = [((i+y)*dis_elem_2(1)+(4/12))]; %discretización por elemento finito sec-2 
            lim_elem_3(i+1,j) = [((i+y)*dis_elem_3(1)+(10/12))]; %discretización por elemento finito sec-3
            lim_elem_4(i+1,j) = [((i+y)*dis_elem_4(1)+(1))]; %discretización por elemento finito sec-4
            end
        end
            %Funciones de interpolación
            syms x %x1 limite inferior x2 limite superior
            x1=1;
            x2=2;
            x3=((x1+x2)/2);
            N = [(((x-x2)*(x-x3))/((x1-x2)*(x1-x3))),(((x-x1)*(x-x3))/((x2-x1)*(x2-x3))),((x-x1)*(x-x2))/((x3-x1)*(x3-x2))];
            
            %calculo de las matrices de Rigidez por elemento finito
            A1 = pi*(-2*x+14)^(2);
            %Primera seccion
            K_1 = cell(k,1);
            kk = zeros(3,3);
            for v=1:k
                x1=lim_elem_1(v,1);
                x2=lim_elem_1(v,2);
                x3= ((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A1*diff(N(i))*diff(N(j))*diff(N(r));
                        kk(i,j) = E*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                        end
                    end
                end
                K_1{v,1}=kk;
            end
            
            %Segunda Seccion
            A2 = pi*(d_1/2)^2; %Area constante de la seccion 2
            K_2 = cell(k,1);
            kk = zeros(3,3);
            for v=1:k
                x1=lim_elem_2(v,1);
                x2=lim_elem_2(v,2);
                x3= ((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A2*diff(N(i))*diff(N(j))*diff(N(r));
                        kk(i,j) = E*int(G,[lim_elem_2(v,1) lim_elem_2(v,2)]);
                        end
                    end
                end
                K_2{v,1}=kk;
            end
            %Tercera seccion
            A3 = pi*(-x+20); %funcion de area para seccion 3
             K_3 = cell(k,1);
            for v=1:k
                x1=lim_elem_3(v,1);
                x2=lim_elem_3(v,2);
                x3= ((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A3*diff(N(i))*diff(N(j))*diff(N(r));
                        kk(i,j) = E*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                        end
                    end
                end
                K_3{v,1}=kk;
            end
            
            %Cuarta seccion
            A4 = pi*(8)^2; %Area de seecion 4
            K_4 = cell(k,1);
            for v=1:k
                x1=lim_elem_4(v,1);
                x2=lim_elem_4(v,2);
                x3= ((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A4*diff(N(i))*diff(N(j))*diff(N(r));
                        kk(i,j) = E*int(G,[lim_elem_4(v,1) lim_elem_4(v,2)]);
                        end
                    end
                end
                K_4{v,1}=kk;
            end
            
                %ENSAMBLAJE DE LAS MATRICES DE RIGIDEZ, RIGIDEZ GLOBAL
                
                %Matriz Global de rigidez seccion 1
                Kglobal_1 = zeros((2*(k-1)+3),(2*(k-1)+3));
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Kglobal_1(i+2*(c-1),j+2*(c-1))=K_1{c}(i,j)+Kglobal_1(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
                
                %Matriz Global de rigidez seccion 2
                Kglobal_2 = zeros((2*(k-1)+3),(2*(k-1)+3));
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Kglobal_2(i+2*(c-1),j+2*(c-1))=K_2{c}(i,j)+Kglobal_2(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
                %Matriz Global de rigidez seccion 3
                Kglobal_3 = zeros((2*(k-1)+3),(2*(k-1)+3));
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Kglobal_3(i+2*(c-1),j+2*(c-1))=K_3{c}(i,j)+Kglobal_3(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
                %MatrizGlobal de rigidez seccion 4
                Kglobal_4 = zeros((2*(k-1)+3),(2*(k-1)+3));
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Kglobal_4(i+2*(c-1),j+2*(c-1))=K_4{c}(i,j)+Kglobal_4(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
               %Matriz Global de rigidez total
          
               %ponemos todas las matrices gobales en un arreglo tipo
               %celula
               CelK = cell(4,1);
               CelK{1} = Kglobal_1;
               CelK{2} = Kglobal_2;
               CelK{3} = Kglobal_3;
               CelK{4} = Kglobal_4;
               Kglobaltotal = zeros( (8*(k-1)+9),(8*(k-1)+9) );
               for c=1:4
                    for i= (1:(2*k+1))
                        for j= (1:(2*k+1))
                            Kglobaltotal(i+2*(k)*(c-1),j+2*(k)*(c-1))= CelK{c}(i,j)+Kglobaltotal(i+(2*k)*(c-1),j+(2*k)*(c-1));
                        end
                    end
               end
               
            %CALCULO DE LAS MATRICES DE INERCIA POR ELEMENTO FINITO
             A1 = pi*(-2*x+14)^(2);
             A2 = pi*(d_1/2)^2;
             A3 = pi*(-x+10);
             A4 = pi*(8)^2;
            %Matriz de inercia seccion 1
            M_1 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(3,3); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_1(v,1);
                x2=lim_elem_1(v,2);
                x3=((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A1*(N(i))*(N(j))*(N(r));
                        MM(i,j) = rho*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                        end
                    end
                end
                M_1{v,1}=MM;
            end
            
            %Matriz de inercia seccion 2
            M_2 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            
           for v=1:k
                x1=lim_elem_2(v,1);
                x2=lim_elem_2(v,2);
                x3=((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A2*(N(i))*(N(j))*(N(r));
                        MM(i,j) = rho*int(G,[lim_elem_2(v,1) lim_elem_2(v,2)]);
                        end
                    end
                end
                M_2{v,1}=MM;
            end
            
            %Matriz de inercia seccion 3
            M_3 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            
            for v=1:k
                x1=lim_elem_3(v,1);
                x2=lim_elem_3(v,2);
                x3=((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A3*(N(i))*(N(j))*(N(r));
                        MM(i,j) = rho*int(G,[lim_elem_3(v,1) lim_elem_3(v,2)]);
                        end
                    end
                end
                M_3{v,1}=MM;
            end
               
            %Matriz de inercia seccion 4  
            M_4 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
           
            for v=1:k
                x1=lim_elem_4(v,1);
                x2=lim_elem_4(v,2);
                x3=((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A4*(N(i))*(N(j))*(N(r));
                        MM(i,j) = rho*int(G,[lim_elem_4(v,1) lim_elem_4(v,2)]);
                        end
                    end
                end
                M_4{v,1}=MM;
            end
           
            %Matriz de emsamblaje seccion 1
            Mglobal_1 = zeros(2*(k-1)+3,2*(k-1)+3);
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Mglobal_1(i+2*(c-1),j+2*(c-1))=M_1{c}(i,j)+Mglobal_1(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 2
            Mglobal_2 = zeros(2*(k-1)+3,2*(k-1)+3);
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Mglobal_2(i+2*(c-1),j+2*(c-1))=M_2{c}(i,j)+Mglobal_2(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 3
            Mglobal_3 = zeros(2*(k-1)+3,2*(k-1)+3);
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Mglobal_3(i+2*(c-1),j+2*(c-1))=M_3{c}(i,j)+Mglobal_3(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 4
            Mglobal_4 = zeros(2*(k-1)+3,2*(k-1)+3);
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Mglobal_4(i+2*(c-1),j+2*(c-1))=M_4{c}(i,j)+Mglobal_4(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
                
           %Ensamblaje de matriz de inercia de toda la barra
           %ponemos todas las matrices gobales en un arreglo tipo
               %celula
               elK = cell(4,1);
               elK{1} = Mglobal_1;
               elK{2} = Mglobal_2;
               elK{3} = Mglobal_3;
               elK{4} = Mglobal_4;
               Mglobaltotal = zeros( (8*(k-1)+9),(8*(k-1)+9) );
               for c=1:4
                    for i= (1:(2*k+1))
                        for j= (1:(2*k+1))
                            Mglobaltotal(i+2*(k)*(c-1),j+2*(k)*(c-1))= elK{c}(i,j)+Mglobaltotal(i+2*(k)*(c-1),j+2*(k)*(c-1));
                        end
                    end
               end
               
         %Imposición de condiciones de frontera de la barra en cuestion
    Kglobaltotal_t = zeros(4*(k-1)+4,4*(k-1)+4);
    Mglobaltotal_t = zeros(4*(k-1)+4,4*(k-1)+4);
    for i=2:(4*(k-1)+5)
        for j=2:(4*(k-1)+5)
            Kglobaltotal_t(i-1,j-1)=Kglobaltotal(i,j);
            Mglobaltotal_t(i-1,j-1)=Mglobaltotal(i,j);
        end
    end          
       
disp('Matriz global de rigidez con condiciones de frontera')
    Kglobaltotal_t
    disp('Matriz global de inercia con condiciones de frontera')
    Mglobaltotal_t
    disp('Modos y Frecuencias naturales')
    
               %Calculo de modos y frecuencias naturales
               MI = inv (Mglobaltotal);
               KM = MI*Kglobaltotal;
              [L,V]= eig(KM)
                %Vector de modo L, Frecuencias naturales V
        
    elseif w == 2
        disp('Cuantos elementos finitos desea por seccion')
        fprintf('\n');
        k=input('cantidad de elementos'); %%cantidad de elementos finitos por barra
        dis_elem_1 = zeros(1,k); %vectores que guardan la longitud de los vectores de elementos finitos
        dis_elem_2 = zeros(1,k); % estos solo crean los vectores
        dis_elem_3 = zeros(1,k);
        dis_elem_4 = zeros(1,k);
        
        disp('Introducir longitudes de la primera seccion de la barra')
        fprintf('\n');
        for u=1:(k-1)
            y=input('Longitud de subsegmento n');
            dis_elem_1(u)=y %se deja ver para que el usuario pueda apreciar que valores va introduciendo
        end
        sum1 = 0; %suma para encontrar el valor del ultimo subsegmento
        for t=1:(k-1)
            sum1=sum1+dis_elem_1(t);
        end
        if sum1 >= L_1
            disp('Las medidas de sus elementos en la seccion 1 excede la longitud total del mismo, reinicie')
        end
        
        dis_elem_1(k)= L_1 - sum1;
        
        disp('Introducir longitudes de la segunda barra seccion de la barra')
        fprintf('\n');
        for u=1:(k-1)
            y=input('Longitud de subsegmento n');
            dis_elem_2(u)=y %se deja ver para que el usuario pueda apreciar que valores va introduciendo
        end
        sum2 = 0; %suma para encontrar el valor del ultimo subsegmento
        for t=1:(k-1)
            sum2=sum2+dis_elem_2(t);
        end
        if sum1 >= L_2
            disp('Las medidas de sus elementos en la seccion 2 excede la longitud total del mismo, reinicie')
        end
        
        dis_elem_1(k)= L_2 - sum2;
        
        disp('Introducir longitudes de la tercera seccion de la barra')
        fprintf('\n');
        for u=1:(k-1)
            y=input('Longitud de subsegmento n');
            dis_elem_3(u)=y %se deja ver para que el usuario pueda apreciar que valores va introduciendo
        end
        sum3 = 0; %suma para encontrar el valor del ultimo subsegmento
        for t=1:(k-1)
            sum3=sum3+dis_elem_3(t);
        end
        if sum3 >= L_3
            disp('Las medidas de sus elementos en la seccion 3 excede la longitud total del mismo, reinicie')
        end
        
        dis_elem_3(k)= L_3 - sum3;
        
        disp('Introducir longitudes de la cuarta seccion de la barra')
        fprintf('\n');
        for u=1:(k-1)
            y=input('Longitud de subsegmento n');
            dis_elem_4(u)=y %se deja ver para que el usuario pueda apreciar que valores va introduciendo
        end
        sum4 = 0; %suma para encontrar el valor del ultimo subsegmento
        for t=1:(k-1)
            sum1=sum1+dis_elem_4(t);
        end
        if sum4 >= L_4
            disp('Las medidas de sus elementos en la seccion 1 excede la longitud total del mismo, reinicie')
        end
        
        dis_elem_4(k)= L_4 - sum4;
      
        lim_elem_1 = zeros(k,2);
            lim_elem_2 = zeros(k,2);
            lim_elem_3 = zeros(k,2);
            lim_elem_4 = zeros(k,2);
        for i=0:(k-1)
            for j=1:2
                y=j-1;
            lim_elem_1(i+1,j) = [((i+y)*dis_elem_1(1))]; %discretización por elemento finito sec-1
            lim_elem_2(i+1,j) = [((i+y)*dis_elem_2(1)+(4/12))]; %discretización por elemento finito sec-2 
            lim_elem_3(i+1,j) = [((i+y)*dis_elem_3(1)+(10/12))]; %discretización por elemento finito sec-3
            lim_elem_4(i+1,j) = [((i+y)*dis_elem_4(1)+(1))]; %discretización por elemento finito sec-4
            end
        end
        %Funciones de interpolación
            syms x %x1 limite inferior x2 limite superior
            x1=1;
            x2=2;
            x3=((x1+x2)/2);
            N = [(((x-x2)*(x-x3))/((x1-x2)*(x1-x3))),(((x-x1)*(x-x3))/((x2-x1)*(x2-x3))),((x-x1)*(x-x2))/((x3-x1)*(x3-x2))];
            
            %calculo de las matrices de Rigidez por elemento finito
            A1 = pi*(-2*x+14)^(2);
            %Primera seccion
            K_1 = cell(k,1);
            kk = zeros(3,3);
            for v=1:k
                x1=lim_elem_1(v,1);
                x2=lim_elem_1(v,2);
                x3= ((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A1*diff(N(i))*diff(N(j))*diff(N(r));
                        kk(i,j) = E*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                        end
                    end
                end
                K_1{v,1}=kk;
            end
            
            %Segunda Seccion
            A2 = pi*(d_1/2)^2; %Area constante de la seccion 2
            K_2 = cell(k,1);
            kk = zeros(3,3);
            for v=1:k
                x1=lim_elem_2(v,1);
                x2=lim_elem_2(v,2);
                x3= ((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A2*diff(N(i))*diff(N(j))*diff(N(r));
                        kk(i,j) = E*int(G,[lim_elem_2(v,1) lim_elem_2(v,2)]);
                        end
                    end
                end
                K_2{v,1}=kk;
            end
            %Tercera seccion
            A3 = pi*(-x+20); %funcion de area para seccion 3
             K_3 = cell(k,1);
            for v=1:k
                x1=lim_elem_3(v,1);
                x2=lim_elem_3(v,2);
                x3= ((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A3*diff(N(i))*diff(N(j))*diff(N(r));
                        kk(i,j) = E*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                        end
                    end
                end
                K_3{v,1}=kk;
            end
            
            %Cuarta seccion
            A4 = pi*(8)^2; %Area de seecion 4
            K_4 = cell(k,1);
            for v=1:k
                x1=lim_elem_4(v,1);
                x2=lim_elem_4(v,2);
                x3= ((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A4*diff(N(i))*diff(N(j))*diff(N(r));
                        kk(i,j) = E*int(G,[lim_elem_4(v,1) lim_elem_4(v,2)]);
                        end
                    end
                end
                K_4{v,1}=kk;
            end
            
                %ENSAMBLAJE DE LAS MATRICES DE RIGIDEZ, RIGIDEZ GLOBAL
                
                %Matriz Global de rigidez seccion 1
                Kglobal_1 = zeros((2*(k-1)+3),(2*(k-1)+3));
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Kglobal_1(i+2*(c-1),j+2*(c-1))=K_1{c}(i,j)+Kglobal_1(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
                
                %Matriz Global de rigidez seccion 2
                Kglobal_2 = zeros((2*(k-1)+3),(2*(k-1)+3));
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Kglobal_2(i+2*(c-1),j+2*(c-1))=K_2{c}(i,j)+Kglobal_2(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
                %Matriz Global de rigidez seccion 3
                Kglobal_3 = zeros((2*(k-1)+3),(2*(k-1)+3));
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Kglobal_3(i+2*(c-1),j+2*(c-1))=K_3{c}(i,j)+Kglobal_3(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
                %MatrizGlobal de rigidez seccion 4
                Kglobal_4 = zeros((2*(k-1)+3),(2*(k-1)+3));
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Kglobal_4(i+2*(c-1),j+2*(c-1))=K_4{c}(i,j)+Kglobal_4(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
               %Matriz Global de rigidez total
          
               %ponemos todas las matrices gobales en un arreglo tipo
               %celula
               CelK = cell(4,1);
               CelK{1} = Kglobal_1;
               CelK{2} = Kglobal_2;
               CelK{3} = Kglobal_3;
               CelK{4} = Kglobal_4;
               Kglobaltotal = zeros( (8*(k-1)+9),(8*(k-1)+9) );
               for c=1:4
                    for i= (1:(2*k+1))
                        for j= (1:(2*k+1))
                            Kglobaltotal(i+2*(k)*(c-1),j+2*(k)*(c-1))= CelK{c}(i,j)+Kglobaltotal(i+(2*k)*(c-1),j+(2*k)*(c-1));
                        end
                    end
               end
               
            %CALCULO DE LAS MATRICES DE INERCIA POR ELEMENTO FINITO
             A1 = pi*(-2*x+14)^(2);
             A2 = pi*(d_1/2)^2;
             A3 = pi*(-x+10);
             A4 = pi*(8)^2;
            %Matriz de inercia seccion 1
            M_1 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            MM = zeros(3,3); %Matrices de inercia del elemento
            for v=1:k
                x1=lim_elem_1(v,1);
                x2=lim_elem_1(v,2);
                x3=((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A1*(N(i))*(N(j))*(N(r));
                        MM(i,j) = rho*int(G,[lim_elem_1(v,1) lim_elem_1(v,2)]);
                        end
                    end
                end
                M_1{v,1}=MM;
            end
            
            %Matriz de inercia seccion 2
            M_2 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            
           for v=1:k
                x1=lim_elem_2(v,1);
                x2=lim_elem_2(v,2);
                x3=((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A2*(N(i))*(N(j))*(N(r));
                        MM(i,j) = rho*int(G,[lim_elem_2(v,1) lim_elem_2(v,2)]);
                        end
                    end
                end
                M_2{v,1}=MM;
            end
            
            %Matriz de inercia seccion 3
            M_3 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
            
            for v=1:k
                x1=lim_elem_3(v,1);
                x2=lim_elem_3(v,2);
                x3=((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A3*(N(i))*(N(j))*(N(r));
                        MM(i,j) = rho*int(G,[lim_elem_3(v,1) lim_elem_3(v,2)]);
                        end
                    end
                end
                M_3{v,1}=MM;
            end
               
            %Matriz de inercia seccion 4  
            M_4 = cell(k,1); %Celula para guardar matrices de inercia de seccion 1
           
            for v=1:k
                x1=lim_elem_4(v,1);
                x2=lim_elem_4(v,2);
                x3=((x1+x2)/2);
                for i=1:3
                    for j=1:3
                        for r=1:3
                        G =A4*(N(i))*(N(j))*(N(r));
                        MM(i,j) = rho*int(G,[lim_elem_4(v,1) lim_elem_4(v,2)]);
                        end
                    end
                end
                M_4{v,1}=MM;
            end
           
            %Matriz de emsamblaje seccion 1
            Mglobal_1 = zeros(2*(k-1)+3,2*(k-1)+3);
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Mglobal_1(i+2*(c-1),j+2*(c-1))=M_1{c}(i,j)+Mglobal_1(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 2
            Mglobal_2 = zeros(2*(k-1)+3,2*(k-1)+3);
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Mglobal_2(i+2*(c-1),j+2*(c-1))=M_2{c}(i,j)+Mglobal_2(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 3
            Mglobal_3 = zeros(2*(k-1)+3,2*(k-1)+3);
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Mglobal_3(i+2*(c-1),j+2*(c-1))=M_3{c}(i,j)+Mglobal_3(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
            %Matriz de emsamblaje seccion 4
            Mglobal_4 = zeros(2*(k-1)+3,2*(k-1)+3);
                for c=1:k
                    for i=1:3
                      for j=1:3
                        Mglobal_4(i+2*(c-1),j+2*(c-1))=M_4{c}(i,j)+Mglobal_4(i+2*(c-1),j+2*(c-1));
                      end
                    end
                end
                
           %Ensamblaje de matriz de inercia de toda la barra
           %ponemos todas las matrices gobales en un arreglo tipo
               %celula
               elK = cell(4,1);
               elK{1} = Mglobal_1;
               elK{2} = Mglobal_2;
               elK{3} = Mglobal_3;
               elK{4} = Mglobal_4;
               Mglobaltotal = zeros( (8*(k-1)+9),(8*(k-1)+9) );
               for c=1:4
                    for i= (1:(2*k+1))
                        for j= (1:(2*k+1))
                            Mglobaltotal(i+2*(k)*(c-1),j+2*(k)*(c-1))= elK{c}(i,j)+Mglobaltotal(i+2*(k)*(c-1),j+2*(k)*(c-1));
                        end
                    end
               end
               
       %Imposición de condiciones de frontera de la barra en cuestion
    Kglobaltotal_t = zeros(4*(k-1)+4,4*(k-1)+4);
    Mglobaltotal_t = zeros(4*(k-1)+4,4*(k-1)+4);
    for i=2:(4*(k-1)+5)
        for j=2:(4*(k-1)+5)
            Kglobaltotal_t(i-1,j-1)=Kglobaltotal(i,j);
            Mglobaltotal_t(i-1,j-1)=Mglobaltotal(i,j);
        end
    end
    disp('Matriz global de rigidez con condiciones de frontera')
    Kglobaltotal_t
    disp('Matriz global de inercia con condiciones de frontera')
    Mglobaltotal_t
    disp('Modos y Frecuencias naturales')
               %Calculo de modos y frecuencias naturales
               MI = inv (Mglobaltotal);
               KM = MI*Kglobaltotal;
              [L,V]= eig(KM)
                %Vector de modo L, Frecuencias naturales V
    else 
         fprintf('\n');
         disp('**El valor introducido no es valido**')
    end
else 
    fprintf('\n');
    disp('**El valor introducido no es valido**') 
end

    
