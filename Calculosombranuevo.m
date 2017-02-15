
%% Rutina para obtener los puntos de sombra e insolacion en cada punto de un dem, en una fecha determinada previamente por el usuario sacada de la rutina UbiSolarMarsGrilla.m (fuente Libro Felic�simo rutina Perfilsol)

% Distancia de radio de influencia para interpolacion
dist=.02;

% Primero se levanta el dem (en caso de no utilizar la rutina anterior)
% Ubicación donde esta ubicado el dem
% path='E:\Publicaciones\Mars habitat\DemEjemplo\';
path='E:\Publicaciones\Mars habitat\DemEjemplo\';
cd([path '\Sombras'])
% Buscamos todos los archivos .tiff que estan en la carpeta
dem=dir([path '*.tif']);

% Nombre del DEM que utilizaremos para guardar los datos
nom2save=strrep(dem.name,'.tif','');

% Cargamos el DEM de interes
[A, R]=geotiffread([path dem.name]);
A=double(A);
% Creo las matrices de latitud y longitud
Longitudes=((R.XLimWorld(1):R.DeltaX:R.XLimWorld(2))-180) * (-1);
Latitudes=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1);

[Lon, Lat]=meshgrid(Longitudes, Latitudes);
tam=size(Lat);

clear Longitudes Latitudes

% calculo de la escala en metros
RadioM=3396190.0;
Perim= 2*pi*RadioM;
MxDeg= Perim/360;
ESCALA=(R.DeltaX)*MxDeg;

% Buscamos los archivos .mat
fecha=dir([path '*.mat']);


for k=1:length(fecha)
    load(fecha(k).name);
    %         load('dem_736696-0417.mat')
    
    % Creo la matriz que alojara los resultados de cada pixel
    % S es el Vector resultado, con las coordenadas x y luz (longitud, latitud y luz)
    % del perfil que estamos buscando
    
    S= zeros(tam(1),tam(2));  
    % % %
    
    tic
    for i=10:20%tam(1)
        for j=10:20%tam(2)
            % variable auxiliar de rutina de sombra
            l=2;
            %             if D.Zh(i,j)>0
            
    Distmax=A/tand(D.Zh(i,j));
    
%     agregar medir la distancia entre el ij  y cada punto de la matriz

trans=Distmax-Distpix2pix;
[idx d]=find(trans>0);
if ~isempty(idx)



    
            %             Cambio de fila para cada punto del perfil (real)
            
                INCF = -cosd(D.Az(i,j))*(R.DeltaX*2);
            %             Cambio de columna para cada punto del perfil (real)
             INCC = sind(D.Az(i,j))*(R.DeltaX*2);
          
            
            
            
            %             Distancia real entre el punto inicial del perfil
            %             y el punto actual en m
            
            DISTANCIA = 0;
            % NP en el n�mero de puntos del perfil
            %             NP = 1;
            
            
            % Definimos la latitud y longitud del punto inicial del perfil
            Latit = Lat(i,j);
            % Definimos la longitud y longitud del punto inicial del perfil
            Longit = Lon(i,j);
            
            % Vector resultado, con las coordenadas x y z (longitud, latitud y altura)
            % del perfil que estamos buscando
            V=[];
            
            %             Generaci�n del perfil topografico en la linea de azimuth
            %             solar.
            %X son las longitudes limites del dem
            xmax=((R.XLimWorld(1)-180)  * (-1));
            xmin=((R.XLimWorld(2)-180)  * (-1));
            % Y son las latitudes limites del dem
            ymin=R.YLimWorld(1);
            ymax=R.YLimWorld(2);
            
            % El primer punto de nuestro perfil es:
            V=[V; Latit Longit A(i,j) DISTANCIA];
            
            while Latit >= ymin  && Latit <= ymax && Longit >= xmin && Longit <= xmax
                
                Latit=Latit-INCF;
                Longit=Longit-INCC;
                % %                    Agregar la interpolacion bilineal de los puntos que est�n
                %             dentro
                [alt, altStd]=interpolacione(Lon(:),Lat(:), A(:), Longit, Latit, dist);
                
                
                DISTANCIA = DISTANCIA+ESCALA;
                 V=[V; Latit Longit alt DISTANCIA];
            end
             DISTANCIA = DISTANCIA+ESCALA;
            alt=A(i,j)+(DISTANCIA*(tand(D.Zh(i,j))));
            V=[V; Latit Longit alt DISTANCIA];
            
                end
    end
    toc
            
            
            %             Calculo de sombras Opcion1
%             %{
            NP= size(V,1);
            DIST1 = sqrt (((V (1,2) - V (NP,2))^2) + ((V (1,1) - V (NP,1))^2));
            
            %             end
            
            if DIST1>0.0078
                TG1 = (V(NP,3)- V(1,3))/(DIST1*R.DeltaX);
                
                while l>=2 && l< NP
                    
                    DIST2 = sqrt (((V (1,2) - V (l,2))^2) + ((V (1,1) - V (l,1))^2));
                    TG2 = (V(l,3)- V(1,3))/(DIST2*R.DeltaX);
                    if TG1<TG2
                        break
                    end
                    l=l+1;
                    
                    clear DIST2 TG2
                end
                if l==NP
                    S(i,j)=1;
                    
                end
                clear l
            end
            clear TG1 DIST1 NP V
            %}
            
            % Calculo de sombras Opcion 2
            %{
            altmax=max(V(:,3));
            Zmax=tand(D.Zh(a,j))/altmax;
            aux=V(:,3)==altmax;
                        
            if Zmax>V(aux, 4)
                S(i,j)=1;
            end
            %}
        end
    end
    toc
    save([nom2save '_' fecha2save '_Sombra' ], 'S')
end

