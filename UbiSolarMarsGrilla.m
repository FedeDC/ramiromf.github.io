%% Rutina para obtener los valores de Azimut y Elevaciï¿½n del Sol sobre el horizonte para un tiempo real (ej. 25 enero 2024 H:M:S) en una determinada Latititud y Longititud en Marte

% Periodo de fechas de interes del estudio
% Fecha inicial
preg='Â¿Fecha inicial? ("HH:MM:SS mmm dd yyyy")  ';
ini=input(preg);
ini=datenum(ini, 'HH:MM:SS mmm dd yyyy');

% Fecha final
preg='Â¿Fecha final? ("HH:MM:SS mmm dd yyyy")  ';
fin=input(preg);
fin=datenum(fin, 'HH:MM:SS mmm dd yyyy');

% Una Hora en tiempo serial
hourinserial = datenum(2010,1,1,1,0,0)-datenum(2010,1,1,0,0,0);

%Periodo de tiempo de interes
timelapse=ini:hourinserial:fin;

clear ini fin hourinserial
%% Fecha J2000 (la diferencia del 1:09 minutos corresponde a la diferencia
% entre TAI-UTC (http://hpiers.obspm.fr/) y TT-UTC (32.184seg)
J2000=datenum('11:58:56 Jan 01 2000', 'HH:MM:SS mmm dd yyyy');

% UbicaciÃ³n donde esta ubicado el dem
path='E:\Publicaciones\Mars habitat\DemEjemplo\';
cd(path)
% Buscamos todos los archivos .tiff que estan en la carpeta
files=dir([path '*.tif']);

% Nombre del DEM que utilizaremos para guardar los datos
nom2save=strrep(files.name,'.tif','');

% Cargamos el DEM de interes
[A, R]=geotiffread([path files.name]);

% Creo las matrices de latitud y longitud
Longitudes=((R.XLimWorld(1):R.DeltaX:R.XLimWorld(2))-180) * (-1);
Latitudes=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1);

[Lon, Lat]=meshgrid(Longitudes, Latitudes);
tam=size(Lat);

clear Longitudes Latitudes

 tic
%% Determinando el periodo de interes
% Our Mars time calcuLatitions will use the parameter ?tJ2000: the elapsed
% time in days since the J2000 epoch, i.e., 12:00 on Jan. 1, 2000 (TT).
for k=1:length(timelapse)
    %% Hasta donde queremos integrar
    % TT=datenum('13:46:31 Jan 03 2004', 'HH:MM:SS mmm dd yyyy');
    TT=timelapse(k);
    % JDTT=juliandate('00:00:00 Jan 06 2000', 'HH:MM:SS mmm dd yyyy');
    %             RECORDAR MODIFICAR EL VALOR 69.184 SEGÚN EL AÑO DE ESTUDIO TAI-UTC=37 para 2017 (http://hpiers.obspm.fr/) + TT-UTC 32.184seg
    % 01/01/1999 32s
    % 01/01 2006 33s
    % 01/01 2009 34s
    % 01/07 2012 35s
    % 01/07 2015 36s
    % 01/01 2017 37s
    %
    %             (http://hpiers.obspm.fr/)
    JDTT=juliandate(TT+(69.184/86400));
    
    % Cantidad de dï¿½as desde el J2000
    dtJ2000=TT-J2000;
    
    %% Determinando parametros marcianos
    
    % Determinamos Mars mean anomaly. (AM2000, eq. 16)
    
    M = 19.3871 + 0.52402073 * dtJ2000;
    
    % Determinamos el angulo de Fiction Mean Sun. (AM2000, eq. 17)
    
    alphaFMS = 270.3871 + 0.524038496 * dtJ2000;
    
    % Determinamos las perturbers. (AM2000, eq. 18)
    % Datos Auxiliares;
    A=[0.0071, 0.0057, 0.0039 , 0.0037, 0.0021, 0.0020, 0.0018];
    tau=[2.2353, 2.7543, 1.1177, 15.7866, 2.1354, 2.4694, 32.8493];
    phi=[49.409, 168.173, 191.837, 21.736, 15.704, 95.528, 49.095];
    PBS=0;
    for n=1:7;
        trnsit = A(n) * cosd(((0.985626 * dtJ2000) / tau(n) )+ phi(n));
        PBS = PBS +trnsit;
    end
    % donde 0.985626ï¿½ = 360ï¿½ / 365.25
    
    % Determinaciï¿½n de la Equation of Center. (Bracketed term in AM2000, eqs. 19 and 20)
    
    % La ecuacion del centro es la true anomaly (v) minus mean anomaly (M).
    
    vminusM = (10.691 + 3.0*10^-7 * dtJ2000) * sind(M) + 0.623 * sind(2*M) + 0.050 * sind(3*M) + 0.005 * sind(4*M) + 0.0005 * sin(5*M) + PBS;
    
    % Determinacion areocentric solar Longititude (). (AM2000, eq. 19)
    
    Ls = alphaFMS + (vminusM);
    
    %% Calculos orbitales del SOL
    % Parametros orbitales
    % Declinaciï¿½n solar;
    delta = asind (0.42565 * sind(Ls)) + 0.25 * sind(Ls);
    
    % Angulo Horario;
    % Determinacion de la Ecuacion del tiempo. (AM2000, eq. 20)
    
    EOT = (2.861 * sind(2*Ls) - 0.071 * sind(4*Ls) + 0.002 * sind(6*Ls) - (vminusM)) * (24/360);
    
    % El resultado de arriba esta en grados. Multiply by (24 h / 360ï¿½) = (1 h / 15ï¿½) to obtain the result in hours.
    
    % Determinacion de el Coordinated Mars Time, i.e., Airy Mean Time. (AM2000, eq. 22, modified)
    % Este es el tiempo solar medio en el meridiano 0 de Marte.
    
    MTC = mod(24 * ( ((JDTT - 2451549.5) / 1.0274912517) + 44796.0 - 0.0009626 ), 24);
    
    % % Determinamos la hora solar local media (Local Mean Solar Time).
    % % The Local Mean Solar Time for a given planetographic Longititude, ?, in
    % % degrees west, is easily determined by offsetting from the mean solar time
    % % on the prime meridian.
    %
    % LMST = MTC - ? (24 h / 360ï¿½) = MTC - ? (1 h / 15ï¿½)
    %
    % C-4. Determine Local True Solar Time. (AM2000, eq. 23)
    %
    % LTST = LMST + EOT (24 h / 360ï¿½) = LMST + EOT (1 h / 15ï¿½)
    
    % Determinacion subsolar Longititude
    
    Longits = ((MTC + EOT) * 15 ) - 180;
    %             Longits = 4.70500;
    
    % Creo la matriz que alojara los resultados de cada pixel
    
    D.Az= nan*ones(tam(1),tam(2));
    D.Z= nan*ones(tam(1),tam(2));
    D.Zh = nan*ones(tam(1),tam(2));
    D.delta = nan*ones(tam(1),tam(2));
    D.DMS = nan*ones(tam(1),tam(2));
    
    % Por cada punto de la grilla del dem obtendremos los parametros de
    % interes.
    for i=1:tam(1)
       
        for j=1:tam(2)
            %% Ahora haremos el recorrido por cada punto de la grilla.
            
            %%  ? is the planetographic Longititude, Decimal degrees.000.0000
            Longit = Lon(i,j);
            % Latititud de interes;
            Latit = Lat(i,j);
            
            % Determinacion de angulo horario
            H= Longit - Longits;
            
            % Determinaciï¿½n de la elevaciï¿½n solar local
            
            % PAra cualquier punto en la superficie de Marte. El angulo zenital es:
            
            Z = acosd (sind(delta)* sind(Latit) + cosd(delta) * cosd(Latit)* cosd(H));
            
            % donde Latit es la Latititud planetografica, Longit es la Longititud planetografica, y H  es la angulo horario, Longit - Longits.
            
            % La elevacion del sol es 90ï¿½ - Z.
            
            Zh = 90 - Z;
            
            % Determinacion del azimuth solar local
            
            % As = (atand (sind(H) / ((cosd(Latit) * tand(delta)) - (sind(Latit) * cosd(H)))))+ 180;
            Az = atan2d (sind(H) , ((cosd(Latit) * tand(delta)) - (sind(Latit) * cosd(H))));
            if Az<0
                Az = Az + 360;
            end
            
            % Determinar la distancia heliocentrica(AM2000, eq. 25, corrected)
            
            DMS = 1.52367934 * (1.00436 - (0.09309 * cosd (M)) - (0.004336 * cosd (2*M)) - (0.00031 * cosd (3*M)) - (0.00003 *cos (4*M)));
            
            
            D.Az(i,j)=Az;
            D.Z(i,j)=Z;
            D.Zh(i,j)=Zh;
            D.DMS(i,j)=DMS;
            D.delta(i,j)=delta;
            
        end
  
        
        clear Longit Latit H Az Z Zh DMS
    end
    fecha2save=strrep(num2str(timelapse(k)),'.','-');
    save([nom2save '_' fecha2save ], 'D')
    clear D  Ls M PBS A TT JDTT alphaFMS tau phi trnsit vminusM delta
end

% Para exportar la tabla a un txt separado por tab (\t)
% j=[Lat(:) Lon(:) D.delta(:)];
% dlmwrite('prueba.txt',j,'delimiter','\t')
%%%%Para transformar el formato string a fecha de nuevo reemplazar ### por el número de string 
%%%%datestr(###);
toc