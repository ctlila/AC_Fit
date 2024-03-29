                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%                            %%%%%%%%%%%%%
                %%%%%%%%%%%%%         FIT AC DATA        %%%%%%%%%%%%%
                %%%%%%%%%%%%%                            %%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   %%%% Implemented by Lila CT in april 2021 %%%%

                                   
%Last modifications are commented %2.5

clear

format longG

state = rng ; 


%% OPTIONAL FEATURES

w = 0.8 ; % weight of chi2 compared to chi1 in the fitting process.

lambda = 0.001 ; %Lagrange parameter for Lasso regularization of parameters with norm l1

%Do you want to plot each temperatures separately or all temperatures on a single figure ?
plot_option = 'single' ; %options : 'single' / 'separated'                                 

%Do you want to plot and fit all the data of your file (false) or only one (true) ? 
only_one = false ;
which = 11 ; %if only one, which one ? (row index)

%Are there rows to ignore at the beginning of the matrix ? 
%If yes specify the first row to be taken into account (frst=2 means the first measure is ignored).
%If no the correct value is 1.
frst = 1; %default value for first row (no row to ignore)


%Initial parameters of the fits + typical patterns of AC curves. You may modify/add patterns
%Generally no change in initial parameters or patterns is needed.

%order : (chiT-chiS)_1 tau_1 alpha_1 chiS  (chiT-chiS)_2 tau_2 alpha_2
param2 = [0.226 5.93e-6 0.192 0.016 0.1 2.95e-4 0.075];
param1 = [0.745 0.0012 0.238 0.012];

%patterns with 1 contrubition
pattern1 = [0.87 0.007 0.156 0.012];
pattern5 = [0 0 0 0];  

%pattern6 = [0 0 0 0]; %uncomment here to add pattern6 (+all the lines containing this instruction)

%patterns with 2 contributions
pattern2 = [0.63 0.0015 0.132 0.5 0.325 0.0 0.133];
pattern3 = [0.28 4.74e-5 0.148 0.012 0.14 0.0025 0.096];
pattern4 = [0.226 5.93e-6 0.192 0.012 0.1 2.95e-4 0.075];

%Please note that no fit is performed on chiS_2 and we just assume
%chiS_1=chiS_2 because the only value obtained with the fit is chiS_1+chiS_2
%so we cannot separate the two contributions. However this is physically correct


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%%% RUN %%%
%%%%%%%%%%%

%% Inputs

%Instrument choice
instr = input("ppms or squid ?",'s');
if isempty(instr)
    instr = 'ppms';
end

%data file
if contains(instr ,'ppms')
    rawfile=input("Name (and path) of your raw data file from PPMS",'s');
    if isempty(rawfile)
        rawfile = 'FSS2-65A2p1_ScanT_0.1T_test.dat';
    end
elseif contains(instr,'squid')
    rawfile=input("Name (and path) of your raw data file from SQUID",'s');
    if isempty(rawfile)
        rawfile = 'KV788_Lurad_acscanH.ac.dat';
    end
end

[filepath,filename,~] = fileparts(rawfile) ;
if isempty(filepath)
    filepath = '.' ;
end

%Molecular weight of your compound and mass of the sample !!!ALL IN mg!!!
MW=input('Molecular Weight of your compound in g/mol');
if isempty(MW)
     MW = 838.81; %enter here a default value for the molecular weight.
end 

mass=input('Mass of the sample in mg');
if isempty(mass)
     mass = 10.2; %=default value for sample mass.
end 

mteflon=input('Mass of teflon in mg');
if isempty(mteflon)
     mteflon = 8.17; %default value for teflon mass.
end 

%importing data
Rawdata = readmatrix(rawfile);
[nbrow,~] = size(Rawdata); %2.5

disp(['---------Data imported from',rawfile,'---------']) %2.5

%% Definitions and preallocations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%      MAIN      %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%fit functions with 1 or 2 contributions for chiprime and chisecond, don't touch

fitpr1 = @(x,xdata)x(4)+x(1)*(1+((2*pi*xdata*x(2)).^(1-x(3)))*sin((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3))));
fitpr2 = @(x,xdata)x(4)+x(1)*(1+((2*pi*xdata*x(2)).^(1-x(3)))*sin((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3)))) + x(4)+x(5)*(1+((2*pi*xdata*x(6)).^(1-x(7)))*sin((pi*x(7))/2))./(1+(2*((2*pi*xdata*x(6)).^(1-x(7)))*sin(pi*(x(7)/2)))+((2*pi*xdata*x(6)).^(2-2*x(7))));

fitsec1= @(x,xdata)(x(1)*((2*pi*xdata*x(2)).^(1-x(3)))*cos((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3))));
fitsec2 = @(x,xdata)(x(1)*((2*pi*xdata*x(2)).^(1-x(3)))*cos((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3))))+(x(5)*((2*pi*xdata*x(6)).^(1-x(7)))*cos((pi*x(7))/2))./(1+(2*((2*pi*xdata*x(6)).^(1-x(7)))*sin(pi*(x(7)/2)))+((2*pi*xdata*x(6)).^(2-2*x(7))));



%Calculates the number of points for each acquisition nbpt %2.5
%Ensure all measure sets have the same 
if contains(instr,'squid','IgnoreCase',true)
    m=max(Rawdata(:,15));
    nbocc=length(find(Rawdata(:,15)==m));
    nbpt=nbrow/nbocc;
elseif contains(instr,'ppms','IgnoreCase',true)
    m=max(Rawdata(:,5));
    nbocc=length(find(Rawdata(:,5)==m));
    nbpt=nbrow/nbocc;
end

%specifications of your data file : the number of measure sets.
nbfit = (nbrow-(frst-1))/nbpt;

if ~(floor(nbfit)==nbfit)
    disp('Wrong entry. Please ensure your .dat file is clean (no additionnal row) and all acquisitions have the same length.')
    return 
end

%Import relevant columns
if contains(instr,'squid','IgnoreCase',true)
for i = 1:nbfit %Reordering the Rawdata matrix to have incresing frequencies for each measure set.  %2.5
    Rawdata((frst+(i-1)*nbpt):(frst-1+i*nbpt),:) = sortrows(Rawdata((frst+(i-1)*nbpt):(frst-1+i*nbpt),:),15);
end
    temperatures = Rawdata(:,4);
    amplitudes = Rawdata(:,14);
    fields = Rawdata(:,3);
    mpr = Rawdata(:,5);
    stddev = (Rawdata(:,6)+Rawdata(:,8))./2;
    fq = Rawdata(:,15);
    msec = Rawdata(:,7);
elseif contains(instr,'ppms','IgnoreCase',true)
for i = 1:nbfit %Reordering the Rawdata matrix to have incresing frequencies for each measure set.  %2.5
    Rawdata((frst+(i-1)*nbpt):(frst-1+i*nbpt),:) = sortrows(Rawdata((frst+(i-1)*nbpt):(frst-1+i*nbpt),:),5);
end
    temperatures = Rawdata(:,3);
    amplitudes = Rawdata(:,6);
    fields = Rawdata(:,4);
    mpr = Rawdata(:,9);
    stddev = Rawdata(:,8);
    fq = Rawdata(:,5);
    msec = Rawdata(:,10);
end

%ScanT or ScanH : Is it a scan in Temperature (option T) or in Magnetic field (option H)? %2.5
if range(temperatures)<range(fields)/100
    TorH='H';
else 
    TorH='T';
end


%Now I just do a little trick to calculate the variable stepH that represent the smallest difference between two consecutivefields in the data. It is only used to build the legend.
%The only restriction is that you cannot have a difference of less than 10 Oe between two fields and less than 0.1K between two temperatures. It should not be too restrictive as it is close to the machine precision. 

if TorH == 'T'   %2.5
    stepT=0.1;
    stepH=100;
else 
    for i = 1:nbfit
    meansH(i)= mean(fields(((i-1)*nbpt+1):(i*nbpt)));
    end
    dH=diff(meansH);
    mH=min(dH);
    stepH=round(mH/10)*10;
    stepT=0.1;
end


%Calculation of chi' and chi" from M' and M" respectively.
chisec = msec./amplitudes.*(MW*10^3/mass);
chipr = ((mpr+(3.7*10^-10 .* mteflon .*amplitudes))./amplitudes).* MW.*10^3./mass + MW.*10^-6./2 ;
%Note that the diamagnetic contribution is only roughly calculated. You can modify this formula to have better accuracy. 

%Get the standard deviations of each measure for weights
std_devstock= stddev;


%options for the optimization algorithm
optionsnonlin = optimoptions('lsqnonlin','MaxFunctionEvaluations',1.2e+05,'MaxIterations',4e4,'FunctionTolerance',1e-10,'StepTolerance',1e-10,'OptimalityTolerance',1e-8,'UseParallel',false,'Display','none'); %2.5 (just the display options)

%lower bounds ans upper bounds of parameters
lb2=[0,0,0,0,0,0,0]; 
lb1 = [0,0,0,0];
ub2=[];
ub1 = [];

%for colormap 
newcolors = jet(nbfit);


%initializing the matrixes (preallocation)                      

choice = zeros(nbfit,1);
xdatastock = zeros(nbpt,nbfit);                       %Storage for each fit xdata (frequences) and y data (chi' and chi")
ydatasecstock = zeros(nbpt,nbfit);
ydataprstock = zeros(nbpt,nbfit);
tempparamwstock = zeros(nbfit, 7);                    %Storage for the temporary parameters of weighted fit
tempparamuwstock = zeros(nbfit, 7);                   %Storage for the temporary parameters of unweighted fit
residualstock = zeros(nbpt,nbfit);                    %Storage for residuals on each point for each fit
resnormstock = zeros(nbfit,1);                        %Storage for norm of residuals for each fit.
errors = cell (1,nbfit);                              %Storage for errors
J = cell(1,nbfit);                                    %Storage for jacobian matrixes reshaped
fitparamstock = zeros(nbfit,23);                      %Storage for fitted parameters and temperature/field for each fit
pointssec = gobjects(1,nbfit);                        %Storage for chi" graph
pointspr = gobjects(1,nbfit);                         %Storage for chi' graph
legendtemp = cell(1,nbfit);                           %Storage for legends
legendfield = cell(1,nbfit);
jacobian = zeros(nbpt);                               %Storage for raw jacobians
outputfile = zeros(nbfit*nbpt,7);                     % x and y data of each curve and each fit + T + H
Column_names = {'d1','errd1','t1','errt1','a1','erra1','chiS1','errchiS1','d2','errd2','t2','errt2','a2','erra2','chiS2','errchiS2','T','1/T','H','ln(t1)','err(ln(t1))','ln(t2)','err(ln(t2))'} ;


if only_one 
    end_index = which+1 ;
    i = which ;
    plot_option = 'single' ;
else 
    i = 1 ;
    end_index = nbfit + 1;
end    

%% Loops

    while i<end_index

        clear xdata ydata xdatalog fitparam W jacobian

        xdata = fq((frst+(i-1)*nbpt):(frst-1+i*nbpt)); %Frequencies

        ydatasec = chisec((frst+(i-1)*nbpt):(frst-1+i*nbpt),1); %Chi"
        
        ydatapr = chipr((frst+(i-1)*nbpt):(frst-1+i*nbpt),1); %Chi'

        Temps = temperatures((frst+(i-1)*nbpt):(frst-1+i*nbpt)); %temperature
        
        Temp = mean(Temps);
        
        rTemp = round(Temp/stepT)*stepT; %round temperature

        Fields = fields((frst+(i-1)*nbpt):(frst-1+i*nbpt)); %magnetic field
        
        Field = mean(Fields);
        
        rField = round(Field/stepH)*stepH ; %round field

                

        
        
        %USER CHOICE ON 1 or 2 contributions (or ignore this acquisition).
        
            figure('Name','temp_data','units','normalized','outerposition',[0.6 0 0.4 1]) ;
            upplot = subplot(2,1,1) ;
            downplot = subplot(2,1,2) ;
            %tempaxes = axes;

            if contains(TorH,'T')
                semilogx(downplot,xdata,ydatasec,'o','DisplayName', 'Data ' + string(rTemp)+' K','color','k');
                hold(downplot,'on')
            else
                semilogx(downplot,xdata,ydatasec,'o','DisplayName', 'Data ' + string(rField)+' Oe','color','k');
                hold(downplot,'on')
            end
            
            ylabel(downplot,'$\chi$" (emu/mol)','Interpreter','latex');
            xlabel(downplot,'Frequence (Hz)', 'Interpreter', 'latex');
            legend(downplot,'Location','Best');


            if contains(TorH,'T')
                semilogx(upplot,xdata,ydatapr,'o','DisplayName',"Data " + string(rTemp)+" K",'color','k');
                hold(upplot,'on')
            else
                semilogx(upplot,xdata,ydatapr,'o','DisplayName', "Data " + string(rField)+' Oe','color','k');
                hold(upplot,'on')
            end
            
            ylabel(upplot,"$\chi$' (emu/mol)",'Interpreter','latex');
            xlabel(upplot,'Frequence (Hz)', 'Interpreter', 'latex');
            legend(upplot,'Location','Best');

            [Cswitch,~] = listdlg('ListString',{'1 contribution', '2 contributions','skip this entry'},'OKString','Enter','PromptString','Your choice ?','SelectionMode','single','ListSize',[160,50]);
            close('temp_data')

        if isempty(Cswitch)
            disp('Execution was interrupted')
            break
        end
        
        if Cswitch == 3
            if i==1
                lastparam1([1 3 5 7]) = param1([1 2 3 4]);
                lastparam2([1 3 5 7 9 11 13]) = param2([1 2 3 4 5 6 7]);
                figure('Name','Chisec','units','normalized','outerposition',[0.3 0.3 0.3 0.6]);
                axchisec = gca ; 
                figchisec = gcf ;        
                figure('Name','Chipr','units','normalized','outerposition',[0 0.3 0.3 0.6]);
                axchipr = gca ; 
                figchipr = gcf ;

            end
            i=i+1 ;
            continue
        end

        if Cswitch == 1 %1 contribution chosen
           fitchoicepr = fitpr1 ;
           fitchoicesec = fitsec1 ;

           choicec = 1 ;
           choice(i,1) = 1;

           obj1 = @(x) [(1-w).*((x(4)+x(1)*(1+((2*pi*xdata*x(2)).^(1-x(3)))*sin((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3)))))-ydatapr) ; w.*(((x(1)*((2*pi*xdata*x(2)).^(1-x(3)))*cos((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3)))))-ydatasec); lambda*(x(1)+x(2)+x(3)+x(4))];


        elseif Cswitch == 2 %2 contributions chosen
            
            fitchoicepr = fitpr2 ;
           fitchoicesec = fitsec2 ;

           choicec = 2;
           choice(i,1) = 2;
           
           obj2 = @(x) [(1-w).*((x(4)+x(1)*(1+((2*pi*xdata*x(2)).^(1-x(3)))*sin((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3)))) + x(4)+x(5)*(1+((2*pi*xdata*x(6)).^(1-x(7)))*sin((pi*x(7))/2))./(1+(2*((2*pi*xdata*x(6)).^(1-x(7)))*sin(pi*(x(7)/2)))+((2*pi*xdata*x(6)).^(2-2*x(7)))))-ydatapr) ; w.*(((x(1)*((2*pi*xdata*x(2)).^(1-x(3)))*cos((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3))))+(x(5)*((2*pi*xdata*x(6)).^(1-x(7)))*cos((pi*x(7))/2))./(1+(2*((2*pi*xdata*x(6)).^(1-x(7)))*sin(pi*(x(7)/2)))+((2*pi*xdata*x(6)).^(2-2*x(7)))))-ydatasec); lambda*(x(5)+x(6)+x(7)); lambda*(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7))];

        end


        %CALCULATE THE WEIGHTED FIT

                    std_dev = std_devstock((frst+(i-1)*nbpt):(frst-1+i*nbpt)); %standard deviation from measures

                    scaled_dev = 1 + ((std_dev-min(std_devstock))./(max(std_dev)-min(std_dev))).*(2) ;

                    W = 1./(scaled_dev.^2); %weights matrix of each point

                    

                   
                    %functions to be optimized with weights
                     wobj1 = @(x) [W.*((1-w).*((x(4)+x(1)*(1+((2*pi*xdata*x(2)).^(1-x(3)))*sin((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3)))))-ydatapr)) ; W.*(w.*(((x(1)*((2*pi*xdata*x(2)).^(1-x(3)))*cos((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3)))))-ydatasec)); lambda*(x(1)+x(2)+x(3)+x(4))];
                     wobj2 = @(x) [W.*((1-w).*((x(4)+x(1)*(1+((2*pi*xdata*x(2)).^(1-x(3)))*sin((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3)))) + x(4)+x(5)*(1+((2*pi*xdata*x(6)).^(1-x(7)))*sin((pi*x(7))/2))./(1+(2*((2*pi*xdata*x(6)).^(1-x(7)))*sin(pi*(x(7)/2)))+((2*pi*xdata*x(6)).^(2-2*x(7)))))-ydatapr)) ; W.*(w.*(((x(1)*((2*pi*xdata*x(2)).^(1-x(3)))*cos((pi*x(3))/2))./(1+(2*((2*pi*xdata*x(2)).^(1-x(3)))*sin(pi*(x(3)/2)))+((2*pi*xdata*x(2)).^(2-2*x(3))))+(x(5)*((2*pi*xdata*x(6)).^(1-x(7)))*cos((pi*x(7))/2))./(1+(2*((2*pi*xdata*x(6)).^(1-x(7)))*sin(pi*(x(7)/2)))+((2*pi*xdata*x(6)).^(2-2*x(7)))))-ydatasec)); lambda*(x(5)+x(6)+x(7)); lambda*(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7))];


                    %optimization with 2 contributions from initial parameters or parameters from last loop

                    if choicec == 2
                        if i==1 || only_one
                            [tempparam2w([1 3 5 7 9 11 13]),resnorm2w,residual2w,exitflag2w,~,~,jacobian2w] = lsqnonlin(wobj2,param2,lb2,ub2,optionsnonlin);
                        else
                            [tempparam2w([1 3 5 7 9 11 13]),resnorm2w,residual2w,exitflag2w,~,~,jacobian2w] = lsqnonlin(wobj2,lastparam2([1 3 5 7 9 11 13]),lb2,ub2,optionsnonlin);
                        end

                        %Optimization with 2 contributions from classic patterns defined earlier

                        [tempparam2wp2([1 3 5 7 9 11 13]),resnorm2wp2,residual2wp2,exitflag2wp2,~,~,jacobian2wp2] = lsqnonlin(wobj2,pattern2,lb2,ub2,optionsnonlin);
                        [tempparam2wp3([1 3 5 7 9 11 13]),resnorm2wp3,residual2wp3,exitflag2wp3,~,~,jacobian2wp3] = lsqnonlin(wobj2,pattern3,lb2,ub2,optionsnonlin);
                        [tempparam2wp4([1 3 5 7 9 11 13]),resnorm2wp4,residual2wp4,exitflag2wp4,~,~,jacobian2wp4] = lsqnonlin(wobj2,pattern4,lb2,ub2,optionsnonlin);
                        
                         %evaluating the best obtained fit (minimum resnorm)
                         [~,bestw] = min([resnorm2w resnorm2wp4 resnorm2wp2 resnorm2wp3 +Inf +Inf +Inf +Inf]) ;

                    end
                    %optimization with 1 contribution from initial parameters or parameters from last loop

                    if choicec==1
                        if i==1 || only_one
                            [tempparam1w([1 3 5 7]),resnorm1w,residual1w,exitflag1w,~,~,jacobian1w] = lsqnonlin(wobj1,param1,lb1,ub1,optionsnonlin);
                        else
                            [tempparam1w([1 3 5 7]),resnorm1w,residual1w,exitflag1w,~,~,jacobian1w] = lsqnonlin(wobj1,lastparam1([1 3 5 7]),lb1,ub1,optionsnonlin);
                        end

                        %Optimization with 1 contribution from classic patterns defined earlier


                          [tempparam1wp1([1 3 5 7 ]),resnorm1wp1,residual1wp1,exitflag1wp1,~,~,jacobian1wp1] = lsqnonlin(wobj1,pattern1,lb1,ub1,optionsnonlin);
                          [tempparam1wp5([1 3 5 7 ]),resnorm1wp5,residual1wp5,exitflag1wp5,~,~,jacobian1wp5] = lsqnonlin(wobj1,pattern5,lb1,ub1,optionsnonlin);
       
                          %uncomment here to add pattern6
                          % [tempparam1wp6([1 3 5 7 ]),resnorm1wp6,residual1wp6,exitflag1wp6,~,~,jacobian1wp6] = lsqnonlin(wobj1,pattern6,lb1,ub1,optionsnonlin);

                         %evaluating the best obtained fit (minimum resnorm)
                         [~,bestw] = min([+Inf +Inf +Inf +Inf resnorm1w resnorm1wp1 resnorm1wp5 +Inf]) ; 
                         %[~,bestw] = min([+Inf +Inf +Inf +Inf resnorm1w resnorm1wp1 resnorm1wp5 resnorm1wp6]) ;     %uncomment here to add pattern6
                    end

                   

                    %saving the best parameters

                    if bestw == 1  
                        
                        fitparamw([1 3 5 7 9 11 13]) = tempparam2w([1 3 5 7 9 11 13]) ;
                        resnormw = resnorm2w;
                        residualw = residual2w;
                        exitflagw = exitflag2w;
                        jacobianw = jacobian2w;
                        
                    elseif bestw == 2
                        
                        fitparamw([1 3 5 7 9 11 13]) = tempparam2wp4([1 3 5 7 9 11 13]) ;
                        resnormw = resnorm2wp4;
                        residualw = residual2wp4;
                        exitflagw = exitflag2wp4;
                        jacobianw = jacobian2wp4;
                    elseif bestw == 3
                        
                        fitparamw([1 3 5 7 9 11 13]) = tempparam2wp2([1 3 5 7 9 11 13]) ;
                        resnormw = resnorm2wp2;
                        residualw = residual2wp2;
                        exitflagw = exitflag2wp2;
                        jacobianw = jacobian2wp2;
                    elseif bestw == 4
                        fitparamw([1 3 5 7 9 11 13]) = tempparam2wp3([1 3 5 7 9 11 13]) ;
                        resnormw = resnorm2wp3;
                        residualw = residual2wp3;
                        exitflagw = exitflag2wp3;
                        jacobianw = jacobian2wp3;
                    elseif bestw == 5
                        fitparamw([1 3 5 7 ]) = tempparam1w([1 3 5 7 ]) ;
                        fitparamw([9 11 13]) = [NaN NaN NaN] ;
                        resnormw = resnorm1w;
                        residualw = residual1w;
                        exitflagw = exitflag1w;
                        jacobianw = jacobian1w;
                    elseif bestw == 6
                        fitparamw([1 3 5 7 ]) = tempparam1wp1([1 3 5 7 ]) ;
                        fitparamw([9 11 13]) = [NaN NaN NaN] ;
                        resnormw = resnorm1wp1;
                        residualw = residual1wp1;
                        exitflagw = exitflag1wp1;
                        jacobianw = jacobian1wp1;
                    elseif bestw == 7
                        fitparamw([1 3 5 7 ]) = tempparam1wp5([1 3 5 7 ]) ; 
                        fitparamw([9 11 13]) = [NaN NaN NaN] ;
                        resnormw = resnorm1wp5;
                        residualw = residual1wp5;
                        exitflagw = exitflag1wp5;
                        jacobianw = jacobian1wp5;
    %                 elseif bestw == 8                                          %uncomment here to add pattern6
    %                     fitparamw([1 3 5 7 ]) = tempparam1wp6([1 3 5 7 ]) ;
    %                     fitparamw([9 11 13]) = [NaN NaN NaN] ;
    %                     resnormw = resnorm1wp6;
    %                     residualw = residual1wp6;
    %                     exitflagw = exitflag1wp6;
    %                     jacobianw = jacobian1wp6;
                    end 

                    tempparamwstock(i,:) = fitparamw([1 3 5 7 9 11 13]);




        %CALCULATE THE UNWEIGHTED FIT
                    if choicec ==2
                        if i==1 || only_one
                            [tempparam2uw([1 3 5 7 9 11 13]),resnorm2uw,residual2uw,exitflag2uw,~,~,jacobian2uw] =  lsqnonlin(obj2,param2,lb2,ub2,optionsnonlin);
                        else
                            [tempparam2uw([1 3 5 7 9 11 13]),resnorm2uw,residual2uw,exitflag2uw,~,~,jacobian2uw] = lsqnonlin(obj2,lastparam2([1 3 5 7 9 11 13]),lb2,ub2,optionsnonlin);
                        end
                        
                        [tempparam2uwp2([1 3 5 7 9 11 13]),resnorm2uwp2,residual2uwp2,exitflag2uwp2,~,~,jacobian2uwp2] = lsqnonlin(obj2,pattern2,lb2,ub2,optionsnonlin);
                        [tempparam2uwp3([1 3 5 7 9 11 13]),resnorm2uwp3,residual2uwp3,exitflag2uwp3,~,~,jacobian2uwp3] = lsqnonlin(obj2,pattern3,lb2,ub2,optionsnonlin);
                        [tempparam2uwp4([1 3 5 7 9 11 13]),resnorm2uwp4,residual2uwp4,exitflag2uwp4,~,~,jacobian2uwp4] = lsqnonlin(wobj2,pattern4,lb2,ub2,optionsnonlin);
                        
                        [~,bestuw] = min([resnorm2w resnorm2uwp4 resnorm2wp2 resnorm2wp3 +Inf +Inf +Inf +Inf]) ;

                    end
                    
                    
                    if choicec == 1
                        if i==1 || only_one
                            [tempparam1uw([1 3 5 7]),resnorm1uw,residual1uw,exitflag1uw,~,~,jacobian1uw] = lsqnonlin(obj1,param1,lb1,ub1,optionsnonlin);
                        else
                            [tempparam1uw([1 3 5 7]),resnorm1uw,residual1uw,exitflag1uw,~,~,jacobian1uw] = lsqnonlin(obj1,lastparam1([1 3 5 7]),lb1,ub1,optionsnonlin);
                        end
                    

                      [tempparam1uwp1([1 3 5 7 ]),resnorm1uwp1,residual1uwp1,exitflag1uwp1,~,~,jacobian1uwp1] = lsqnonlin(obj1,pattern1,lb1,ub1,optionsnonlin);
                      [tempparam1uwp5([1 3 5 7 ]),resnorm1uwp5,residual1uwp5,exitflag1uwp5,~,~,jacobian1uwp5] = lsqnonlin(obj1,pattern5,lb1,ub1,optionsnonlin);
                      
                      %uncomment here to add pattern6
                      %[tempparam1uwp6([1 3 5 7 ]),resnorm1uwp6,residual1uwp6,exitflag1uwp6,~,~,jacobian1uwp6] = lsqnonlin(obj1,pattern6,lb1,ub1,optionsnonlin);
                      
                      [~,bestuw] = min([+Inf +Inf +Inf +Inf resnorm1uw resnorm1uwp1 resnorm1uwp5 +Inf]) ;
                      %[~,bestuw] = min([+Inf +Inf +Inf +Inf resnorm1uw resnorm1uwp1 resnorm1uwp5 +Inf]) ;  %uncomment here to add pattern6
                    end

                    %choosing 

                    if bestuw == 1  
                        fitparamuw([1 3 5 7 9 11 13]) = tempparam2uw([1 3 5 7 9 11 13]) ;
                        resnormuw = resnorm2uw;
                        residualuw = residual2uw;
                        exitflaguw = exitflag2uw;
                        jacobianuw = jacobian2uw;
                    elseif bestuw == 2
                        fitparamuw([1 3 5 7 9 11 13]) = tempparam2uwp4([1 3 5 7 9 11 13]) ;
                        resnormuw = resnorm2uwp4;
                        residualuw = residual2uwp4;
                        exitflaguw = exitflag2uwp4;
                        jacobianuw = jacobian2uwp4;
                    elseif bestuw == 3
                        fitparamuw([1 3 5 7 9 11 13]) = tempparam2uwp2([1 3 5 7 9 11 13]) ;
                        resnormuw = resnorm2uwp2;
                        residualuw = residual2uwp2;
                        exitflaguw = exitflag2uwp2;
                        jacobianuw = jacobian2uwp2;
                    elseif bestuw == 4
                        fitparamuw([1 3 5 7 9 11 13]) = tempparam2uwp3([1 3 5 7 9 11 13]) ;
                        resnormuw = resnorm2uwp3;
                        residualuw = residual2uwp3;
                        exitflaguw = exitflag2uwp3;
                        jacobianuw = jacobian2uwp3;
                    elseif bestuw == 5
                        fitparamuw([1 3 5 7 ]) = tempparam1uw([1 3 5 7 ]) ;
                        fitparamuw([9 11 13]) = [NaN NaN NaN] ;
                        resnormuw = resnorm1uw;
                        residualuw = residual1uw;
                        exitflaguw = exitflag1uw;
                        jacobianuw = jacobian1uw;
                    elseif bestuw == 6
                        fitparamuw([1 3 5 7 ]) = tempparam1uwp1([1 3 5 7 ]) ;
                        fitparamuw([9 11 13]) = [NaN NaN NaN] ;
                        resnormuw = resnorm1uwp1;
                        residualuw = residual1uwp1;
                        exitflaguw = exitflag1uwp1;
                        jacobianuw = jacobian1uwp1;
                    elseif bestuw == 7
                        fitparamuw([1 3 5 7 ]) = tempparam1uwp5([1 3 5 7 ]) ;
                        fitparamuw([9 11 13]) = [NaN NaN NaN] ;
                        resnormuw = resnorm1uwp5;
                        residualuw = residual1uwp5;
                        exitflaguw = exitflag1uwp5;
                        jacobianuw = jacobian1uwp5;
    %                 elseif bestuw == 8                                             %uncomment here to add pattern6
    %                     fitparamuw([1 3 5 7 ]) = tempparam1uwp6([1 3 5 7 ]) ;
    %                     fitparamuw([9 11 13]) = [NaN NaN NaN] ;
    %                     resnormuw = resnorm1uwp6;
    %                     residualuw = residual1uwp6;
    %                     exitflaguw = exitflag1uwp6;
    %                     jacobianuw = jacobian1uwp6;
                    end 


                    tempparamuwstock(i,:) = fitparamuw([1 3 5 7 9 11 13]);


        %USER CHOICE ON WEIGHTED OR UNWEIGHTED FIT


            figure('Name','temp_fig','units','normalized','outerposition',[0.6 0 0.4 1]) ;
            upplot = subplot(2,1,1) ;
            downplot = subplot(2,1,2) ;
            %tempaxes = axes;

           
            if contains(TorH,'T')
                semilogx(downplot,xdata,ydatasec,'o','DisplayName', 'Data ' + string(rTemp)+'K','color','k');
                hold(downplot,'on')

            else
                semilogx(downplot,xdata,ydatasec,'o','DisplayName', 'Data ' + string(rField)+'Oe','color','k');
                hold(downplot,'on')
            end
            semilogx(downplot,xdata,fitchoicesec(fitparamw([1 3 5 7 9 11 13]),xdata),'b-','color','g','DisplayName','Weighted fit'); %plot fit weighted
            hold(downplot,'on')
            semilogx(downplot,xdata,fitchoicesec(fitparamuw([1 3 5 7 9 11 13]),xdata),'b-','color','r','DisplayName','Not Weighted fit'); %plot fit unweighted
            hold(downplot,'on')
            ylabel(downplot,'$\chi$" (emu/mol)','Interpreter','latex');
            xlabel(downplot,'Frequence (Hz)', 'Interpreter', 'latex');
            legend(downplot,'Location','Best');


            if contains(TorH,'T')
                semilogx(upplot,xdata,ydatapr,'o','DisplayName',"Data " + string(rTemp)+"K",'color','k');
                hold(upplot,'on')
            else
                semilogx(upplot,xdata,ydatapr,'o','DisplayName', "Data " + string(rField)+'Oe','color','k');
                hold(upplot,'on')
            end
            semilogx(upplot,xdata,fitchoicepr(fitparamw([1 3 5 7 9 11 13]),xdata),'b-','color','g','DisplayName','Weighted fit'); %plot fit weighted
            hold(upplot,'on')
            semilogx(upplot,xdata,fitchoicepr(fitparamuw([1 3 5 7 9 11 13]),xdata),'b-','color','r','DisplayName','Not Weighted fit'); %plot fit unweighted
            hold(upplot,'on')
            ylabel(upplot,"$\chi$' (emu/mol)",'Interpreter','latex');
            xlabel(upplot,'Frequence (Hz)', 'Interpreter', 'latex');
            legend(upplot,'Location','Best');







            [Wswitch,~] = listdlg('ListString',{'Weighted', 'Not weighted'},'OKString','Enter','PromptString','Your choice ?','SelectionMode','single','ListSize',[160,50]);
            close('temp_fig')

        if isempty(Wswitch)
            disp('Execution was interrupted')
            break
        end
        
        if Wswitch == 1 %WEIGHTED CHOSEN
            fitparam([1 3 5 7 9 11 13]) = fitparamw([1 3 5 7 9 11 13]);

            residual = residualw;

            resnormstock(i,:) = resnormw ;

            
            if choicec == 2
                 for h = 1:7
                    for n = 1:nbpt
                         jacobian(n,h) = (jacobianw(n,h)+jacobianw(nbpt+n,h))/2 ; %jacobian = average of chi' and chi"
                    end
                 end
            elseif choicec == 1
                for h = 1:4
                    for n = 1:nbpt
                         jacobian(n,h) = (jacobianw(n,h)+jacobianw(nbpt+n,h))/2 ;
                    end
                end
            end
            
            

        elseif Wswitch == 2 %NOT WEIGHTED CHOSEN
            fitparam([1 3 5 7 9 11 13]) = fitparamuw([1 3 5 7 9 11 13]);

            residual = residualuw;

            resnormstock(i,:) = resnormuw ;

            if choicec == 2
                 for h = 1:7
                    for n = 1:nbpt
                         jacobian(n,h) = (jacobianuw(n,h)+jacobianuw(nbpt+n,h))/2 ;
                    end
                 end
            elseif choicec == 1
                for h = 1:4
                    for n = 1:nbpt
                         jacobian(n,h) = (jacobianuw(n,h)+jacobianuw(nbpt+n,h))/2 ;
                    end
                end
            end
        end
        
        
        fitparamstock(i,[1 3 5 7 9 11 13]) = fitparam([1 3 5 7 9 11 13]) ;
        
        for n = 1:nbpt
         residualstock(n,i) = (residual(n)+residual(nbpt+n))/2 ; %residual = average of chi' and chi" residual
        end
        
      
        
        %For errors display

        if choice(i,1) == 2 
            newjac = reshape(jacobian,nbpt,7) ;
            J{i} = full(newjac);

            errors{i} = nlparci(fitparam([1 3 5 7 9 11 13]),residualstock(:,i),'jacobian',J{i}); %95% confidence intervals for d1

            for p  = 1:7
                fitparamstock(i,2*p)= (errors{i}(p,2)-errors{i}(p,1))/2 ; %(highest value - lowest value) / 2 to have a +/- display
            end

        elseif choice(i,1) == 1
            newjac = reshape(jacobian,nbpt,4) ;
            J{i} = full(newjac);

            errors{i} = nlparci(fitparam([1 3 5 7]),residualstock(:,i),'jacobian',J{i}); %95% confidence intervals for d1

            for p  = 1:4
                fitparamstock(i,2*p) = (errors{i}(p,2)-errors{i}(p,1))/2 ;
                fitparamstock(i,2*p+8) = NaN ;
            end
        end



         %Add three columns with temperatures and 1/T and H


        fitparamstock(i,17)= Temp;
        fitparamstock(i,18)=1/fitparamstock(i,17);
        fitparamstock(i,19)=Field;


         %%display parameters  %2.6
            if i==1

                figparam = uifigure('Name','Parameters','units','normalized','Position',[0.001 0.01 0.6 0.3]); %2.6
                uit = uitable(figparam,'units','normalized','Position',[0 0 1 1],'ColumnName',Column_names);
            end
            
            uit.Data(i,:)=fitparamstock(i,:);

        %%%%PLOTS%%%%


        if (i==1 || only_one) && contains(plot_option,'single')
            figure('Name','Chisec','units','normalized','outerposition',[0.3 0.4 0.3 0.58]);
            axchisec = gca ;
            figchisec = gcf ;        
            figure('Name','Chipr','units','normalized','outerposition',[0 0.4 0.3 0.58]);
            axchipr = gca ; 
            figchipr = gcf ;
        end

        %plot chi"

        if contains(plot_option,'separated') 
            if contains(TorH,'T')
                figure('Name','Chisec_'+string(rTemp)+'K','units','normalized','outerposition',[0.3 0.3 0.3 0.6]);
            elseif contains(TorH,'H')
                figure('Name','Chisec_'+string(rField)+'Oe','units','normalized','outerposition',[0 0.3 0.3 0.6]);
            end
            axchisec = gca ;
        end
        
        if contains(TorH,'T')
                pointssec(i) = semilogx(axchisec,xdata,ydatasec,'o','DisplayName', string(rTemp)+'K','color',newcolors(i,:));
                hold(axchisec,'on')
        else
                pointssec(i) = semilogx(axchisec,xdata,ydatasec,'o','DisplayName', string(rField)+'Oe','color',newcolors(i,:));
                hold(axchisec,'on')
        end
        
        yfitsec = fitchoicesec(fitparamstock(i,[1 3 5 7 9 11 13]),xdata);

        curvesec = semilogx(axchisec,xdata,yfitsec,'b-','color',newcolors(i,:));
        %'DisplayName', 'fit '+ string(rTemp)+'K',
        hold(axchisec,'on')
        xtsec = xticks(axchisec) ;
        xticklabels(axchisec,num2cell(xtsec));
        ylabel(axchisec,'$\chi$" (emu/mol)','Interpreter','latex');
        xlabel(axchisec,'Frequence (Hz)', 'Interpreter', 'latex');

        title(axchisec,'$\chi$"','Interpreter','latex');

        %plot chi '    

        if contains(plot_option,'separated') 
            if contains(TorH,'T')
                figure('Name','Chipr_'+string(rTemp)+'K','units','normalized','outerposition',[0 0.3 0.3 0.6]);
            elseif contains(TorH,'H')
                figure('Name','Chipr_'+string(rField)+'Oe','units','normalized','outerposition',[0 0.3 0.3 0.6]);
            end
            axchipr = gca ; 
        end

        if contains(TorH,'T')
                pointspr(i) = semilogx(axchipr,xdata,ydatapr,'o','DisplayName', string(rTemp)+'K','color',newcolors(i,:));
                hold(axchipr,'on')
                legendtemp{i} = string(rTemp)+' K';
        elseif contains(TorH,'H')
                pointspr(i) = semilogx(axchipr,xdata,ydatapr,'o','DisplayName', string(rField)+'Oe','color',newcolors(i,:));
                hold(axchipr,'on')
                legendfield{i} = string(rField)+' Oe';  
        end
        
        %creating a file with the x and y values of the fits
        
        yfitpr = fitchoicepr(fitparamstock(i,[1 3 5 7 9 11 13]),xdata);

        curvepr = semilogx(axchipr,xdata,yfitpr,'b-','color',newcolors(i,:));
        hold(axchipr,'on')
        xtpr = xticks(axchipr) ;
        xticklabels(axchipr,num2cell(xtpr));
        ylabel(axchipr,"$\chi$' (emu/mol)",'Interpreter','latex');
        xlabel(axchipr,'Frequence (Hz)', 'Interpreter', 'latex');
        
           

        title(axchipr,"$\chi$'",'Interpreter','latex');



        if contains(plot_option,'separated') && contains(TorH,'T')
            legend(pointspr(i),string(rTemp)+' K','Location','Best'); 
            legend(pointssec(i),string(rTemp)+' K','Location','Best'); 
        elseif contains(plot_option,'separated') && contains(TorH,'H')
            legend(pointspr(i), string(rField)+' Oe','Location','Best');
            legend(pointssec(i),string(rField)+' Oe','Location','Best'); 
        end

        %File with frequence, chi'_data, chi'_fit , chi"_data, chi"_fit, H and T 
        
        outputfile(((i-1)*nbpt+1:i*nbpt),1) = xdata;
        outputfile(((i-1)*nbpt+1:i*nbpt),2) = ydatapr;
        outputfile(((i-1)*nbpt+1:i*nbpt),3) = yfitpr;
        outputfile(((i-1)*nbpt+1:i*nbpt),4) = ydatasec;
        outputfile(((i-1)*nbpt+1:i*nbpt),5) = yfitsec;
        outputfile(((i-1)*nbpt+1:i*nbpt),6) = Temps;        
        outputfile(((i-1)*nbpt+1:i*nbpt),7) = Fields;
        
        
        %Now save the selected parameters to serve as initial parameters in next iteration

        if choice(i,1)==1
            lastparam1([1 3 5 7]) = fitparamstock(i,[1 3 5 7]) ;
            lastparam2([1 3 5 7]) = fitparamstock(i,[1 3 5 7]) ;
            lastparam2([9 11 13]) = param2([5 6 7]) ;
        end

        if choice(i,1) == 2
            lastparam2([1 3 5 7 9 11 13]) = fitparamstock(i,[1 3 5 7 9 11 13]) ;
            if fitparamstock(i,1)>fitparamstock(i,9)
                lastparam1([1 3 5 7])= fitparamstock(i,[1 3 5 7]);
            else 
                lastparam1([1 3 5 7])= fitparamstock(i,[9 11 13 7]);
            end 
        end   


        %Fix chiS2 as ChiS1 (they were already fitted together, just copying the value)

        fitparamstock(i,[15 16]) = fitparamstock(i,[7 8]);

        
        %Adding columns with ln(tau) and its error for both contributions
        fitparamstock(i,20) = log(fitparamstock(i,3)); %ln(tau)
        fitparamstock(i,21) = fitparamstock(i,4)/fitparamstock(i,3) ; %err of ln(tau) = err(tau)/tau
        
        if choicec == 2
            fitparamstock(i,22) = log(fitparamstock(i,11)); %ln(tau)
            fitparamstock(i,23) = fitparamstock(i,12)/fitparamstock(i,11) ; %err of ln(tau) = err(tau)/tau
        elseif choicec == 1
            fitparamstock(i,22) = NaN; 
            fitparamstock(i,23) = NaN ; 
        
        end
        
        
        
        i=i+1;
    end
  
%% SAVES
    
%Now clear the ignored sets and saving files

Column_names = {'d1','errd1','t1','errt1','a1','erra1','chiS1','errchiS1','d2','errd2','t2','errt2','a2','erra2','chiS2','errchiS2','T','1/T','H','ln(t1)','err(ln(t1))','ln(t2)','err(ln(t2))'} ;

Out_columns = {'Frequence','Chi''_data','Chi''_fit','Chi''''_data','Chi''''_fit','T','H'} ;

if contains(plot_option,'single')

    seclegend = gobjects(1,0);
    prlegend = gobjects(1,0);

    for j = 1:nbfit
        if isgraphics(pointssec(j))
            seclegend(end+1) = pointssec(j);
            prlegend(end+1) = pointspr(j);
        end
    end    


    if contains(TorH,'T')
        legendtemp = legendtemp(~cellfun('isempty',legendtemp));
        legend(seclegend,legendtemp, 'NumColumns',2,'Location','Best');
        legend(prlegend,legendtemp, 'NumColumns',2,'Location','Best');
    elseif contains(TorH,'H')
        legendfield = legendfield(~cellfun('isempty',legendfield));
        legend(seclegend,legendfield, 'NumColumns',2,'Location','Best');
        legend(prlegend,legendfield, 'NumColumns',2,'Location','Best');
    end
    
    if only_one
        savefile = questdlg("Do you want to save this fit ?");
    else
        savefile = questdlg("Do you want to save these fits ?");
    end

    if contains(savefile,'No') || contains(savefile,'Cancel')
        disp('No file was saved. To save the fits please run the last section of the code %%SAVES (ctrl+enter).') %2.5
    
    elseif contains(savefile,'Yes') && ~only_one
        set(figchipr, 'units','normalized','outerposition',[0 0 1 1])
        set(figchisec,'units','normalized','outerposition',[0 0 1 1])
        status = mkdir(strcat(filepath,'\',filename)) ;
        exportgraphics(axchipr,strcat(filepath,'\',filename,'\',filename,'_chiprime.png'));
        exportgraphics(axchisec,strcat(filepath,'\',filename,'\',filename,'_chisecond.png'));
        params = fitparamstock;
        params = params(any(params~=0,2), : );
        table = array2table(params,'VariableNames',Column_names) ;
        out = outputfile ;
        out = out(any(out~=0,2),:);
        tableout = array2table(out,'VariableNames',Out_columns);        
        writetable(tableout,strcat(filepath,'\',filename,'\',filename,'_Curves.txt'),'Delimiter','\t') ;
        writetable(table,strcat(filepath,'\',filename,'\',filename,'_Parameters.txt'),'Delimiter','\t') ;
        disp(['The files were saved successfully in directory ',filepath,'\',filename,'\']);
        
    elseif contains(savefile,'Yes') && only_one
        set(figchipr, 'units','normalized','outerposition',[0 0 1 1])
        set(figchisec,'units','normalized','outerposition',[0 0 1 1])
        status = mkdir(strcat(filepath,'\',filename)) ;
        params = fitparamstock;
        params = params(any(params~=0,2), : );
        table = array2table(params,'VariableNames',Column_names) ;
        out = outputfile ;
        out = out(any(out~=0,2),:);
        tableout = array2table(out,'VariableNames',Out_columns);        
        if contains(TorH,'T')
            exportgraphics(axchipr,strcat(filepath,'\',filename,'\',filename,'_chiprime',string(rTemp),'K.png'));
            exportgraphics(axchisec,strcat(filepath,'\',filename,'\',filename,'_chisecond',string(rTemp),'K.png'));
            writetable(tableout,strcat(filepath,'\',filename,'\',filename,'_Curves',string(rTemp),'K.txt'),'Delimiter','\t') ;
            writetable(table,strcat(filepath,'\',filename,'\',filename,'_Parameters',string(rTemp),'K.txt'),'Delimiter','\t') ;
        elseif contains(TorH,'H')
            exportgraphics(axchipr,strcat(filepath,'\',filename,'\',filename,'_chiprime',string(rField),'Oe.png'));
            exportgraphics(axchisec,strcat(filepath,'\',filename,'\',filename,'_chisecond',string(rField),'Oe.png'));
            writetable(tableout,strcat(filepath,'\',filename,'\',filename,'_Curves',string(rField),'Oe.txt'),'Delimiter','\t') ;
            writetable(table,strcat(filepath,'\',filename,'\',filename,'_Parameters',string(rField),'Oe.txt'),'Delimiter','\t') ;
        end
        disp(['The files were saved successfully in directory ',filepath,'\',filename,'\']);
    end
    

end



if contains(plot_option,'separated')
    savefile = questdlg("Do you want to save these fits ?");

    if contains(savefile,'Yes')
        status = mkdir(strcat(filepath,'\',filename)) ;
        params = fitparamstock;
        params = params(any(params~=0,2), : );
        table = array2table(params,'VariableNames',Column_names);
        out = outputfile ;
        out = out(any(out~=0,2),:);
        tableout = array2table(out,'VariableNames',Out_columns);        
        writetable(tableout,strcat(filepath,'\',filename,'\',filename,'_Curves.txt'),'Delimiter','\t') ;        
        writetable(table,strcat(filepath,'\',filename,'\',filename,'_Parameters.txt'),'Delimiter','\t') ;
        FigList = findobj('Type', 'figure');
        for iFig = 1:length(FigList)
          FigHandle = FigList(iFig);
          set(FigHandle, 'units','normalized','outerposition',[0 0 1 1])
          FigName   = get(FigHandle, 'Name');
          exportgraphics(FigHandle,strcat(filepath,'\',filename,'\',filename,'_',FigName,'.png'));       
        end
        disp(['The files were saved successfully in directory ',filepath,'\',filename,'\']);
    end
end


%%%% Thank you Ilyes for your precious help %%%%
