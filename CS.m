%%Developed by Alvaro Martinez Guimera - 2021%%
clear
clc

%Shuffle random number generator seeding
rng('shuffle')


%Parameter range
ParamVals=linspace(0.01, 0.1 , 10);

%Timecourse settings
Timecourse_Duration=400;

%Species abundance
SpeciesAbundanceScale=50;

%Create SimBiology model (once model is made this part of the code (line 11 to 63) can be commenetd out and the model imported by unhashing line 70)

%
%Define model object
modelObj = sbiomodel ('Model_1'); %Define model name

%Create compartment
compObj = addcompartment(modelObj, 'cell');

%Create species
speciesObj1 = addspecies (compObj, 'A');
speciesObj2 = addspecies (compObj, 'Aact');
speciesObj3 = addspecies (compObj, 'B');
speciesObj4 = addspecies (compObj, 'Bact');
speciesObj5 = addspecies (compObj, 'C');
speciesObj6 = addspecies (compObj, 'Cact');
speciesObj7 = addspecies (compObj, 'NegReg');
speciesObj8 = addspecies (compObj, 'Stimulus');
speciesObj9 = addspecies (compObj, 'CS');
speciesObj10 = addspecies (compObj, 'Nil');
speciesObj11 = addspecies (compObj, 'Ainhib');
speciesObj12 = addspecies (compObj, 'Bihhib');
speciesObj13 = addspecies (compObj, 'Cinhib');

%Set initial abundances
set (speciesObj1, 'InitialAmount',SpeciesAbundanceScale);
set (speciesObj2, 'InitialAmount',0);
set (speciesObj3, 'InitialAmount',SpeciesAbundanceScale);
set (speciesObj4, 'InitialAmount',0); 
set (speciesObj5, 'InitialAmount',SpeciesAbundanceScale);
set (speciesObj6, 'InitialAmount',0); 
set (speciesObj7, 'InitialAmount',0); 
set (speciesObj8, 'InitialAmount',0); 
set (speciesObj9, 'InitialAmount',1); 
set (speciesObj10, 'InitialAmount',1); 
set (speciesObj11, 'InitialAmount',0); 
set (speciesObj12, 'InitialAmount',0); 
set (speciesObj13, 'InitialAmount',0);  

% Add a kinetic parameters
k1=ParamVals(randi(numel(ParamVals)));
k2=ParamVals(randi(numel(ParamVals)));
k3=ParamVals(randi(numel(ParamVals)));
Kinhib=ParamVals(randi(numel(ParamVals)));
Kregen=(ParamVals(randi(numel(ParamVals)))/2);
Kcs1=0;
Kcs2=0;
Kcs3=0;
Ksynth=(min(ParamVals))/60;
Kdeg=Ksynth/SpeciesAbundanceScale;
Knr=(ParamVals(randi(numel(ParamVals))))/2;
Knrd=(ParamVals(randi(numel(ParamVals))))/2;

parameter1 = addparameter(modelObj, 'k1',k1 , 'ConstantValue', false);
parameter2 = addparameter(modelObj, 'k1r', 2*k1, 'ConstantValue', false);
parameter3 = addparameter(modelObj, 'k2', k2, 'ConstantValue', false);
parameter4 = addparameter(modelObj, 'k2r', 2*k2, 'ConstantValue', false);
parameter5 = addparameter(modelObj, 'k3', k3, 'ConstantValue', false);
parameter6 = addparameter(modelObj, 'k3r', 2*k3, 'ConstantValue', false);
parameter7 = addparameter(modelObj, 'Kinhib', Kinhib, 'ConstantValue', false);
parameter8 = addparameter(modelObj, 'Kcs1', Kcs1, 'ConstantValue', false);
parameter9 = addparameter(modelObj, 'Kcs2', Kcs2, 'ConstantValue', false);
parameter10 = addparameter(modelObj, 'Kcs3', Kcs3, 'ConstantValue', false);
parameter11 = addparameter(modelObj, 'Ksynth', Ksynth, 'ConstantValue', false);
parameter12 = addparameter(modelObj, 'Kdeg', Kdeg, 'ConstantValue', false);
parameter13 = addparameter(modelObj, 'Knr', Knr, 'ConstantValue', false);
parameter14 = addparameter(modelObj, 'Knrd', Knrd, 'ConstantValue', false);
parameter15 = addparameter(modelObj, 'Kregen', Kregen, 'ConstantValue', false);

%Add reactions
reaction1 = addreaction(modelObj, 'A + Stimulus -> Aact + Stimulus'); 
reaction2 = addreaction(modelObj, 'Aact -> A'); 
reaction3 = addreaction(modelObj, 'Aact + B -> Bact + Aact'); 
reaction4 = addreaction(modelObj, 'Bact -> B'); 
reaction5 = addreaction(modelObj, 'Bact + C -> Cact + Bact'); 
reaction6 = addreaction(modelObj, 'Cact -> C'); 
reaction7 = addreaction(modelObj, 'Nil -> A + Nil'); 
reaction8 = addreaction(modelObj, 'Nil -> B + Nil');
reaction9 = addreaction(modelObj, 'Nil -> C + Nil');
reaction10 = addreaction(modelObj, 'A + Nil -> Nil'); 
reaction11 = addreaction(modelObj, 'B + Nil -> Nil');
reaction12 = addreaction(modelObj, 'C + Nil -> Nil');
reaction13 = addreaction(modelObj, 'NegReg + Aact -> Ainhib + NegReg'); %Model 1
%reaction13= addreaction(modelObj, 'NegReg + Bact -> Binhib + NegReg'); %Model 2
%reaction13 = addreaction(modelObj, 'NegReg + Cact -> Cinhib + NegReg'); %Model 3
reaction14 = addreaction(modelObj, 'Cact -> Cact + NegReg'); 
reaction15 = addreaction(modelObj, 'NegReg + Nil -> Nil'); 
reaction16 = addreaction(modelObj, 'CS + A -> CS + Aact'); %set to zero?
reaction17 = addreaction(modelObj, 'CS + B -> CS + Bact'); %set to zero?
reaction18 = addreaction(modelObj, 'CS + C -> CS + Cact'); %set to zero?
reaction19 = addreaction(modelObj, 'Ainhib -> A'); %regeneration
reaction20 = addreaction(modelObj, 'Binhib  -> B'); %regeneration
reaction21 = addreaction(modelObj, 'Cinhib  -> C'); %regeneration
%reaction22 = addreaction(modelObj, 'NegReg + Cact -> Cinhib + NegReg'); %Model 3

%Add kinetic Laws
kineticLaw1 = addkineticlaw(reaction1,'MassAction');
kineticLaw2 = addkineticlaw(reaction2,'MassAction');
kineticLaw3 = addkineticlaw(reaction3,'MassAction');
kineticLaw4 = addkineticlaw(reaction4,'MassAction');
kineticLaw5 = addkineticlaw(reaction5,'MassAction');
kineticLaw6 = addkineticlaw(reaction6,'MassAction');
kineticLaw7 = addkineticlaw(reaction7,'MassAction');
kineticLaw8 = addkineticlaw(reaction8,'MassAction');
kineticLaw9 = addkineticlaw(reaction9,'MassAction');
kineticLaw10 = addkineticlaw(reaction10,'MassAction');
kineticLaw11 = addkineticlaw(reaction11,'MassAction');
kineticLaw12 = addkineticlaw(reaction12,'MassAction');
kineticLaw13 = addkineticlaw(reaction13,'MassAction');
kineticLaw14 = addkineticlaw(reaction14,'MassAction');
kineticLaw15 = addkineticlaw(reaction15,'MassAction');
kineticLaw16 = addkineticlaw(reaction16,'MassAction');
kineticLaw17 = addkineticlaw(reaction17,'MassAction');
kineticLaw18 = addkineticlaw(reaction18,'MassAction');
kineticLaw19 = addkineticlaw(reaction19,'MassAction');
kineticLaw20 = addkineticlaw(reaction20,'MassAction');
kineticLaw21 = addkineticlaw(reaction21,'MassAction');
%kineticLaw22 = addkineticlaw(reaction22,'MassAction');

%Provide parameter names
kineticLaw1.ParameterVariableNames = 'k1';  %Through the names the kinetic laws can then find the appropriate parameter (scoped to the model)
kineticLaw2.ParameterVariableNames = 'k1r';
kineticLaw3.ParameterVariableNames = 'k2';
kineticLaw4.ParameterVariableNames = 'k2r';
kineticLaw5.ParameterVariableNames = 'k3';
kineticLaw6.ParameterVariableNames = 'k3r';
kineticLaw7.ParameterVariableNames = 'Ksynth';
kineticLaw8.ParameterVariableNames = 'Ksynth';
kineticLaw9.ParameterVariableNames = 'Ksynth';
kineticLaw10.ParameterVariableNames = 'Kdeg';
kineticLaw11.ParameterVariableNames = 'Kdeg';
kineticLaw12.ParameterVariableNames = 'Kdeg';
kineticLaw13.ParameterVariableNames = 'Kinhib';
kineticLaw14.ParameterVariableNames = 'Knr';
kineticLaw15.ParameterVariableNames = 'Knrd';
kineticLaw16.ParameterVariableNames = 'Kcs1';
kineticLaw17.ParameterVariableNames = 'Kcs2';
kineticLaw18.ParameterVariableNames = 'Kcs3';
kineticLaw19.ParameterVariableNames = 'Kregen';
kineticLaw20.ParameterVariableNames = 'Kregen';
kineticLaw21.ParameterVariableNames = 'Kregen';
%kineticLaw22.ParameterVariableNames = 'Kinhib';

%Add stimulus event into simulation	
event1 = addevent(modelObj,'time>=200','Stimulus = 20'); 
event2 = addevent(modelObj,'time>=220','Stimulus = 0');
 
%Save Model as
sbiosaveproject('Model_1', 'modelObj') %Simbiology project
%sbmlexport(modelObj, 'Model1')  %SBML file

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%If model has already been created and saved then comment out all previous lines
%uncomment model import line 65 below.
sbioloadproject('Model_1')

%Number of parameter permutations
NumParamPerm=11000;

%Number of parameters in the model
NumOfParams=numel(modelObj.Parameters);

%Number of events in the model
NumOfEvents=numel(modelObj.Events);

Results=[];

for Parameter_Permutations= 1:NumParamPerm

	%Sample new parameter values Parameters
	k1=ParamVals(randi(numel(ParamVals)));
	k2=ParamVals(randi(numel(ParamVals)));
	k3=ParamVals(randi(numel(ParamVals)));
	Kinhib=ParamVals(randi(numel(ParamVals)));
    Kregen=(ParamVals(randi(numel(ParamVals)))/2);
    Knr=(ParamVals(randi(numel(ParamVals))))/2;
    Knrd=(ParamVals(randi(numel(ParamVals))))/2;
    
    %Set new parameter values
	modelObj.Parameters(1).value=k1;
	modelObj.Parameters(2).value=(k1*2);
	modelObj.Parameters(3).value=k2;
	modelObj.Parameters(4).value=(k2*2);
	modelObj.Parameters(5).value=k3;
	modelObj.Parameters(6).value=(k3*2);
	modelObj.Parameters(7).value=Kinhib;
    modelObj.Parameters(13).value=Knr;
    modelObj.Parameters(14).value=Knrd;
    modelObj.Parameters(15).value=Kregen;

	%In-loop Data
	Control_val=[];
	CS_val=[];
	Diff=[];
    Simulation_Output=[];
    Control_Output=[];

	%Perform a 'Control' deterministic simulation to compare subsequent analysis output to

	%Perform deterministic simulation
	cs = getconfigset(modelObj,'active');
	cs.SolverType = 'ode45';  %ODE solver - ODE45 non-stiff, ode45 & ODE23 = stiff
	cs.SolverOptions.AbsoluteTolerance= 1.0e-12;
	cs.SolverOptions.RelativeTolerance= 1.0e-6;
	cs.StopTime = Timecourse_Duration; %Simulation stop time
	%cs.SolverOptions.LogDecimation = 200;   %how frequently you want to record the output of a stochastic simulation (ex. every 200 ime units)
	cs.CompileOptions.UnitConversion = false;  %No unit conversion
	Simulation_Output=sbiosimulate(modelObj); %Simulate model
	%sbioplot(Simulation_Output);  %plot simulation output

	%Extract control data 
	Control_Output=Simulation_Output.Data(:,1:end-NumOfParams); %Columns correspond to species in the order they appear in the modelObj. The last column indices corresponding to the
	%number of parameters are not included since when parameters are defined as non-constant (so they can be changed by an Event) then they are plotted as variables (even though they
	%might remain constant)
    
    %Plot Control Timecourse
    %{
    Aact=Simulation_Output.Data(MinTimeIndex:end,2);
    Bact=Simulation_Output.Data(MinTimeIndex:end,4);
    Cact=Simulation_Output.Data(MinTimeIndex:end,6);
    NegReg=Simulation_Output.Data(MinTimeIndex:end,7);
    Stimulus=Simulation_Output.Data(MinTimeIndex:end,8);
    Time=Simulation_Output.Time(MinTimeIndex:end,1);
    plot(Time,Aact,'Linewidth', 4)
    xlabel('Time (AU)')
    ylabel('Species abundance (AU)')
    set(gca,'FontWeight','bold','fontsize',30)
    hold all
    plot(Time,Bact,'Linewidth', 4)
    plot(Time,Cact,'Linewidth', 4)
    plot(Time, NegReg,'Linewidth', 4)
    plot(Time, Stimulus,'Linewidth', 4)
    legend('A_a_c_t','B_a_c_t','C_a_c_t','NegReg','Stimulus')
    %}

	%Define relevant data range in species and time
	EventTime=round(Timecourse_Duration/2);
	[number, MinTimeIndex ] = min( abs( Simulation_Output.Time-EventTime ) );%Identify index closes to Event time in order not to include peaks/troughs arising from system equilibration
 	
    %Check for peaks
 	RelevantData=Simulation_Output.Data(MinTimeIndex:end,6); %INDEX FOR Cact DATA
	[PeakMagnitude,location,width,prominence]=findpeaks(RelevantData,'MinPeakProminence',1); %Minimum peak magnitude of 1 a.u
	Peak=max(PeakMagnitude);
	if numel(Peak) == 0   %if no peaks are found
		Peak = 0;
	end
	Control_val=[Control_val;Peak];  %Extract control data  
    
	%Perform simulation under constitutive signal for each constitutive signal feeding

	for ConstSig=1:3
        cs2=[];
        CS_Simulation_Output=[];
        CS_Output=[];
		if ConstSig == 1    %For upstream CS, sample a random paramter value and set other two CS reactions to rate 0
	    	Kcs1 = ParamVals(randi(numel(ParamVals)));
			Kcs2 = 0;
			Kcs3 = 0;
			modelObj.Parameters(8).value=(Kcs1);
			modelObj.Parameters(9).value=(Kcs2);
			modelObj.Parameters(10).value=(Kcs3);
		elseif ConstSig== 2 %For middlestream CS, sample a random paramter value and set other two CS reactions to rate 0
	    	Kcs1 = 0;
	    	Kcs2 = ParamVals(randi(numel(ParamVals)));
			Kcs3 = 0;
			modelObj.Parameters(8).value=(Kcs1);
			modelObj.Parameters(9).value=(Kcs2);
			modelObj.Parameters(10).value=(Kcs3);
		elseif ConstSig== 3 %For downstream CS, sample a random paramter value and set other two CS reactions to rate 0
			Kcs1 = 0;
			Kcs2 = 0;
	    	Kcs3 = ParamVals(randi(numel(ParamVals)));
			modelObj.Parameters(8).value=(Kcs1);
			modelObj.Parameters(9).value=(Kcs2);
			modelObj.Parameters(10).value=(Kcs3);
        end

        %Simulate effect of constitutive signal
		cs2 = getconfigset(modelObj,'active');
		cs2.SolverType = 'ode45';  %ODE solver - ODE45 non-stiff, ode45 & ODE23 = stiff
		cs2.SolverOptions.AbsoluteTolerance= 1.0e-12;
		cs2.SolverOptions.RelativeTolerance= 1.0e-6;
		cs2.StopTime = Timecourse_Duration; %Simulation stop time
		%cs2.SolverOptions.LogDecimation = 200;   %how frequently you want to record the output of a stochastic simulation (ex. every 200 ime units)
		cs2.CompileOptions.UnitConversion = false;  %No unit conversion
		CS_Simulation_Output=sbiosimulate(modelObj); %Simulate model
		%sbioplot(Simulation_Output);  %plot simulation output
             
        %Plot CS Timecourse
        %{
        Aact=CS_Simulation_Output.Data(MinTimeIndex:end,2);
        Bact=CS_Simulation_Output.Data(MinTimeIndex:end,4);
        Cact=CS_Simulation_Output.Data(MinTimeIndex:end,6);
        NegReg=CS_Simulation_Output.Data(MinTimeIndex:end,7);
        Stimulus=CS_Simulation_Output.Data(MinTimeIndex:end,8);
        Time=CS_Simulation_Output.Time(MinTimeIndex:end,1);
        plot(Time,Aact,'Linewidth', 4)
        xlabel('Time (AU)')
        ylabel('Species abundance (AU)')
        set(gca,'FontWeight','bold','fontsize',30)
        hold all
        plot(Time,Bact,'Linewidth', 4)
        plot(Time,Cact,'Linewidth', 4)
        plot(Time, NegReg,'Linewidth', 4)
        plot(Time, Stimulus,'Linewidth', 4)
        legend('A_a_c_t','B_a_c_t','C_a_c_t','NegReg','Stimulus')
        %}


		%Extract constitutive signal data 
		CS_Output=CS_Simulation_Output.Data(:,1:end-NumOfParams); %Columns correspond to species in the order they appear in the modelObj. The last column indices corresponding to the
		%number of parameters are not included since when parameters are defined as non-constant (so they can be changed by an Event) then they are plotted as variables (even though they
		%might remain constant)

		%Define relevant data range in species and time
		EventTime=round(Timecourse_Duration/2);
		[number, MinTimeIndex ] = min( abs( CS_Simulation_Output.Time-EventTime ) );%Identify index closes to Event time in order not to include peaks/troughs arising from system equilibration
        
 		%Check for peaks
 		RelevantData=CS_Simulation_Output.Data(MinTimeIndex:end,6); %INDEX FOR Cact DATA
		[PeakMagnitude,location,width,prominence]=findpeaks(RelevantData,'MinPeakProminence',1); %Minimum peak magnitude of 1 a.u
		Peak=max(PeakMagnitude);
		if numel(Peak) == 0
			Peak = 0;
		end
		CS_val=[CS_val;Peak]; %Extract data for the effect of each constitutive signal
        
        if ConstSig == 3
            modelObj.Parameters(10).value=0; %Reset Kcs3 to zero
        end
	end   
	
	Diff=[(CS_val(1)/Control_val(1))*100; (CS_val(2)/Control_val(1))*100;  (CS_val(3)/Control_val(1))*100];	%Express changes in Cact as a percentage of the value for the control simulation
    Results=[Results Diff];	 %Extract results
	
end      
				 		      

xlswrite('Model_1_Results.xlsx',Results); %Save results


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot heatmap

%Create xticks and yticks vectors for figures
%xtks=1:numel(ParamNumParam);
ytks=1:3;
    
%Create Figures   
figure(1)
imagesc(Results)
caxis([75 125])   %colourjet mapping to data to highlight changes over 25%
%in responsiveness
%xticks(xtks)
%xtickangle(45)
CS_Names={'Upstream CS' ; 'Middlestream CS' ; 'Downstream CS'};
yticks(ytks)
yticklabels(CS_Names)
set(gca,'FontWeight','bold','fontsize',30)
colorbar
colormap(jet)
    
%Save Figure
%h=figure(TP);
%savefig(['Sensitivity_Timepoint_' num2str(Time)])
  

