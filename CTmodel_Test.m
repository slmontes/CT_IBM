%% CT model 
%% TEST

clear all; 
clc; 
close all;

  N = 10000;                  %total population. This is N=Nf+Nm
  EqPeriod_Days = 20*365;     % days to simulate for equilibrium period
  Test_Days = 10*365;         % days to apply test/recall

  VERBOSE = true;             %Show process in command window (set to false
                              % when setting up a loop for repetitions)
  LOW_MEMORY = false;         %New model has only been tested with 
                              %LOW_MEMORY = false        

        % load preset of simulation parameters from external file
          load('params_het_CT.mat','params');
            
          params.RESTRICT_RATE = 30;

        % Initial parameters that can be added or modifed to the preset:
	      params.ENABLE_nonAMR_RECALL=false;
	      params.ENABLE_nonAMR_TRACE=false;

	      params.HR_SCREEN_EFFICACY = 0;
          params.SCREEN_SEEKED_HIGHRISK = false;
          params.SCREEN_NumContacts_HIGHRISK = false;
          params.TRACE_PARNTERS_NumContacts_HIGHRISK = false;
          params.HR_PSI = 0;
          params.Re_Screening_Period = 92; %92 for ~3 months
                                           %182 for ~6 months
                                           %365 for ~12 months
          params.ALPHA = 2.8;
          params.ALPHAM = 2.5;

	      params.P_SYMPTOMS = 0.1;

	      params.GAMMA= 0;

          params.PSI = 0.0; 
          params.BETA = [0.0022 0];   %Only using a single strain of CT, 
                                      % thus, second strain is set to 0
          params.p0=[0.04,0,0];

  %% Initialise model (create new model object)
  
        hetero_CT_model=CT_IBM_2021(N,params,[],VERBOSE,LOW_MEMORY);

  %% Run simulation for equilibrium_Days # of days 
  
        hetero_CT_model.simulate(EqPeriod_Days);

	%% run simulation for Test_Days # of days

        %Change parameters to introduce intervention 
        % Choose only one set for screening and comment the other:
        
        %Follow-up of patients seeking attention:
% 	    hetero_CT_model.SCREEN_SEEKED_HIGHRISK = true;
%         hetero_CT_model.HR_SCREEN_EFFICACY = 0.6; % 60% of probabilities to be screened after seeking treatment
						                          % Only works like this for HR_Seeking_Screening!!!
        
        hetero_CT_model.SCREEN_NumContacts_HIGHRISK = true;
        hetero_CT_model.HR_SCREEN_EFFICACY = 0.0012;  % This is a rate that works as GAMMA 
                                                       % Only works like this for
                                                       % SCREEN_NumContacts_HIGHRISK!!!
                                                       
        % Choose if the interventio will include partner tracing or not:
        hetero_CT_model.TRACE_PARNTERS_NumContacts_HIGHRISK = true;
        hetero_CT_model.HR_PSI = 0.4;
	    hetero_CT_model.PRESCREEN_TRACED = true;

        hetero_CT_model.simulate(Test_Days);
        
        %% Pre-defined plots
        
        hetero_CT_model.plots
        
                