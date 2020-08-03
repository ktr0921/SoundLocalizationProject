%% Environment

clc;
clear;

Ts = 1e-6;
Tf = 1e-6 * 500000 - 1e-6;
configuration; % Configuration

% create the observation info; [20, 1] vector
observationInfo = rlNumericSpec([NR 1], 'LowerLimit', -inf * ones(20, 1), 'UpperLimit', inf * ones(20, 1));
observationInfo.Name = 'observations';
% action Info
actionInfo = rlNumericSpec([NR 1], 'LowerLimit', -inf * ones(20, 1), 'UpperLimit', inf * ones(20, 1));
actionInfo.Name = 'u_Input';
% define environment
env = rlFunctionEnv(observationInfo, actionInfo, 'ddpg_Step_Func', 'ddpg_Reset_Func');

%% Critic
% number of neurons
L = 20;

% Define two structures, statePath and actionPath
% imageInputLayer([height, width, channels])
% fullyConnectedLayer(OutputSize) & additionLayer(NumInputs)

%    depthConcatenationLayer(2,'Name','concat')
%    fullyConnectedLayer(L,'Name','fc3')
%    reluLayer('Name','relu3')
%criticNetwork = connectLayers(criticNetwork,'action','concat/in2');

statePath = [
    imageInputLayer([NR 1 1],'Name','observation','Normalization','none')
    additionLayer(2,'Name','add')
    fullyConnectedLayer(1,'Name','CriticOutput')
    ];
actionPath = [
    imageInputLayer([20 1 1],'Name','action','Normalization','none')
    ];

% Define critic network
% Map the only layer (statePath) with layerGraph(), add layer with
% addLayers(), and connect 'fc5' with  'add' as second input 'in2'
criticNetwork = layerGraph(statePath);
criticNetwork = addLayers(criticNetwork, actionPath);
criticNetwork = connectLayers(criticNetwork,'action','add/in2');

% Visualise the critic network
% plot(criticNetwork)

% Set options, including learn rate, gradient threshold and regularization
criticOptions = rlRepresentationOptions('LearnRate',1e-3, ...
    'GradientThreshold',1,'L2RegularizationFactor',1e-4);
% Create critic network with 'observation' and 'action'
critic = rlRepresentation(criticNetwork,observationInfo,actionInfo, ...
    'Observation',{'observation'},'Action',{'action'},criticOptions);
disp(1);
%% Actor
% This is exactly the same as critic network, but without additional layer
actorNetwork = [
    imageInputLayer([20 1 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(L,'Name','ActorOutput')
    ];
% Visualise the actor network
% plot(layerGraph(actorNetwork))
%     scalingLayer('Name','ActorScaling1','Scale',2.5,'Bias',-0.5)
actorOptions = rlRepresentationOptions('LearnRate',1e-4,...
    'GradientThreshold',1,'L2RegularizationFactor',1e-4);
actor = rlRepresentation(actorNetwork,observationInfo,actionInfo,...
    'Observation',{'observation'},'Action',{'ActorOutput'},actorOptions);
disp(2);
%% DDPG
% Similar arguments
agentOptions = rlDDPGAgentOptions(...
    'SampleTime',Ts,...
    'TargetSmoothFactor',1e-3,...
    'ExperienceBufferLength',1e6,...
    'DiscountFactor',0.99,...
    'MiniBatchSize',64);
agentOptions.NoiseOptions.Mean = 0;
agentOptions.NoiseOptions.Variance = 1;
agentOptions.NoiseOptions.VarianceDecayRate = 1e-5;
agent = rlDDPGAgent(actor,critic,agentOptions);

%%
maxepisodes = 5000;
maxsteps = ceil(Tf/Ts);
trainOpts = rlTrainingOptions(...
    'MaxEpisodes',maxepisodes,...
    'MaxStepsPerEpisode',maxsteps,...
    'ScoreAveragingWindowLength',5,...
    'Verbose',false,...
    'Plots','training-progress',...
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',5,...
    'SaveAgentCriteria','EpisodeReward',...
    'SaveAgentValue',5);

%%
% Train the agent.
trainingStats = train(agent,env,trainOpts);





