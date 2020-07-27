%% Environment

clc;
clear;

Ts = 1e-6;
Tf = 1e-6 * 5000 - 1e-6;

% create the observation info
observationInfo = rlNumericSpec([20 1], 'LowerLimit', -inf * ones(20, 1), 'UpperLimit', inf * ones(20, 1));
observationInfo.Name = 'observations';

% action Info¡
actionInfo = rlNumericSpec([20 1], 'LowerLimit', -inf * ones(20, 1), 'UpperLimit', inf * ones(20, 1));
actionInfo.Name = 'u_Input';

% define environment
env = rlFunctionEnv(observationInfo,actionInfo,'ddpg_Step_Func','ddpg_Reset_Func');
%% Critic
% number of neurons
L = 20;

% Define two structures, statePath and actionPath
% imageInputLayer([height, width, channels])
% fullyConnectedLayer(OutputSize) & additionLayer(NumInputs)
statePath = [
    imageInputLayer([20 1 1],'Name','observation','Normalization','none')
    depthConcatenationLayer(2,'Name','concat')
    fullyConnectedLayer(L,'Name','fc3')
    reluLayer('Name','relu3')
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
criticNetwork = connectLayers(criticNetwork,'action','concat/in2');

% Visualise the critic network
% plot(criticNetwork)

% Set options, including learn rate, gradient threshold and regularization
criticOptions = rlRepresentationOptions('LearnRate',1e-1, ...
    'GradientThreshold',1,'L2RegularizationFactor',1e-4);
% Create critic network with 'observation' and 'action'
critic = rlRepresentation(criticNetwork,observationInfo,actionInfo, ...
    'Observation',{'observation'},'Action',{'action'},criticOptions);
%% Actor
% This is exactly the same as critic network, but without additional layer
actorNetwork = [
    imageInputLayer([20 1 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(L,'Name','f')
    reluLayer('Name','r')
    fullyConnectedLayer(L,'Name','fc1')
    ];
% Visualise the actor network
% plot(layerGraph(actorNetwork))
%     scalingLayer('Name','ActorScaling1','Scale',2.5,'Bias',-0.5)
actorOptions = rlRepresentationOptions('LearnRate',1e-2,...
    'GradientThreshold',1,'L2RegularizationFactor',1e-4);
actor = rlRepresentation(actorNetwork,observationInfo,actionInfo,...
    'Observation',{'observation'},'Action',{'fc1'},actorOptions);

%% DDPG
% Similar arguments
agentOptions = rlDDPGAgentOptions(...
    'SampleTime',Ts,...
    'TargetSmoothFactor',1e-3,...
    'ExperienceBufferLength',1e6,...
    'DiscountFactor',0.99,...
    'MiniBatchSize',64);

agentOptions.NoiseOptions.InitialAction = 1e10;
agentOptions.NoiseOptions.Variance = 25;
% 6.25e10;
agentOptions.NoiseOptions.VarianceDecayRate = 1;
agent = rlDDPGAgent(actor,critic,agentOptions);

%%
maxepisodes = 5000000;
maxsteps = ceil(Tf/Ts);
trainOpts = rlTrainingOptions(...
    'MaxEpisodes',maxepisodes,...
    'MaxStepsPerEpisode',maxsteps,...
    'ScoreAveragingWindowLength',5,...
    'Verbose',false,...
    'Plots','training-progress',...
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',100,...
    'SaveAgentCriteria','EpisodeReward',...
    'SaveAgentValue',100);

%%
% Train the agent.
trainingStats = train(agent,env,trainOpts);
%%
% simOptions = rlSimulationOptions('MaxSteps',500);
% experience = sim(env,agent,simOptions);





