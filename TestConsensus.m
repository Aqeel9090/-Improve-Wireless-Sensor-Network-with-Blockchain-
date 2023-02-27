clc,clear all, close all;
%addpath(genpath('+leach')) %% Adding path to working directory
%addpath('+leach');
addpath(genpath([pwd,'/','SoftwareKeysRSA']));
addpath(genpath([pwd,'/','leach']));
addpath(genpath([pwd,'/','blockChain']));
addpath(genpath([pwd,'/','ElectricalCircuit']));
addpath(genpath([pwd,'/','Compression']));
%% Run the Electrical Circuit
[keys, compressedData] = runCircuit();
%h = double(cell2mat(keys.HashHex(1)));  % Convert Hash to double
%% Define Block in the Chain
% nonce = rand(1); % it is random used one time for each comm
% timestamp = datestr(datetime);
% previousHash = '0';
% blockInitials = struct('index', 1, 'timestamp',timestamp, 'data',345566, 'nonce',nonce,'hash',1234453,'previous_hash',previousHash);
% block = myBlock(blockInitials);

%% Create sensor nodes, Set Parameters and Create Energy Model 
%%%%%%%%%%%%%%%%%%%%%%%%% Initial Parameters %%%%%%%%%%%%%%%%%%%%%%%
%     n=20;                                  %Number of Nodes in the field
%     rmax = 500; % Round max
noOfNodes=100;
rounds=100;
XYarea=100;
[Area,Model]=setParameters(noOfNodes,XYarea,rounds);     		%Set Parameters Sensors and Network

%%%%%%%%%%%%%%%%%%%%%%%%% configuration Sensors %%%%%%%%%%%%%%%%%%%%
CreateRandomSen(Model,Area);            %Create a random scenario
load Locations                          %Load sensor Location
Sensors=ConfigureSensors(Model,noOfNodes,X,Y);
%% Adding Block
message = 11234;
previousHash = 0;
callHashStartTime = tic;
[~, ~,currentHash] = DataHash(num2str(message),tic, 'SHA-256', 'ascii');
nonce = rand(1); % it is random used one time for each comm
timestamp = datestr(datetime);
blockInitials = struct('index', 1, 'timestamp',timestamp,...
    'data',message, 'nonce',nonce,...
    'hash',currentHash,'previous_hash',previousHash);
block = myBlock(blockInitials);
previousHash = currentHash;

%% Add to ledger
Sensors(1).ledger.data = block.data;
Sensors(1).ledger.currentHash = block.hash;
Sensors(1).ledger.previousHash = block.previous_hash;



