function [ResultsPerRound,sensorsPerCHperRound,CHperRound, maxBlkSizePerRound, minBlkSizePerRound, TotalExtraFileSizeinBytes,LifeExpectancy,Sensors,CompRatio]=simulateWSN(method,Compression,noOfNodes,XYarea,rounds)
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
%     rmax = 500;                             % Round max
    [Area,Model]=setParameters(noOfNodes,XYarea,rounds);     		%Set Parameters Sensors and Network

    %%%%%%%%%%%%%%%%%%%%%%%%% configuration Sensors %%%%%%%%%%%%%%%%%%%%
    CreateRandomSen(Model,Area);            %Create a random scenario
    load Locations                          %Load sensor Location
    Sensors=ConfigureSensors(Model,noOfNodes,X,Y);


    %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%%
    countCHs=0;         %counter for CHs
    flag_first_dead=0;  %flag_first_dead
    deadNum=0;          %Number of dead nodes

    initEnergy=0;       %Initial Energy
    for i=1:noOfNodes
          initEnergy=Sensors(i).E+initEnergy;
    end

    SRP=zeros(1,Model.rmax);    %number of sent routing packets
    RRP=zeros(1,Model.rmax);    %number of receive routing packets
    SDP=zeros(1,Model.rmax);    %number of sent data packets 
    RDP=zeros(1,Model.rmax);    %number of receive data packets 

    Sum_DEAD=zeros(1,Model.rmax);
    CLUSTERHS=zeros(1,Model.rmax);
    AllSensorEnergy=zeros(1,Model.rmax);

    %% %%%%%%%%%%%%%%%%%%%%%%% Start Simulation %%%%%%%%%%%%%%%%%%%%%%%%%
    global srp rrp sdp rdp
    srp=0;          %counter number of sent routing packets
    rrp=0;          %counter number of receive routing packets
    sdp=0;          %counter number of sent data packets 
    rdp=0;          %counter number of receive data packets 
    % Base Station (BS) is the node n+1

        BaseStation = noOfNodes+1;
    
    %Sink broadcast start message to all nodes
    Sender=BaseStation;     %Sink
    Receiver=1:noOfNodes;   %All nodes
    Sensors=SendReceivePackets(Sensors,Model,Sender,'Hello',Receiver,1,1);

    % All sensor send location information to Sink .
     Sensors=disToSink(Sensors,Model);
    % Sender=1:n;     %All nodes
    % Receiver=BaseStation;   %Sink
    % Sensors=SendReceivePackets(Sensors,Model,Sender,'Hello',Receiver);

    %Save metrics
    SRP(1)=srp;
    RRP(1)=rrp;  
    SDP(1)=sdp;
    RDP(1)=rdp;

    %% Main loop program
    L = height(keys);
    rounds = Model.rmax; 
    previousHash = 0;
    LifeExpectancy = 0;
    HashCallTimeForRound = zeros(rounds,1);
    HashExecTimeForRound = zeros(rounds,1);
    HashComplexityForRound= zeros(rounds,1);
    EncryptComplexityForRound=zeros(rounds,1);
    CHElectTime  = zeros(rounds,1);
    sendToCHtime = zeros(rounds,1);
    sendToBSTime = zeros(rounds,1);
    directToBS = zeros(rounds,1);
    TotalComplexity= zeros(rounds,1);
    blockSize = [];
    maxBlkSizePerRound =[];
    minBlkSizePerRound = [];
    CHperRoundAll = [];
    s = struct('CHperRound',0,'sensorsPerCH',[]);
    sensorsPerCHperRoundAll = repmat(s,[rounds 1]);
    for r=1:1:rounds
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
        %This section Operate for each epoch   
        member=[];              %Member of each cluster in per period
        countCHs=0;             %Number of CH in per period
        %counter for bit transmitted to Bases Station and Cluster Heads
        srp=0;          %counter number of sent routing packets
        rrp=0;          %counter number of receive routing packets
        sdp=0;          %counter number of sent data packets to sink
        rdp=0;          %counter number of receive data packets by sink
        %initialization per round
        SRP(r+1)=srp;
        RRP(r+1)=rrp;  
        SDP(r+1)=sdp;
        RDP(r+1)=rdp;   
    %    pause(0.001)    %pause simulation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        Sensors=resetSensors(Sensors,Model);
        %allow to sensor to become cluster-head. LEACH Algorithm  
        AroundClear=10;
        if(mod(r,AroundClear)==0) 
            for i=1:1:noOfNodes
                Sensors(i).G=0;
            end
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% plot sensors %%%%%%%%%%%%%%%%%%%%%%%
    %    deadNum=ploter(Sensors,Model,TopolgyFigure);
        %Save r'th period When the first node dies
        if (deadNum>=1)      
            if(flag_first_dead==0)
                first_dead=r;
                flag_first_dead=1;
            end  
        end

    %% %%%%%%%%%%%%%%%%%%%%% cluster head election %%%%%%%%%%%%%%%%%%%
        %Selection Candidate Cluster Head Based on LEACH Set-up Phase
        [TotalCH,Sensors]=SelectCH(Sensors,Model,r); 
            
        
        CHElectTimeStart=tic;
        %Broadcasting CHs to All Sensor that are in Radio Rage CH.
        for i=1:length(TotalCH)
            Sender=TotalCH(i).id;
            SenderRR=Model.RR;
            Receiver=findReceiver(Sensors,Model,Sender,SenderRR);   
            Sensors=SendReceivePackets(Sensors,Model,Sender,'Hello',Receiver,r,1);
        end 
        
        CHElectTime(r,1) = toc(CHElectTimeStart);
        %Sensors join to nearest CH 
        Sensors=JoinToNearestCH(Sensors,Model,TotalCH);
        % Calculate size of cluster and CHs
        CHperRoundAll(r,1) = length(TotalCH);
        sensorsPerCHperRoundAll(r).CHperRound = CHperRoundAll(r,1);
        for i=1:length(TotalCH)
            memberNode =0;
            kk=1;
            for j = 1:length(Sensors)
                if(isequal(Sensors(j).MCH, TotalCH(i).id))
                    sensorsPerCHperRoundAll(r).sensorsPerCH(kk) = Sensors(j).id;
                    memberNode = memberNode+1;
                    kk=kk+1;
                end
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%% end of cluster head election phase %%%%%%

    %%%%%%%%%%%%%%%%%%%%%%% plot network status in end of set-up phase 

       % ploter(Sensors,Model,figure1);                  %Plot sensors
        for i=1:noOfNodes

            if (Sensors(i).type=='N' && Sensors(i).dis2ch<Sensors(i).dis2sink && ...
                    Sensors(i).E>0)

                XL=[Sensors(i).xd ,Sensors(Sensors(i).MCH).xd];
                YL=[Sensors(i).yd ,Sensors(Sensors(i).MCH).yd];
                %hold on
    %            figure(TopolgyFigure),line(XL,YL)

            end

        end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% steady-state phase %%%%%%%%%%%%%%%%%
        
        if(Compression)
           NumPacket=Model.NumPacket/2;
        else
            NumPacket=Model.NumPacket;
        end
            

        sendToCHstart = tic;
%         for i=1:1:1%NumPacket 
% 
%             %Plotter     
%     %        deadNum=ploter(Sensors,Model,TopolgyFigure);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All sensor send data packet to  CH 
% %             for j=1:length(TotalCH)
% % 
% %                 Receiver=TotalCH(j).id;
% %                 Sender=findSender(Sensors,Model,Receiver); 
% % 
% %     %             nonce = rand(1); % it is random used one time for each comm
% %     %             timestamp = datestr(datetime);
% %     %             blockInitials = struct('index', 1+r, 'timestamp',timestamp, 'data',compressedData.CompressedData(r), 'nonce',nonce,'hash',1234453,'previous_hash',0);
% %     %             block = myBlock(blockInitials);    
% %                 Sensors=SendReceivePackets(Sensors,Model,Sender,'Data',Receiver,r,1);
% % 
% %             end
% 
%         end

    sendToCHtime(r,1)=toc(sendToCHstart);


        %% CH to Base Station 
        % Send Data packet from CH to the base station after Data aggregation
        % 1- Spread the public key from Base Station to all CHs in this round
      % Generate keys
              [Modulus, PublicExponent, PrivateExponent] = GenerateKeyPair;
        sendToBSStart=tic;
        
    for packet = 1:NumPacket
        HashCallTime= zeros(length(TotalCH),1);
        HashExecTime = zeros(length(TotalCH),1);
        HashComplexity = zeros(length(TotalCH),1);
        EncryptComplexity = zeros(length(TotalCH),1);
        TotalComplexityTime = zeros(length(TotalCH),1);
        for i=1:length(TotalCH)
                % Get and Compress data
                originalMessage = compressedData.CompressedData;
                [message,CompRatio] = compressChaotic(originalMessage);
                % Get  Hash
                callHashStartTime = tic;
                %%% Works only for Windows
                %[HashCallTime(i,1), HashExecTime(i,1),currentHash] = getHash256(message,callHashStartTime);
                %%% Otherwise, and for Non-Windows OS use this function:
                if(method == 1)
                    [HashCallTime(i,1), HashExecTime(i,1),currentHash] = DataHash(num2str(originalMessage),callHashStartTime, 'SHA-256', 'ascii');
                else
                    [HashCallTime(i,1), HashExecTime(i,1),currentHash] = DataHash(num2str(message),callHashStartTime, 'SHA-256', 'ascii');
                end
                HashComplexity(i,1) = HashExecTime(i,1);

                callEncryptStart=tic;
                 [callEncryptTime(i,1),EncryptExecTime(i,1),DataTosend] = Encrypt(Modulus, PublicExponent, num2str(message),callEncryptStart);
                 EncryptComplexity(i,1) = EncryptExecTime(i,1);

                 Receiver=BaseStation;               %Sink
                 Sender=TotalCH(i).id;       %CH 
             
            if(method == 1)  % Software without compression
                CompRatioSoftware = 1;
                nonce = rand(1); % it is random used one time for each comm
                timestamp = datestr(datetime);
                blockInitials = struct('index', 1+r, 'timestamp',timestamp,...
                    'data',originalMessage, 'nonce',nonce,...
                    'hash',currentHash,'previous_hash',previousHash);
                block = myBlock(blockInitials);    
                previousHash = currentHash;
                
%                 disp(['ROUND  ', num2str(r)]);
%                 disp(['No of Nodes: ', num2str(noOfNodes)]);
%                 disp(['Size of Sensors: ', num2str(length(Sensors))]);
%                 
                Sensors = SendReceivePackets(Sensors,Model,Sender,block,Receiver,r,CompRatioSoftware);
    %        Sensors=SendReceivePackets(Sensors,Model,Sender,'Data',Receiver);

                % Update Energy for Cluster heads : 
                % According to 
                % Energy consumption and execution time consumption of embedded
                % systems By Callo et. al. (page 12 table3), 
                % Each Joule of energy consumption =Execution time*0.0556 
                TotalComplexityTime(i,1) = HashComplexity(i,1) + EncryptComplexity(i,1);
                % Sensors(Sender).E = Sensors(Sender).E - (TotalComplexityTime(i,1)*0.0556); 
                SoftwareHashEnergy = 68.710526*1e-9 * length(de2bi(round(abs(originalMessage))));
                Sensors(Sender).E = Sensors(Sender).E - SoftwareHashEnergy; 
%                 blockLength = length(block.data) + length(block.hash) +...
%                     length(block.index)+...
%                     length(block.timestamp) + ...
%                     length(block.nonce)+ length(block.previous_hash);
                blockLength = getBlockSize(block);
                blockSize = [blockSize,blockLength];
                
                %% Consensus
                if(r==1)
                    Sensors = addToLedger(Sensors,block);
                else
                     % consensus
                     [result, percentageAgree] =Consensus(Sensors,block,TotalCH, 51);
                    if(result)
%                         disp(['Round ', num2str(r) ,'. Block Accepted, precentage Agree=', num2str(percentageAgree)]);
                    else
%                         disp(['Round ', num2str(r) ,'. Block Denied, it is an attack, precentage Agree=', num2str(percentageAgree)]);
                    end
                end

            elseif(method == 0)  % Hardware
                nonce = rand(1); % it is random used one time for each comm
                timestamp = datestr(datetime);
                [~,~,currentHash] = SecurityKeysHash(message);
                blockInitials = struct('index', 1+r, 'timestamp',timestamp,...
                    'data',message, 'nonce',nonce,...
                    'hash',currentHash,'previous_hash',previousHash);
                block = myBlock(blockInitials);
                
                previousHash = currentHash;
                Sensors = SendReceivePackets(Sensors,Model,Sender,block,Receiver,r,CompRatio);
                %% Consensus
                if(r==1)
                    Sensors =addToLedger(Sensors,block);
                else
                     % Add to ledger using consensus
                     [result, percentageAgree] =ConsensusHardware(Sensors,block,TotalCH, 51);
                    if(result)
%                         disp(['Round ', num2str(r) ,'. Block Accepted, precentage Agree=', num2str(percentageAgree)]);
                    else
%                         disp(['Round ', num2str(r) ,'. Block Denied, it is an attack, precentage Agree=', num2str(percentageAgree)]);
                    end
                end
%%
    %        Sensors=SendReceivePackets(Sensors,Model,Sender,'Data',Receiver);

                % Sensors(Sender).E = Sensors(Sender).E - Model.Emp; 
                % Or using the FPGA implementation test on Dec 2020 which
                % has the following results:
                % For a message of N = 3800 bits
                % compression time T = 7.405 ns
                % and FPGA power consumption was P = 10.078 Watt.
                % Hence, energy E = P x T   
                % In average, Eav = P/N  x T/N   average energy per bit.
                % Eav = 0.00265 x  0.001949 x 10^(-9)
                % Eav = 2.65 x 10^(-3) x  1.949 x 10^(-12)
                % Eav  = 1.9639e-11  Joule/Bit
%                 blockLength = length(block.data) + length(block.hash) +...
%                     length(block.index)+...
%                     length(block.timestamp) + ...
%                     length(block.nonce)+ length(block.previous_hash);
%                 blockSize = [blockSize,blockLength];
                
                
                %HW_Energy = 1.9639e-11 * blockLength;
%                 HW_Energy = 1.9639e-11 *
%                 length(de2bi(round(abs(message))));% OLD
                HW_Energy = 1.5530e-12 * length(de2bi(round(abs(message))));
%                 RSAHardwareEnergy = 0.082*2e-9;
%                 HW_Energy = HW_Energy + RSAHardwareEnergy;
                Sensors(Sender).E = Sensors(Sender).E - HW_Energy- 2.32*1e-9 * length(de2bi(round(abs(originalMessage)))); 
                blockLength = getBlockSize(block);
                blockSize = [blockSize,blockLength];
                              
            elseif(method == 2) % software with compression
                nonce = rand(1); % it is random used one time for each comm
                timestamp = datestr(datetime);
                blockInitials = struct('index', 1+r, 'timestamp',timestamp,...
                    'data',message, 'nonce',nonce,...
                    'hash',currentHash,'previous_hash',previousHash);
                block = myBlock(blockInitials);    
                previousHash = currentHash;             
                Sensors = SendReceivePackets(Sensors,Model,Sender,block,Receiver,r,CompRatio);
    %        Sensors=SendReceivePackets(Sensors,Model,Sender,'Data',Receiver);
                % Update Energy for Cluster heads : 
                % According to 
                % Energy consumption and execution time consumption of embedded
                % systems By Callo et. al. (page 12 table3), 
                % Each Joule of energy consumption =Execution time*0.0556 
                TotalComplexityTime(i,1) = HashComplexity(i,1) + EncryptComplexity(i,1);
                % Sensors(Sender).E = Sensors(Sender).E - (TotalComplexityTime(i,1)*0.0556); 
                compressionEnergy = 3.58e-10* length(de2bi(round(abs(originalMessage))));
                SoftwareHashEnergy = 68.710526*1e-9 * length(de2bi(round(abs(message))));
                Sensors(Sender).E = Sensors(Sender).E - (SoftwareHashEnergy+compressionEnergy); 
%                 blockLength = length(block.data) + length(block.hash) +...
%                     length(block.index)+...
%                     length(block.timestamp) + ...
%                     length(block.nonce)+ length(block.previous_hash);
                blockLength = getBlockSize(block);
                blockSize = [blockSize,blockLength];
                
                %% Consensus
                if(r==1)
                    Sensors = addToLedger(Sensors,block);
                else
                     % consensus
                     [result, percentageAgree] =Consensus(Sensors,block,TotalCH, 51);
                    if(result)
%                         disp(['Round ', num2str(r) ,'. Block Accepted, precentage Agree=', num2str(percentageAgree)]);
                    else
%                         disp(['Round ', num2str(r) ,'. Block Denied, it is an attack, precentage Agree=', num2str(percentageAgree)]);
                    end
                end
                
            end

        end 
    end
        sendToBSTime(r,1) = toc(sendToBSStart);
       %% Find Complexity for each round         
       HashCallTimeForRound(r,1) = mean(HashCallTime);
        % Rempve NaN
        ind = find(isnan(HashCallTimeForRound));
        if (ind==1)
            HashCallTimeForRound(1,1)=max(HashCallTimeForRound);
        else
            HashCallTimeForRound(ind) = arrayfun(@(x) nanmean(HashCallTimeForRound(x-1)), ind);
        end

        HashExecTimeForRound(r,1) = mean(HashExecTime);
        % Rempve NaN
        ind = find(isnan(HashExecTimeForRound));
         if (ind==1)
            HashExecTimeForRound(1,1)=max(HashExecTimeForRound);
        else
            HashExecTimeForRound(ind) = arrayfun(@(x) nanmean(HashExecTimeForRound(x-1)), ind);
        end
        % Complexity is only sum because at each i it the program calls
        % Hash function one time only. Therefore, multiplication part is 
        % negligible. Only multiplied by 1 ptq=1 call.
        % Reference:A Novel Function Complexity-Based Code Migration Policy 
        % for Reducing Power Consumption. By: Hayeon Choi, Youngkyoung Koo, 
        % and Sangsoo Park
        HashComplexityForRound(r,1) = sum(HashComplexity);
         % Rempve NaN
        ind = find(isnan(HashComplexityForRound));
        if (ind==1)
            HashComplexityForRound(1,1)=max(HashComplexityForRound);
        else
            HashComplexityForRound(ind) = arrayfun(@(x) nanmean(HashComplexityForRound(x-1)), ind);
        end
        EncryptComplexityForRound(r,1) = sum(EncryptComplexity);
         % Rempve NaN
        ind = find(isnan(EncryptComplexityForRound));
         if (ind==1)
            EncryptComplexityForRound(1,1)=max(EncryptComplexityForRound);
        else
            EncryptComplexityForRound(ind) = arrayfun(@(x) nanmean(EncryptComplexityForRound(x-1)), ind);
        end
    %% % send data packet directly from other nodes(that aren't in each cluster) to Sink
        if(isempty(TotalCH))
            originalMessage = compressedData.CompressedData;
             [message,CompRatio] = compressChaotic(originalMessage);
            if(method == 0)
                [~,~,currentHash] = SecurityKeysHash(message);
            elseif(method == 1)
                [~, ~,currentHash] = DataHash(num2str(originalMessage),callHashStartTime, 'SHA-256', 'ascii');
                message = originalMessage;
                
            else
                [~, ~,currentHash] = DataHash(num2str(message),callHashStartTime, 'SHA-256', 'ascii');
            end
             
            nonce = rand(1); % it is random used one time for each comm
                timestamp = datestr(datetime);
                blockInitials = struct('index', 1+r, 'timestamp',timestamp,...
                    'data',message, 'nonce',nonce,...
                    'hash',currentHash,'previous_hash',previousHash);
                block = myBlock(blockInitials);    
                previousHash = currentHash;
        end
        directToBSstart= tic;
        for i=1:noOfNodes
            if(isequal(Sensors(i).MCH,Sensors(BaseStation).id))
                Receiver=BaseStation;               %Sink
                Sender=Sensors(i).id;       %Other Nodes 
                Sensors=SendReceivePackets(Sensors,Model,Sender,block,Receiver,r,1);
            end
        end

    directToBS(r,1)=toc(directToBSstart);
    %% STATISTICS

        Sum_DEAD(r+1)=deadNum;

        SRP(r+1)=srp;
        RRP(r+1)=rrp;  
        SDP(r+1)=sdp;
        RDP(r+1)=rdp;

        CLUSTERHS(r+1)=countCHs;
%% Subtract RSA key energy from all nodes
% Only in the first round -- REMOVED in JUNE 2022
%         RSAHardwareEnergy = 0;
%         RSASoftwareEnergy = 2.35066 * 1e-7;
%         if(r==1)
%             for i=1:noOfNodes
%                 if(Sensors(i).E >0)
%                     if(method)
%                         Sensors(i).E = Sensors(i).E- RSASoftwareEnergy;
%                     else
%                         Sensors(i).E = Sensors(i).E- RSAHardwareEnergy;
%                     end
%                 end
%             end
%         end
        alive=0;
        SensorEnergy=0;
        for i=1:noOfNodes
            if Sensors(i).E>0
                alive=alive+1;
                SensorEnergy=SensorEnergy+Sensors(i).E;
            else
                Sensors(i).alive = 0;
            end
            
        end
%         noOfNodes = alive;
%         
%        Model.n = alive;
        %% Energy Calculations

        AliveSensors(r,1)=alive; %#ok

        SumEnergyAllSensor(r+1,1)=SensorEnergy; %#ok

        AvgEnergyAllSensor(r+1,1)=SensorEnergy/alive; %#ok

        AvgConsumedEnergy(r+1,1)=(initEnergy-SumEnergyAllSensor(r+1))/noOfNodes; %#ok
        TotalConEn(r+1,1)=(initEnergy-SumEnergyAllSensor(r+1)); %#ok
        
        % Setting up block size
        maxBlkSizePerRound(r+1,1) = max(blockSize);
        minBlkSizePerRound(r+1,1) = min(blockSize);
        
        
        En=0;
        for i=1:noOfNodes
            if Sensors(i).E>0
                En=En+(Sensors(i).E-AvgEnergyAllSensor(r+1,1))^2;
            end
        end

        EnergyVariance(r+1,1)=En/alive; %#ok
        EnergyVariance = EnergyVariance(2:end,1);

       %dead
       if(alive==0)  % Life expectancy is 100%(noOfNodes-alive)

           LifeExpectancy=r;  
           disp(['100% of nodes are dead at round ', num2str(r)]);
           break;
       end
       
    end % for r=0:1:rmax
    %% Calculate Size Differences

    % dirInfo = dir(dirName);  %# Where dirName is the directory name where the
    %                          %#   file is located
    % index = strcmp({dirInfo.name},fileName);  %# Where fileName is the name of
    %                                           %#   the file.
    % fileSize = dirInfo(index).bytes;  %# The size of the file, in bytes
    hashFile = dir('getHash256.m');
    hashFileSize = hashFile.bytes;

    encryptFile = dir('SoftwareKeysRSA/Encrypt.m');
    encryptFileSize = encryptFile.bytes;

    compressFile = dir('Compression/compressChaotic.m');
    compressFileSize = compressFile.bytes;

    decompressFile = dir('Compression/Decomp.m');
    decompressFileSize = decompressFile.bytes;

    encrypt2File = dir('SoftwareKeysRSA/ModularExponentiation.m');
    encrypt2FileSize = encrypt2File.bytes;

    % Total Extra complexity in ms
    TotalExtraTimeComplexity = 1000*(HashComplexityForRound + EncryptComplexityForRound);
      % Total Complexity of the system in ms
        TotalComplexity1 = 1000*(CHElectTime + sendToCHtime + sendToBSTime + directToBS);

     if(method==0) % if hardware is used, subtract the hardware functions
        TotalComplexity = TotalComplexity1 - TotalExtraTimeComplexity;
%         if(TotalComplexity2<0)
%             TotalComplexity=TotalComplexity1;
%         else
%             TotalComplexity = TotalComplexity2;
%         end
        TotalExtraFileSizeinBytes = hashFileSize + encryptFileSize ...
                                    + encrypt2FileSize+ compressFileSize...
                                    + decompressFileSize;
    else
        TotalExtraFileSizeinBytes = 0;
        TotalComplexity = TotalComplexity1;
    end
   
 
%% Arrange Results
%maxBlkSizePerRound = max(blockSize);
%minBlkSizePerRound = min(blockSize);

CHperRound = CHperRoundAll;
sensorsPerCHperRound = sensorsPerCHperRoundAll(2:end,1);


TotalEnergy=SumEnergyAllSensor(2:end,1);
TotalConsumedEnergy=TotalConEn(2:end,1);
AvgConsumedEnergy = AvgConsumedEnergy(2:end,1);
AvgEnergyAllSensor=AvgEnergyAllSensor(2:end,1);

TotalComplexity = TotalComplexity(1:length(TotalEnergy));
ResultsPerRound = table(TotalEnergy,AvgEnergyAllSensor,...
                    TotalConsumedEnergy,AvgConsumedEnergy,EnergyVariance,...
                    TotalComplexity,AliveSensors);
                
% ResultsPerRound = table(TotalEnergy,AvgEnergyAllSensor,...
%                     TotalConsumedEnergy,AvgConsumedEnergy,EnergyVariance,...
%                     TotalComplexity,AliveSensors,CHperRound,sensorsPerCHperRound,...
%                     maxBlkSizePerRound,minBlkSizePerRound);

%disp(['Total Extra File Size: ', num2str(TotalExtraFileSize),'  Bytes']);



