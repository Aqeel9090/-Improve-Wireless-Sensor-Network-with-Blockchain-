function [result, percentageAgree]=ConsensusHardware(Sensors,block, CH, consensusRate)
%% The function checks if the block is OK to be added to the ledger.
% Consensus is used among all sensors. 
% 
% 
prevHash = block.previous_hash;
% curHash = block.hash;
% data = block.data;
votes = zeros(1,length(CH));

for i=1:length(Sensors) % the last one is the base station
%     if((Sensors(i).ledger.previousHash == 0) & ...
%             (Sensors(i).ledger.data == 0) & ...
%             (Sensors(i).ledger.currentHash == 0))
%         Sensors(i).ledger.data = data;
%         Sensors(i).ledger.currentHash = curHash;
%         Sensors(i).ledger.previousHash = prevHash;
%         %disp('addToledger function first Round adding ledger');
%         
%     else % Check prevoius information in ledger
%         %disp('addToledger function voting');

        oldData = Sensors(i).ledger.data;
        oldPrevHash = Sensors(i).ledger.previousHash;
        oldHash = Sensors(i).ledger.currentHash;
        % For software
        %[~, ~,predOldHash] = DataHash(num2str(oldData),tic, 'SHA-256', 'ascii');
        % Using hardware 
        [~,~,predOldHash] = SecurityKeysHash(oldData);
        if(strcmp(predOldHash, oldHash) && strcmp(oldHash,prevHash))
        %if(strcmp(oldHash, predOldHash))
        %if(1)
            votes(i) = 1;
        else
            votes(i)=0;
        end
%     end
    % Deduct Hash energy from all sensors
    originalMessage = block.data;
    SoftwareHashEnergy = 2.32*1e-9 * length(de2bi(round(abs(originalMessage))));
    Sensors(i).E = Sensors(i).E-SoftwareHashEnergy;
end
%% Check votes
totalAgree = sum(votes);
percentageAgree = round(100*totalAgree/length(Sensors));
if(percentageAgree >= consensusRate)
    result = 1;
    % Add To ledger
    addToLedger(Sensors,block);
else
    result = 0;
end
       
   
        