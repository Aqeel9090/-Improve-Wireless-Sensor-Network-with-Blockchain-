function Sensors=addToLedger(Sensors,block)
for i=1:length(Sensors)
        Sensors(i).ledger.data = block.data;
        Sensors(i).ledger.currentHash = block.hash;
        Sensors(i).ledger.previousHash = block.previous_hash;
end
    
end