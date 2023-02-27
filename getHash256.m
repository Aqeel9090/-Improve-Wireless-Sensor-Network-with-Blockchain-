function [callHashTime,HashExecTime, hashStr] = getHash256(txtMessage,callHashStartTime)
callHashTime = toc(callHashStartTime);
HashExecBeginTime = tic;
% Available options are 'SHA1', 'SHA256', 'SHA384', 'SHA512', 'MD5'
%algorithm = 'SHA256';   

% SHA1 category
%hasher = System.Security.Cryptography.HashAlgorithm.Create(algorithm);  % DEFAULT

% SHA2 category
hasher = System.Security.Cryptography.HashAlgorithm.Create('SHA256');  
%hasher = System.Security.Cryptography.HashAlgorithm.Create('SHA384');  
%hasher = System.Security.Cryptography.HashAlgorithm.Create('SHA512');

% SHA3 category:   Does not appear to be supported

% MD5 category
%hasher = System.Security.Cryptography.HashAlgorithm.Create(algorithm);

% GENERATING THE HASH:
%str = 'Now is the time for all good men to come to the aid of their country';
str = txtMessage;
hash_byte = hasher.ComputeHash( uint8(str) );  % System.Byte class
hash_uint8 = uint8( hash_byte );               % Array of uint8
hash_hex = dec2hex(hash_uint8);                % Array of 2-char hex codes


% Generate the hex codes as 1 long series of characters
hashStr = str([]);
nBytes = length(hash_hex);
for k=1:nBytes
    hashStr(end+1:end+2) = hash_hex(k,:);
end
%HashExecEndTime = toc;
HashExecTime = toc(HashExecBeginTime);

%fprintf(1, '\n\tThe %s hash is: "%s" [%d bytes]\n\n', algorithm, hashStr, nBytes);


% SIZE OF THE DIFFERENT HASHES:
%       SHA1:  20 bytes = 20 hex codes =  40 char hash string
%     SHA256:  32 bytes = 32 hex codes =  64 char hash string
%     SHA384:  48 bytes = 48 hex codes =  96 char hash string
%     SHA512:  64 bytes = 64 hex codes = 128 char hash string
%        MD5:  16 bytes = 16 hex codes =  32 char hash string