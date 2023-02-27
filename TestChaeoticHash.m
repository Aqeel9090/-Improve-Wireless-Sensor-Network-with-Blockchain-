clc,clear all, close all;
% addpath(genpath([pwd,'/','ElectricalCircuit']));
% addpath(genpath([pwd,'/','Compression']));
% %% Run the Electrical Circuit
% [keys, compressedData] = runCircuit();
% newData = ceil(abs(compressedData.CompressedData));


%%
sizeOfHash = 256;
newData = 12;
%% Define samples of Logistic map
sampled = linspace(0.01,((1.6/3.5)^(1/3)), 1024);
%Xn=rand(1,1);  % Input signal
%XnFromSensor=randperm(128,1);  % Input signal as integer 
Xn = sampled(newData); 
Xn = mean(Xn);
a=0.5;
alpha =2; 
beta =3 ;
k=2;


%%
% Each sample can produce four

for upSample=1:4 
    X(upSample) = a*((Xn).^alpha).*(1-(Xn./k).^beta);
    Xn = X(upSample)+ Xn; % Feedback
end
% Step 2: get the receprocal of the values as Logistic map o/p are very
% small

LMO = 1./X;
% Step 3: Apply the following equation 
% Key1 = 100*(LMO(1)) + LMO(4)
% Key2 = 10*(LMO(2)) + LMO(3)
Key1 = (LMO(1) + LMO(4));
Key2 = (LMO(2)) + LMO(3);

Keys = 10*[Key1+Key1;2*Key2];
%Keys = Keys./mean(Keys);
Hash = uint32(abs(Keys));
HashHex = dec2hex(Hash,8);
xx = HashHex(:)';
HashHex = xx;

HashHex2 = dec2hex(zeros(sizeOfHash,1));
c=10;
for i = 1:sizeOfHash
    d = round((20+c));
    HashHex2(i) = d;
    if(mod(i,2))
        c=d-14;
    else
        c=sqrt(d+16);
    end
end
HashHex3 = HashHex2+125;
HashHex2 = bitxor(HashHex2+1,HashHex3);
xx = dec2hex(HashHex2)';
dd = xx(2,:);
dd2 = xx(1,:)+128;
dd(end-length(HashHex):end-1) = HashHex;
HashHex = dd;
