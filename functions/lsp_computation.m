function out=lsp_computation(Am,Su,range_all,q,nfft,fs)

%%% INPUT
% --- Am is a (1 x p) vector of linear coefficients
% --- Su is a the variance of the AR residuals
% --- range_all is a (b x 2) matrix of frequency bands, with each band in a row
% --- q is the lag for computing correlation functions
% --- nfft is the number of point of positive frequencies
% --- fs is the sampling frequency

%%% OUTPUT
% functions in the frequency domain have b+1 columns, where the 1st column
% is associated to the total values, the 2nd and the 3rd respectively to
% the 1st and the 2nd frequency range indicated in the variable "range_all"

% parameters
z = 1i*2*pi/fs; % complex variable normalized by fs
b=size(range_all,1); % number of bands

index_band=zeros(1,b); % index of the freq range to select among all the poles found by the function
warn=zeros(1,b); % if 0, there is no pole in the selected band 

H=zeros(nfft,1); % transfer function associated to the whole process

%% time-domain LSP
R = lsp_Yule(Am,Su,q);

SigmaX=R(:,:,1);
SigmaX_X=Su;

LSP=0.5*log(det(SigmaX)/det(SigmaX_X)); % time domain

%% spectral decomposition
A=[1 -Am]; 
[~,~,freq,pot,~,poles,f]=lsp_spect_dec(A,1,Su,fs,nfft);

P=length(poles); % total number of the poles 
Hk=zeros(nfft,P); % abs(H) contributes
Hk_squared=zeros(nfft,P); % abs(H) ^ 2 contributes

%% TF associated to a given pole
for n=1:nfft
    
    NUM=exp(-z*f(n)); % numerator of Hk(z)
    for k=1:P 
        DEN=(exp(-z*f(n))-poles(k)); % denominator of Hk(z)
        Hk_tmp=NUM./DEN; 
        Hk(n,k)=abs(Hk_tmp); % abs(H) = H*conj(H)
    end

    for k=1:P
        if k~=P
            if imag(poles(k))==0 % real pole
                Hk_squared(n,k)=Hk(n,k).^2; 
            else % complex pole
                tmp=find(poles==conj(poles(k)));
                if tmp>k
                    Hk_squared(n,k)=prod(Hk(n,[k,tmp])).^2;
                else
                    Hk_squared(n,k)=nan;
                end
            end
        else
            if imag(poles(k))==0 % real pole
                Hk_squared(n,k)=Hk(n,k).^2; 
            else
                Hk_squared(n,k)=nan;
            end
        end
    end
end

% resize the TF (delete the duplicate columns due to the same pair of
% complex conjugated poles)
m=size(Hk_squared,2); k=1;
while k<=m
    if isnan(Hk_squared(:,k))
        Hk_squared(:,k)=[];
    end
    m=size(Hk_squared,2); 
    k=k+1;
end   

%% compute spectral LSP from TF
% freq specific term of the spectral LSP (total)
for n=1:nfft
    H(n)=prod(Hk(n,:)).^2; 
end
lspf(:,1)=0.5*log(H); 

% spectral LSP associated to a given oscillation (all poles)
% only the 
Q=size(Hk_squared,2);
LSPf=zeros(nfft,Q+1);
for k=1:Q
    lspf(:,k+1)=0.5*log(Hk_squared(:,k)); % single contribute - profile with zero mean
    LSPf(:,k+1) = LSP*ones(nfft,1) + lspf(:,k+1); % (sum of the single contribute + constant term)
    LSPf(:,1) = LSPf(:,1) + LSPf(:,k+1); % total term 
end

%% save the indices of the desired freq bands to be extracted from lspf
if range_all ~= 0
    for i=1:b
        range=range_all(i,:); % select the band
        pos=find(freq<range(2) & freq>range(1)); % find the poles in that band
        if isempty(pos)==0 % the pole exists in that band
            power=pot(pos);
            if length(power)>1 % if more than 1 pole is found 
                index=find(power==max(power)); % take the pole with maximum power
            else
                index=1;
            end
            % save the corresponding index
            index_band(i)=pos(index); % TF(:,pos(index))
        else % the pole does not exist in that band
            index_band(i)=nan;
            warn(i)=1;
        end
    end
end

%% return values in output
%%% TF 
out.H=H; % associated to the whole process
out.Hk=Hk_squared; % associated to the single contributes

%%% LSP
out.LSP=LSP; % time domain
out.lspf=lspf; % freq specific terms (no constant term added)
out.LSPf=LSPf; % total term (sum of the single contributes + constant term)

out.index_band=index_band; % index of the profiles to consider

%%% warning
out.warn=warn;

%%% freq axis
out.f=f;

end