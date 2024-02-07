%% MONOVARIATE SPECTRAL DECOMPOSITION
% Baselli 1997

function [Sk,Sm,freq2,pot2,pot2n,poli,f]=lsp_spect_dec(A,C,varw,fc,nfft)

%% poles, residuals and powers
poli=roots(A);
p=size(poli,1);
q=length(C);

res=NaN*ones(p,1);
for k=1:p
    pk=poli(k);
    zeta=NaN*ones(q,1); zeta1=zeta;
    z=pk; for kk=1:q, zeta(kk)=z^(-(kk-1)); end, Cz=C*zeta;
    z1=pk^(-1); for kk=1:q, zeta1(kk)=z1^(-(kk-1)); end, Cz1=C*zeta1;
   
    tmp=poli;tmp(k)=[];
    res(k)=(Cz*varw*Cz1)/(pk*prod(pk-tmp)*prod(1/pk-poli));
end
pot=real(res);

poli2=[];pot2=[];res2=[];
for n=1:p % poles with frequency between 0 and 0.5 (exclude conjugates)
    polo=poli(n); potenza=pot(n); residuo=res(n);
    if imag(polo)>=0
        poli2=[poli2; polo];
        res2=[res2; residuo];
        if imag(polo)==0
            pot2=[pot2; potenza];
        else
            pot2=[pot2; 2*potenza];
        end
    end
end
p2=length(poli2);

%% spectral components
f = (0:nfft-1)*(fc/(2*nfft)); f=f';
z = exp(1i*2*pi*f/fc); % vector of complex exponentials

poli2_=1./poli2;
res2_=-res2;
Sk=NaN*ones(nfft,p2); 
for k=1:p2 % (A.9) - Baselli 1997
    pk=poli2(k);resk=res2(k);
    pk_=poli2_(k);resk_=res2_(k);
    Sk(:,k)=(resk*pk./(z-pk)+resk_*pk_./(z-pk_));
    if imag(poli2(k))==0
        Sk(:,k)=real(Sk(:,k));
    else
        Sk(:,k)=2*real(Sk(:,k));
    end
end

freq2=angle(poli2).*(fc/(2*pi));
pot2n=100*pot2./sum(pot2);
Sm=sum(Sk,2);


