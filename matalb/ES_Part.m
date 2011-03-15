function [f1,f2,f3,f4,f5,f6,f7]=ES_Part(filename,NumParticle,L)

%data = load(filename);
%t = data(1,1);
%NumParticle = 32;
%if (NumParticle == 32)
    %Data = data(:,1:6);
% else 
%     Data = data(:,3:8);    
% end

%N = -NumParticle/2 + 1:NumParticle/2;
%N=1:NumParticle;
%N = (1 - NumParticle) / 2 * dk : dk : (NumParticle - 1) / 2 * dk;

%kx = (2 * pi .* N  / L)';
%ky = (2 * pi .* N  / L)';
%kz = (2 * pi .* N  / L)';

% origin_kx = (1 + NumParticle) / 2;
% origin_ky = (1 + NumParticle) / 2;
% origin_kz = (1 + NumParticle) / 2;

% origin_kx = 0.0;
% origin_ky = 0.0;
% origin_kz = 0.0;


%x = reshape(Data(:,1),NumParticle,NumParticle,NumParticle);
%y = reshape(Data(:,2),NumParticle,NumParticle,NumParticle);
%z = reshape(Data(:,3),NumParticle,NumParticle,NumParticle);
 %k1 = 2 * pi / L * (x * (NumParticle - 1) / L - (NumParticle) /2);
 %k2 = 2 * pi / L * (y * (NumParticle - 1) / L - (NumParticle) /2);
 %k3 = 2 * pi / L * (z * (NumParticle - 1) / L - (NumParticle) /2);
% k1 =  2 * pi * (x * (NumParticle - 1) / L - (NumParticle + 1) /2);
% k2 =  2 * pi * (y * (NumParticle - 1) / L - (NumParticle + 1) /2);
% k3 =  2 * pi * (z * (NumParticle - 1) / L - (NumParticle + 1) /2);
% k1 = 2 * pi * (x * NumParticle / L + (1 - NumParticle) / 2);
% k2 = 2 * pi * (y * NumParticle / L + (1 - NumParticle) / 2);
% k3 = 2 * pi * (z * NumParticle / L + (1 - NumParticle) / 2);
% Ke = 0.5 * sum(sum(sum(u.^2 + v.^2 + w.^2)));
% 
% % u = fftshift(u);
% % v = fftshift(v);
% % w = fftshift(w);
data = load(filename);
Data = data(:,1:6);

u = reshape(Data(:,4),NumParticle,NumParticle,NumParticle);
v = reshape(Data(:,5),NumParticle,NumParticle,NumParticle);
w = reshape(Data(:,6),NumParticle,NumParticle,NumParticle);

dk = 2 * pi / L;
kx = zeros(NumParticle,NumParticle,NumParticle);
ky = zeros(NumParticle,NumParticle,NumParticle);
kz = zeros(NumParticle,NumParticle,NumParticle);
kk = zeros(NumParticle,NumParticle,NumParticle);

for l = 1:NumParticle,
    for m = 1:NumParticle,
        for n = 1:NumParticle,
            kx(n,m,l) = n - NumParticle/2 - 1/2;
            ky(n,m,l) = m - NumParticle/2 - 1/2;
            kz(n,m,l) = l - NumParticle/2 - 1/2;
            kk(n,m,l) = sqrt(kx(n,m,l)^2 + ky(n,m,l)^2 + kz(n,m,l)^2);
        end
    end
end

U = fftn(u);
V = fftn(v);
W = fftn(w);
U = fftshift(U);
V = fftshift(V);
W = fftshift(W);

[U_incomp, V_incomp, W_incomp, U_comp, V_comp, W_comp] = comp_incomp(kx,ky,kz,U,V,W);

UST = (U .* conj(U) + V .* conj(V) + W .* conj(W)) / ((NumParticle)^3);
UST_incomp = (U_incomp .* conj(U_incomp) + V_incomp .* conj(V_incomp) + W_incomp .* conj(W_incomp)) / ((NumParticle)^3);
UST_comp = (U_comp .* conj(U_comp) + V_comp .* conj(V_comp) + W_comp .* conj(W_comp)) / ((NumParticle)^3);


 E = zeros(1,NumParticle*2);
 E_comp = zeros(1,NumParticle*2);
 E_incomp = zeros(1,NumParticle*2);
 index = 1;
 index_num = 0;
 index_comp_num = 0;
 index_incomp_num = 0;



for k=dk:dk:dk*NumParticle/2,
    index_num = 0;
    index_comp_num = 0;
    index_incomp_num = 0;
    for l = 1:NumParticle, 
        for m = 1:NumParticle,
            for n = 1:NumParticle,
                low_shell = k - 0.5*dk;
                high_shell = k + 0.5*dk;
                if ((kk(n,m,l)>= low_shell) && (kk(n,m,l)< high_shell))
                    E(index) = E(index) + 0.5 * UST(n,m,l);
                    E_incomp(index) = E_incomp(index) + 0.5 * UST_incomp(n,m,l);
                    E_comp(index) = E_comp(index) + 0.5 * UST_comp(n,m,l);
                    index_num  = index_num + 1;
                    index_incomp_num = index_incomp_num + 1;
                    index_comp_num = index_comp_num + 1;
                end
            end
        end
    end
    N_num(index) = index_num ;
    N_comp_num(index) = index_comp_num;
    N_incomp_num(index) = index_incomp_num;
    %Ek(index) = 4 * pi * k * k * E(index) / N_num(index);
    %if (index==1) Ek(index) = E(index);
    %end
    index = index + 1;
end

k=dk:dk:dk*NumParticle/2;
E = E(1:length(k)) ./ N_num(1:length(k));
E_comp = E_comp(1:length(k)) ./ N_comp_num(1:length(k));
E_incomp = E_incomp(1:length(k)) ./ N_incomp_num(1:length(k));
Ek = 4 * pi * k .* k .* E  ;
Ek_comp = 4 * pi * k .* k .* E_comp  ;
Ek_incomp = 4 * pi * k .* k .* E_incomp;
KK = 0.5*sum(UST(:));

%nnnnn_com = sum(sum(N_com_num))
%nnnnn_incom = sum(sum(N_incom_num))
%nnnnn_total = sum(sum(N_num))

f1 = k;
f2 = E;
f3 = Ek;
f4 = E_comp;
f5 = Ek_comp;
f6 = E_incomp;
f7 = Ek_incomp;
end