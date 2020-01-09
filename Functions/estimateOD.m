function [alpha_est,mixture_ODest,Cfimi_ODest,Ylipolytica_ODest,d] = estimateOD(Cfimi_ODavg,Ylipolytica_ODavg,mixture_ODavg)

wavelength = 300:50:1000;
N_C = size(Cfimi_ODavg,2); % number of C. fimi samples
N_Y = size(Ylipolytica_ODavg,2); % number of Y. lipolytica samples

d_alpha = 0.01;
a = 0:d_alpha:1; % estimated alpha

%% Mixture

% Calculate Euclidean Distance
d = zeros(N_C,N_Y,numel(a));
% d = zeros(N_C,N_Y);
for n_y = 1:N_Y % Y. lipolytica
    for n_c = 1:N_C % C. fimi
        Mest = Cfimi_ODavg(:,n_c)*a + Ylipolytica_ODavg(:,n_y)*(1-a);
        d(n_c,n_y,:) = sqrt(nansum(bsxfun(@minus,mixture_ODavg,Mest).^2)); % sqrt(sum(mixture_ODavg - Mest))
%         Mest = Cfimi_ODavg(:,n_c) + Ylipolytica_ODavg(:,n_y);
%         d(n_c,n_y) = sqrt(nansum(bsxfun(@minus,mixture_ODavg,Mest).^2)); % sqrt(sum(mixture_ODavg - Mest))
    end
end

% Estimate Alpha and ODs
[Cfimi_ODest_idx,Ylipolytica_ODest_idx,a_est_idx] = ind2sub(size(d),find(d == min(min(min(d)))));
alpha_est = a(a_est_idx); % alpha estimate
Cfimi_ODest = alpha_est*Cfimi_ODavg(:,Cfimi_ODest_idx);
Ylipolytica_ODest = (1-alpha_est)*Ylipolytica_ODavg(:,Ylipolytica_ODest_idx);
mixture_ODest = Cfimi_ODest + Ylipolytica_ODest;

% Estimate ODs
% [Cfimi_ODest_idx,Ylipolytica_ODest_idx] = find(d == min(d(:)));
% Cfimi_ODest = Cfimi_ODavg(:,Cfimi_ODest_idx);
% Ylipolytica_ODest = Ylipolytica_ODavg(:,Ylipolytica_ODest_idx);
% mixture_ODest = Cfimi_ODest + Ylipolytica_ODest;

end

