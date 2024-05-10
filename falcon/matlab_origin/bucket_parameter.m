function [psi1,psi2,bucket_height,trapping_fraction,bucket_area,bunching,sinaverage]=bucket_parameter(psir)
% calculate the bucket parameters with psir

psi1 = pi-psir; % the largr endpoint of the separatrix
pond = @(psi,psir) cos(psi1)+psi1*sin(psir)-(cos(psi)+psi*sin(psir)); %gamma_{psi}=gamma_{r}
psi2 = fsolve(@(psi) pond(psi,psir), [-pi:psi1*0.99]); % solution of the equation:gamma_{psi}=gamma_{r}
psi2 = psi2(1); % the smaller endpoint of the separatrix
bucket_height = sqrt(cos(psir)-(pi/2-psir)*sin(psir)); % bucket_height
bucket_area = (1-sin(psir))/(1+sin(psir)); % bucket area approx. by S.Y.Lee
trapping_fraction = (psi1-psi2)/2/pi; % trapping fraction
potential = @(psi,psir) cos(psi)+psi*sin(psir); % ponderomotive potential equation
bucket_sep = @(psi) sqrt((potential(psi,psir)-potential(psi1,psir))/2);
bucket_bun = @(psi) bucket_sep(psi).*exp(1i*psi);
bucket_sin = @(psi) bucket_sep(psi).*sin(psi);
bunching = abs(integral(bucket_bun,psi2,psi1)/integral(bucket_sep,psi2,psi1));
sinaverage = abs(integral(bucket_sin,psi2,psi1)/integral(bucket_sep,psi2,psi1));
end





