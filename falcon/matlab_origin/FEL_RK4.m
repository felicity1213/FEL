function [newphasespace,newevalue]= FEL_RK4(phasespace,evalue,parameter,awvalue)

gammar_sq = parameter.lambdau/(2*parameter.lambda0)*(1+awvalue^2);

% RK-4 for intergration
k1theta = parameter.stepsize*(parameter.ku*(1-(gammar_sq./phasespace(:,2).^2)));
k1gamma=parameter.stepsize*(parameter.chi2*(awvalue./phasespace(:,2)).*...
    real(evalue*exp(1i*phasespace(:,1))));

k2theta=parameter.stepsize*(parameter.ku*(1-(gammar_sq./(phasespace(:,2)+0.5*k1gamma).^2)));
k2gamma=parameter.stepsize*(parameter.chi2*(awvalue./(phasespace(:,2)+0.5*k1gamma)).*...
    real(evalue*exp(1i*(phasespace(:,1)+0.5*k1theta))));

k3theta=parameter.stepsize*(parameter.ku*(1-(gammar_sq./(phasespace(:,2)+0.5*k2gamma).^2)));
k3gamma=parameter.stepsize*(parameter.chi2*(awvalue./(phasespace(:,2)+0.5*k2gamma)).*...
    real(evalue*exp(1i*(phasespace(:,1)+0.5*k2theta))));

k4theta=parameter.stepsize*(parameter.ku*(1-(gammar_sq./(phasespace(:,2)+k3gamma).^2)));
k4gamma=parameter.stepsize*(parameter.chi2*(awvalue./(phasespace(:,2)+k3gamma)).*...
    real(evalue*exp(1i*(phasespace(:,1)+k3theta))));

newphasespace(:,1)=phasespace(:,1)+1/6*(k1theta+2*k2theta+2*k3theta+k4theta);
newphasespace(:,2)=phasespace(:,2)+1/6*(k1gamma+2*k2gamma+2*k3gamma+k4gamma);

% Predictor-Corrector method for the field

f1=parameter.chi1*awvalue*...
     mean(exp(-1i*phasespace(:,1))./phasespace(:,2));
f1star=parameter.chi1*awvalue*...
     mean(exp(-1i*newphasespace(:,1))./newphasespace(:,2));

newevalue=evalue-parameter.stepsize/2.*(f1+f1star);





