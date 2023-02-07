function [Rt RaExp RaPeak RmemExp RmemPeak CmExp CmPeak] =  basicParams_calc(TimeChan,CurrChan,VoltChan,Disp)
%sampling interval
 SI = TimeChan(end,1)/length(TimeChan);

%num traces
NumTraces = size(CurrChan,2);

%find where the step occurs
PreVolt = mean(VoltChan(10:20,1));
LowVoltChan = VoltChan(:,1) - PreVolt;
LowIndices = find(LowVoltChan < -5);

StartStepIndex = LowIndices(2);
EndStepIndex = LowIndices(end);

%find number of neg or pos steps
NegSteps = [];
PosSteps = [];
for yy = 1:NumTraces
    if PreVolt > VoltChan(StartStepIndex+10,yy) + 0.5;
        NegSteps(end+1) = yy;
    elseif PreVolt < VoltChan(StartStepIndex+10,yy) - 0.5;
        PosSteps(end+1) = yy;
    end
end

% ivals are two values spanning the initial baseline (-60mV holding)
ivals1 = [300 300+round(0.02/SI)];
ivals2 = [800 900];%[EndStepIndex-round(0.05/SI) EndStepIndex];
% fvals are two values spanning the initial part of the voltage step but after the capacitive transient
fvals1 = [800 900]; %[StartStepIndex+round(0.06/SI) StartStepIndex+round(0.11/SI)];
fvals2 = [1600 1800];%[length(TimeChan(:,1))-round(0.05/SI) length(TimeChan(:,1))-10];

%% calculate Rt, Ra_exp, Ra_peak, Rm_exp, Rm_peak, Cm_exp, Cm_peak


Rtmat = [];
RaExpmat = [];
RaPeakmat = [];
RmemExpmat = [];
RmemPeakmat = [];
CmExpmat = [];
CmPeakmat = [];

for swp = NegSteps
    %first transient
    MinIndex = find(CurrChan(StartStepIndex:StartStepIndex+round(0.075/SI),swp) == min(CurrChan(StartStepIndex:StartStepIndex+round(0.075/SI),swp)));
    PeakIndex = StartStepIndex+MinIndex(1)-1;
    edata = CurrChan(PeakIndex:PeakIndex+round(0.075/SI),swp) - mean(CurrChan(PeakIndex+round(0.075/SI):PeakIndex+round(0.085/SI),swp));
    
    
    
    
    ftype = fittype('a*exp(b*x)+c','independent','x','coefficients',{'a','b','c'});
    fopt = fitoptions('Method','NonLinearLeastSquares','Lower',[-Inf -Inf -Inf],'Upper',[Inf Inf Inf], 'StartPoint', [-400 -200 0]);
    
    
    t = TimeChan(PeakIndex:PeakIndex+round(0.075/SI),swp) - TimeChan(PeakIndex,swp);
    f = fit(t,edata(:,1),ftype,fopt);
    expfit = f.a*exp(f.b*t)+f.c;
    expfitcalc = expfit - f.c;
    
    tau = (-1)*(1/f.b); % the time constant
    
    q1 = SI*sum(expfitcalc)/10^12; % q1 in coulombs
    
    
    deltaI = mean((CurrChan(fvals1(1):fvals1(2),swp)) - mean(CurrChan(ivals1(1):ivals1(2),swp)))/(10^12);
    deltaV = mean((VoltChan(fvals1(1):fvals1(2),swp)) - mean(VoltChan(ivals1(1):ivals1(2),swp)))/(10^3);
    
    Rt1 = (deltaV/deltaI); % total resistance in ohms
    
    q2 = tau*deltaI; % q2 in coulombs
    qt = q1 + q2;
    
    RaExp1 = tau*deltaV/qt; % access resistance in ohms calculated using the exponential
    RmemExp1 = Rt1-RaExp1; % membrane resistance
    CmExp1 = (qt*Rt1)/(deltaV*RmemExp1); % membrane capacitance in farads
    
    %alternative Ra calculation
    delI = (CurrChan(PeakIndex,swp) - median(CurrChan(ivals1(1):ivals1(2),swp)))/(10^12);
    RaPeak1 = deltaV/delI; % access resistance in ohms claculated using the peak
    RmemPeak1 = Rt1-RaPeak1; % membrane resistance
    
    Rth1 = (RaPeak1*RmemPeak1)/(RaPeak1+RmemPeak1); 
    
    CmPeak1 = tau/Rth1; %(qt*Rt1)/(deltaV*RmemPeak1); % membrane capacitance in farads
    



    if Disp
        figure()
        plot(t,edata(:,1));
        hold on
        plot(t,expfit,'r');
        pause
    end
    
    Rtmat(end+1) = Rt1;
    RaExpmat(end+1) = RaExp1;
    RaPeakmat(end+1) = RaPeak1;
    RmemExpmat(end+1) = RmemExp1;
    RmemPeakmat(end+1) = RmemPeak1;
    CmExpmat(end+1) = CmExp1;
    CmPeakmat(end+1) = CmPeak1;
    
    
    
    %2nd transient
    MaxIndex = find(CurrChan(EndStepIndex:EndStepIndex+round(0.075/SI),swp) == max(CurrChan(EndStepIndex:EndStepIndex+round(0.075/SI),swp)));
    PeakIndex = EndStepIndex+MaxIndex(1)-1;
    edata = CurrChan(PeakIndex:PeakIndex+round(0.075/SI),swp) - mean(CurrChan(PeakIndex+round(0.075/SI):PeakIndex+round(0.085/SI),swp));
    
    
    
    ftype = fittype('a*exp(b*x)+c','independent','x','coefficients',{'a','b','c'});
    fopt = fitoptions('Method','NonLinearLeastSquares','Lower',[-Inf -Inf -Inf],'Upper',[Inf Inf Inf],'StartPoint',[400 -200 0]);
    
    
    t = TimeChan(PeakIndex:PeakIndex+round(0.075/SI),swp) - TimeChan(PeakIndex,swp);
    f = fit(t,edata(:,1),ftype,fopt);
    expfit = f.a*exp(f.b*t)+f.c;
    expfitcalc = expfit - f.c;
    tau = (-1)*(1/f.b); % the time constant
    
    q1 = SI*sum(expfitcalc)/10^12; % q1 in coulombs
    
    
    deltaI = mean((CurrChan(fvals2(1):fvals2(2),swp)) - mean(CurrChan(ivals2(1):ivals2(2),swp)))/(10^12);
    deltaV = mean((VoltChan(fvals2(1):fvals2(2),swp)) - mean(VoltChan(ivals2(1):ivals2(2),swp)))/(10^3);
    
    Rt2 = (deltaV/deltaI); % total resistance in ohms
    
    q2 = tau*deltaI; % q2 in coulombs
    qt = q1 + q2;
    
    RaExp2 = tau*deltaV/qt; % access resistance in ohms calculated using the exponential
    RmemExp2 = Rt2-RaExp2; % membrane resistance
    CmExp2 = (qt*Rt2)/(deltaV*RmemExp2); % membrane capacitance in farads
    
    %alternative Ra calculation
    delI = (CurrChan(PeakIndex,swp) - median(CurrChan(ivals2(1):ivals2(2),swp)))/(10^12);
    RaPeak2 = deltaV/delI; % access resistance in ohms claculated using the peak
    RmemPeak2 = Rt2-RaPeak2; % membrane resistance
    
    Rth2 = (RaPeak2*RmemPeak2)/(RaPeak2+RmemPeak2); 
    
    CmPeak2 = tau/Rth2; %(qt*Rt2)/(deltaV*RmemPeak2); % membrane capacitance in farads




    if Disp
        figure()
        plot(t,edata(:,1));
        hold on
        plot(t,expfit,'r');
        pause
    end
    
    Rtmat(end+1) = Rt2;
    RaExpmat(end+1) = RaExp2;
    RaPeakmat(end+1) = RaPeak2;
    RmemExpmat(end+1) = RmemExp2;
    RmemPeakmat(end+1) = RmemPeak2;
    CmExpmat(end+1) = CmExp2;
    CmPeakmat(end+1) = CmPeak2;
    
end

Rt = median(Rtmat)/1e6; %MOhms
RaExp = median(RaExpmat)/1e6;
RaPeak = median(RaPeakmat)/1e6;
RmemExp = median(RmemExpmat)/1e6;
RmemPeak = median(RmemPeakmat)/1e6;
CmExp = median(CmExpmat)*1e12; %pF
CmPeak = median(CmPeakmat)*1e12; 

end
