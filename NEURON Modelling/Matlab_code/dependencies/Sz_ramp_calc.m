function [rev inputres memrev cond] = Sz_ramp_calc(input_current,input_volt,Ra,display,series_r_factor)
%this also corrects for series R


%correcting for seriesR
input_volt = input_volt - input_current*(Ra*series_r_factor/10^9);


base_inds = [2975 3325];
gaba_inds = [5495 5845];



retract_num = 150;

basecurrent = input_current(base_inds(1):base_inds(2));
gabacurrent = input_current(gaba_inds(1):gaba_inds(2)) - input_current(base_inds(1):base_inds(2)); 
%volt = input_volt(base_inds(1):base_inds(2));
volt = input_volt(gaba_inds(1):gaba_inds(2));
try
    inputres = inputResCalc(volt,basecurrent);
catch
    inputres = NaN;
end
try
    memrev = memRevCalc(volt,basecurrent);
catch
    memrev = NaN;
end

try
    cond = condCalc(volt,gabacurrent);
catch
    cond = NaN;
end




%working out reversal using actual crossing
revposi = find((gabacurrent < 2) & (gabacurrent > -2));
rev = mean(volt(revposi)); % the chloride reversal

if display && ~isnan(rev)
    subplot(2,1,1)
    plot(input_current);
    subplot(2,1,2)
    plot(volt,gabacurrent);
    hold on
    plot(rev,0,'r+','MarkerSize',16);
    hold off
end

%working out reversal using linear fit
if (isnan(rev) && (gabacurrent(end) < 0))
    rev = rampLinearFit(volt(end-retract_num:end),gabacurrent(end-retract_num:end));
    if display && ~isnan(rev)
        subplot(2,1,1)
        plot(input_current);
        subplot(2,1,2)
        plot(volt,gabacurrent);
        hold on 
        plot(rev,0,'r+','MarkerSize',6);
        plot(volt(end-retract_num:end),gabacurrent(end-retract_num:end),'g');
        hold off      
    end
elseif (isnan(rev) && (gabacurrent(1) > 0))
    rev = rampLinearFit(volt(50:150),gabacurrent(50:150));
end






