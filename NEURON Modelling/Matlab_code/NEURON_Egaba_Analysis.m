%%This script analyzes NEURON data output after running 'mosinit.hoc'%%
% The script must be run from within the 'data' directory containing the
% dat files output by running the NEURON simulation
% please ensure the functions in the 'dependencies' folder are in the path
% Joseph Raimondo 2/2/2023

GABAe_low = -85;
GABAe_high = -35;
GABAe_step = 5;
Numsimsteps = (GABAe_high - GABAe_low)/GABAe_step;
SetEgabas = [GABAe_low:GABAe_step:GABAe_high-GABAe_step];

display = 0;
series_r_factor = 0.9;

baserampInd = [2930 3325];
pufframpInd = [5450 5845];

prefix = 'puff';
%puffdists = {'0' '50' '100' '150' '200'};
puffdists = {'spread'};
puffdistslabel = {};
for ll = 1:size(puffdists,2)

    for jj = 1:Numsimsteps
        vvec = textread([prefix '_' puffdists{ll} '_v_r00' num2str(jj-1) '.dat']);
        vmemvec = textread([prefix '_' puffdists{ll} '_vmem_r00' num2str(jj-1) '.dat']);
        ivec = textread([prefix '_' puffdists{ll} '_i_r00' num2str(jj-1) '.dat']);
        time = textread([prefix '_' puffdists{ll} '_t_r00' num2str(jj-1) '.dat']);
        ivec = ivec*1000; %convert to pA


        %calculating basic parameters
        [Rt(jj) RaExp(jj) RaPeak(jj) RmemExp(jj) RmemPeak(jj) CmExp(jj) CmPeak(jj)] = basicParams_calc(time(1:2000)/1000,ivec(1:2000),vvec(1:2000),display);

        %with Series R correction
        [rev(jj) inputres(jj) memrev(jj) cond(jj)] = Sz_ramp_calc(ivec,vvec,RaPeak(jj)*1e6,display,series_r_factor);

        %with no Series R correction
        [revnoRs(jj) inputresnoRs(jj) memrevnoRs(jj) condnoRs(jj)] = Sz_ramp_calc(ivec,vvec,RaPeak(jj)*1e6,display,0);

        %using true mem voltage
        %%[revcorr(jj) inputrescorr(jj) memrevcorr(jj) condcorr(jj)] = Sz_ramp_calc(ivec,vmemvec,0,display,0);

        figure(1)
        subplot(3,1,1)
        hold on
        plot(time,ivec);
        ylabel('Current')
        subplot(3,1,2)
        hold on
        plot(time,vmemvec);
        ylabel('Soma true Vm')
        subplot(3,1,3)
        plot(time,vvec);
        ylabel('Command Vm')
        %xlim([0 1500])
        %ylim([-1500 500])

    end


    %actual Egaba on x axis
    figure()
    plot(SetEgabas,rev,'-o')
    xlabel('Actual Egaba')
    ylabel('Measured Egaba')
    hold on
    plot(SetEgabas,revnoRs,'-ro')
    %plot(SetEgabas,revcorr,'-go')
    plot(SetEgabas,SetEgabas,'-k')
    ylim([-100 50])
    xlim([-85 -40])
    legend('90% Rs correction','0% Rs correction')
    title(['NBQX: Dist from Soma = ' puffdists{ll} ' um'])
    puffdistslabel{ll} = ['Dist from Soma = ' puffdists{ll} ' um'];

    %measured Egaba on x axis
    figure()
    plot(rev,SetEgabas,'-o')
    ylabel('Actual Egaba')
    xlabel('Measured Egaba')
    hold on
    plot(revnoRs,SetEgabas,'-ro')
    %plot(SetEgabas,revcorr,'-go')
    plot(SetEgabas,SetEgabas,'-k')
    %xlim([-100 50])
    %ylim([-85 -40])
    %line([-81.1 -81.1],[-85 -40],'Color','red','LineStyle','--') %anaethetised
    %line([-63.3 -63.3],[-85 -40],'Color','red','LineStyle','--') %awake
    line([-78.1 -78.1],[-85 -40],'Color','red','LineStyle','--') %NBQX

    legend('90% Rs correction','0% Rs correction')
    title(['NBQX: Dist from Soma = ' puffdists{ll} ' um'])
    puffdistslabel{ll} = ['Dist from Soma = ' puffdists{ll} ' um'];

    %%work out the errors
    p = polyfit(rev,SetEgabas,4);
    testEgabas = [-85:40];
    errorvals(:,ll) = polyval(p,testEgabas)-testEgabas;


    disp(['Ra: ' num2str(mean(RaPeak))]);
    disp(['Rmem: ' num2str(mean(RmemPeak))]);
    disp(['Rt: ' num2str(mean(Rt))]);
    disp(['Cmem: ' num2str(mean(CmPeak))]);
    disp(['GABA Cond: ' num2str(mean(cond))]);

end
c = colormap(turbo(size(puffdists,2)));
figure()
for jj = 1:size(puffdists,2)
    plot(testEgabas,errorvals(:,jj),'Color',c(jj,:))
    hold on
end
%legend(puffdistslabel);
xlabel('Measured Egaba')
ylabel('Error in Egaba estimate (mV)')
ylim([-10 10])
xlim([-85 -40])

%save('Awake','testEgabas','errorvals')