workdir = '/path/to/the/data/folder';

strain_mag = 10^(-7); %magnitude of strain
num_omega = 15;  %number of different frequencies that are used in the simulation
dt = 0.01; %time step
filedir = workdir;

results_dir =[filedir,'-fig'];   %make a folder to save the analysis results
if exist(results_dir,'dir')
    sprintf('folder exists')
else
    mkdir(results_dir)
end


Gelastic = zeros(1,1);
omega_list = zeros(1,1);

for j = 0:num_omega-1
    %%%%%%%%%%%%%%%%%%%%
    %set up the sampling parameters
    %%%%%%%%%%%%%%%%%%%%%
    T = 2^(num_omega-3)/2^j;    %length of the period used in the simulation
    filename = [filedir,sprintf('/num_period_t_%04.8f.dat',T)];
    Num_period = importdata(filename);
    Num_period_trans = Num_period-min(10,Num_period-1); % number of periods during which the response is regarded as transient so thrown away
    Num_period_steady = Num_period - Num_period_trans;   % number of periods during which the response is regarded as being in steady state
    Num_samp = 25;   % number of sampling points in one period
    Num = Num_period*Num_samp;   %total number of sampling points in one time series data
    Num_samp_trans = Num_samp*Num_period_trans;  %number of sampling points during transient state
    Num_samp_steady = Num_samp*(Num_period-Num_period_trans);    %total number of sampling points of steady state

    omega = 2*pi/T;    %frequency of the shear
    omega_list(j+1) = omega;   %record the list of all frequencies used in the simulation
    time = linspace(0,T*(Num_period-Num_period_trans),Num_samp_steady+1);
    time = time(1:end-1);   %an array of the sampling points in time during steady state
    strain_time = strain_mag*sin(omega*time);
    stress_ls = zeros(Num_samp_steady,4);
    
    filename = [filedir,sprintf('/stress_t_%04.8f.dat',T)];
    m = importdata(filename);
    stress_ls = m.data(Num_samp_trans+1:end,2:5);   %read the response stress tensor
    

    stress_elastic = -(stress_ls(:,2));   %shear stress
    figure(1)
    plot(time,stress_elastic,'c*-')
    hold on
    plot(time,strain_time,'k-')
    legend({'stress elastic','applied strain'});
    xlabel('time')
    ylabel('\sigma')
    filename = [results_dir,'/stress_curve-',num2str(j)];    %plot the stress and strain curve in time
    saveas(gcf,filename);
    close all
    
%%%%%%%%%%%%%%%%%%%%%%%%
%perform Fourier transform on stress signal
%%%%%%%%%%%%%%%%%%%%%%%%%%

    stress_elastic_fft = fft(stress_elastic(1:end))/Num_samp_steady;    %Fourier transform of shear stress
    strain_fft=fft(strain_time(1:end))/Num_samp_steady;                 %Fourier transform of strain
    q = 2*pi/(T*(Num_period_steady))*(0:floor(Num_samp_steady/2));
    figure(2)
    plot(q,abs(stress_elastic_fft(1:floor(Num_samp_steady/2)+1)))  %plot the sress signal in frequency space
    hold on
    plot(q,abs(strain_fft(1:floor(Num_samp_steady/2)+1)),'r')  %plot the strain signal in frequency space
    Gelastic(j+1) = stress_elastic_fft(Num_period_steady+1)/strain_fft(Num_period_steady+1);   %compute the complex moduli
    filename = [results_dir,'/FFT-',num2str(j)];
    saveas(gcf,filename);
    close all
%     pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot storage and loss moduli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loglog(omega_list,abs(real(Gelastic)),'bo-')
hold on
loglog(omega_list,abs(imag(Gelastic)),'b*-')   

xlabel('\omega')
ylabel('modulus')
legend('storage','loss')
legend('location','south')
filename = [results_dir,'/modulus'];
saveas(gcf,filename);
close all


M = [omega_list;real(Gelastic);imag(Gelastic)];
csvwrite([results_dir,'/modulusVSfreq.csv'],M');


