workdir = '/path/to/the/data/folder';

strain_mag = 10^(-7); %magnitude of strain
num_omega = 15;  %number of different frequencies that are used in the simulation

filedir = workdir;

fileID = fopen([filedir,'/num_period_list.out']);
formatSpec = '%d';
num_period_list = fscanf(fileID,formatSpec);  %read the number of periods of each frequency in the simulation

results_dir =[filedir,'-fig'];  %make a folder to save the analysis results
if exist(results_dir,'dir')
    sprintf('folder exists')
else
    mkdir(results_dir)
end

for j = 0:num_omega-1
    Num_period = num_period_list(j+1);  %number of periods
    Num_samp = 25;   % number of sampling points in one periods
    Num_period_perblock = 3; %divide the response into several blocks, each block contains this number of period
    Num = Num_period*Num_samp; %total number of sampling points
    stressxy_transient = zeros(1,floor(Num_period/Num_period_perblock));  %create an array to save the dominant mode of shear stress within each block
    T = 2^(num_omega-3)/2^j; %length of the period used in the simulation
    omega = 2*pi/T; %frequency of the shear
    time = linspace(0,T*Num_period,Num+1);
    time = time(1:end-1); %an array of the sampling points in time
    strain_time = strain_mag*sin(omega*time);  %an array of the sampling points in strain
    stress_ls = zeros(Num,4); %an array to save the stress tensor
    
    filename = [filedir,sprintf('/stress_t_%04.8f.dat',T)];
    m = importdata(filename);
    stress_ls = m.data(:,2:5); %read the stress tensor

    stress_response = -(stress_ls(:,2)); %shear stress



    for k = 1:floor(Num_period/Num_period_perblock)
        stressxy_fft=fft(stress_response(1+(k-1)*Num_period_perblock*Num_samp:k*Num_period_perblock*Num_samp))/(Num_period_perblock*Num_samp); %fourier transform of shear stress within the block
        strain_fft=fft(strain_time(1+(k-1)*Num_period_perblock*Num_samp:k*Num_period_perblock*Num_samp))/(Num_period_perblock*Num_samp); %fourier transfomr of shear strain within the block
        stressxy_transient(k) = stressxy_fft(Num_period_perblock+1); %dominant mode of shear stress within each block
%         plot(abs(strain_fft(1:floor(Num_period_perblock*Num_samp/2)+1)))
%         figure(3)
%         plot(abs(stressxy_fft(1:floor(Num_period_perblock*Num_samp/2)+1)))
%         pause;
    end
    plot(abs(stressxy_transient),'k')  %plot to see how the dominant mode of shear stress converges to a steady state
    title(['T=',num2str(T,'%.3f')])
    filename = [results_dir,'/T',num2str(j)];
    saveas(gcf,filename);
    disp(j)    
end

