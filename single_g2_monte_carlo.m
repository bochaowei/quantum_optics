
%%%%divided by two, assume Rabi frequency is very high and take average for thr signal
 %% calculate
 tic
 tlife=26.2347*10^-9;
 gamma=1/tlife;
 hbar=1.0545718*10^-34;
 dt=0.1*10^-9;%time step
 treal=0;
 x0=-2.5*10^-6;%% start from 0.5um away from fiber
 v=300;
 %%%%%%%%%%initial 30um laser then 90um space then 30um laser
 field=20*10^-6;
 timeout=0.05;
 omega0=4.6*gamma;
 transit_t=(field)/300;
 numberofatoms=round(timeout/transit_t*0.16);%%%%%this number sets how many atoms at the collecting volume
 photons=zeros(numberofatoms,1);
 signal=zeros(numberofatoms,2);
 signalindex=1;
 signal2=zeros(numberofatoms,2);
 signal2index=1;
 %omega0=100*2*pi*10^6;
 %numberofatoms=1;
 %omega1=50*2*pi*10^6;
 m=1.443*10^-25;
 k=1.38*10^-23;
 T=273.15+80;
 px=0:0.5:1000;
 p= 0.5*m^2/k^2/T^2*px.^3.*exp(-0.5*m/k/T.*px.^2);
 p2=4*pi*px.^2*(m/2/pi/k/T)^(1.5).*exp(-0.5*m/k/T.*px.^2);
 maxwell_v=randpdf(p,px,[numberofatoms,1]);
 vapor_v=randpdf(p2,px,[numberofatoms,1]);

 
 for i=1:numberofatoms
     atomin=rand*timeout;
     %atomin=timeout*i/numberofatoms; %%%%this code change to pulse
     %excitation of quantum dot
     %v=abs(normrnd(295,20))+5;
     v=maxwell_v(i);
     %v=300;
     ee=[1,0;0,0]; %excited state is (1;0)
     gg=[0,0;0,1];
     c=[0,0;1,0]; %atomic lowering operator
     Ha=0; %zero detuning
     Hdamp=-1i*hbar/2*gamma*(c'*c);
     psi=[0;1];
     x=x0;
     tatom=0;
     tatom_emission=0;
%     sublevel=rand;
%      if sublevel<0.2
%          omega0=omega1;
%      elseif sublevel<=0.6 && sublevel>=0.2
%          omega0=omega1*sqrt(4/15)/sqrt(3/10);
%      elseif sublevel>0.6  
%          omega0=omega1*sqrt(1/6)/sqrt(3/10);
%      end
     while x<field+0.2*10^-6
        omega=omega0*exp(-(x-field/2)^2/(62*10^-6)^2);
        Hl=hbar/2*omega*(c+c');
        H=Ha+Hl+Hdamp;
        %=expm(-1i*H*dt/hbar);
        U0=eye(2)-1i*H*dt/hbar;
        psi=U0*psi;
        nor=psi'*psi;
        prob=1-nor;
        %%equivalent to monte carlo wave simulation, 
     %%%%%%%%the jump probability is dropping according to the master
     %%%%%%%%equation which agrees with monte carlo wave.
     %prob=dt*initialscat;
        if rand<=prob
            tatom_emission=tatom;
            psi=[0;1];
            if (x>0 && x<field)
            photons(i)=photons(i)+1;
            coupprob=1;%%%calculate the coupling efficiency at this x
            %%%%%the first fiber for the first photon
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rand<=coupprob
                if rand<=0.5
                    signal(signalindex,1)=atomin+tatom;
                    signal(signalindex,2)=1;
                    signalindex=signalindex+1;
                else
                    signal2(signal2index,1)=atomin+tatom;
                    signal2(signal2index,2)=1;
                    signal2index=signal2index+1;  
                end
                
            end
            %%%%%%%%%%%%%%%%%%%
            end
        else
            psi=psi/sqrt(psi'*psi); 
        end
     tatom=tatom+dt;
     x=x+v*dt;
     end
 end
 %%%delete the zeros and sortrows
 signal = signal(any(signal,2),:);
 signal=sortrows(signal);
 signal2 = signal2(any(signal2,2),:);
 signal2=sortrows(signal2);
 
 
histogram(photons,'BinLimits',[0,26],'Normalization','probability')
xlabel('Number of photons per atom')
ylabel('Probability')
 %ylim([-0.1,1.1])
%  ch1=length(signal);
%  length(signal)/numberofatoms*100
%  ch2=length(signal2)
%  length(signal)/numberofatoms*100
%  figure(2)
% plot(signal2(:,1),signal2(:,2),'*') 
%  ylim([-0.1,1.1])
fileID = fopen('fluctuation1_monte_t.txt','w');
fprintf(fileID,'%.12f\n',signal(:,1)');
fclose(fileID);

fileID = fopen('fluctuation2_monte_t.txt','w');
fprintf(fileID,'%.12f\n',signal2(:,1)');
fclose(fileID);

 toc

 
 
 
 
 
 