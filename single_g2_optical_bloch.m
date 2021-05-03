
tlife=26.2347*10^-9;
gamma=1/tlife;%%%%divided by two, assume Rabi frequency is very high and take average for thr signal
 %% calculate
 tic
 dt=0.3*10^-9;% time step
 t0=0.2*10^-6; %%%start from 0.2us
 treal=0;
 x0=-2*10^-6;
 v=300;
 noise=300;
 
 %%%%%%%%%%initial 30um laser then 90um space then 30um laser
 field=20*10^-6;
 timeout=0.1;
 transit_t=(field)/300;
 numberofatoms=round(timeout/transit_t*0.1);%%%%%this number sets how many atoms at the collecting volume
 
 photons=zeros(numberofatoms,1);
 signal=zeros(numberofatoms,2);
 signalindex=1;
 signal2=zeros(numberofatoms,2);
 signal2index=1;
 
 omega0=2.2*gamma;
 %omega0=100*2*pi*10^6;
 %numberofatoms=1;
 %omega1=50*2*pi*10^6;
 
 for i=1:numberofatoms
     atomin=rand*timeout;
     %atomin=timeout*i/numberofatoms; %%%%this code change to pulse
     %excitation of quantum dot
     v=abs(normrnd(300,90));
     %v=300;
     pee=0;
     pgg=1;
     peg_pge=0;
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
     while x<field+0.5*10^-6
        omega=omega0*exp(-(x-field/2)^2/(62*10^-6)^2);
        
        %%%%size of the beam is w=62um
        %prob=dt*gamma*1/2*(1-cos(omega*(tatom-tatom_emission)));
        prob=dt*gamma*pee;
        %%equivalent to monte carlo wave simulation, 
     %%%%%%%%the jump probability is dropping according to the master
     %%%%%%%%equation which agrees with monte carlo wave.
     %prob=dt*initialscat;
        if rand<=prob
            tatom_emission=tatom;
            pee=0;
            pgg=1;
            peg_pge=0;
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
         end
     tatom=tatom+dt;
     x=x+v*dt;
     pee=pee+dt*(-gamma*pee+1i*omega/2*peg_pge);
     pgg=pgg+dt*(gamma*pee-1i*omega/2*peg_pge);
     peg_pge=peg_pge+dt*(-gamma/2*peg_pge+1i*omega*(pee-pgg));%%%transverse relaxation is set to be gamma/2
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
fileID = fopen('fluctuation1.txt','w');
fprintf(fileID,'%.12f\n',signal(:,1)');
fclose(fileID);

fileID = fopen('fluctuation2.txt','w');
fprintf(fileID,'%.12f\n',signal2(:,1)');
fclose(fileID);

 toc

 
 
 
 
 
 