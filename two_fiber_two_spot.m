
tlife=26*10^-9;
gamma=1/tlife;%%%%divided by two, assume Rabi frequency is very high and take average for thr signal
 %% calculate
 tic
 dt=0.3*10^-9;% time step
 t0=0.2*10^-6; %%%start from 0.2us
 treal=0;
 x0=0;
 v=300;

 
 %%%%%%%%%%initial 30um laser then 90um space then 30um laser
 field=100*10^-6;
 sep_dis=500*10^-6;
 timeout=0.02;
 transit_t=(field)/v;
 numberofatoms=round(timeout/transit_t*0.2);%%%%%on average one atom at the collecting volume
 
 
 signal=zeros(numberofatoms,2);
 signalindex=1;
 signal2=zeros(numberofatoms,2);
 signal2index=1;
 
 omega=2.8*gamma;
 
 for i=1:numberofatoms
     atomin=rand*timeout;
     v=abs(normrnd(100,30));
     %v=300;
     x=x0;
     tatom=0;
     tatom_emission=0;
     enter=0;
     while x<field+sep_dis+field+1*10^-6
        if (x>0 && x<field) || (x>field+sep_dis && x<field+sep_dis+field)
        prob=dt*gamma*1/2*(1-cos(omega*(tatom-tatom_emission)));%%equivalent to monte carlo wave simulation, 
     %%%%%%%%the jump probability is dropping according to the master
     %%%%%%%%equation which agrees with monte carlo wave.
     %prob=dt*initialscat;
        if rand<=prob
            tatom_emission=tatom;
            coupprob=1;%%%calculate the coupling efficiency at this x
            %%%%%the first fiber for the first photon
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
         end
         end
     tatom=tatom+dt;
     x=x+v*dt;
     if x>field+sep_dis && enter==0
         tatom_emission=tatom;
         enter=1;
     end
     end
 end
 %%%delete the zeros and sortrows
 signal = signal(any(signal,2),:);
 signal=sortrows(signal);
 signal2 = signal2(any(signal2,2),:);
 signal2=sortrows(signal2);
 
 
 plot(signal(:,1),signal(:,2),'*',signal2(:,1),signal2(:,2),'.') 
 ylim([-0.1,1.1])
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

 
 
 
 
 
 