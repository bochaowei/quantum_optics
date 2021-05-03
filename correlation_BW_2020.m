Atic
%% 
name='3.2 112 4uw s 1-0';
filename=[name,'.txt'];
mydata=load(filename);
ch0=mydata(:,1);ch1=mydata(:,2);
%%%%%%%%%%%%%% if ch0 and ch1 are not same count rate(which is usually the case, then 
%%%%we need to clean these zeros of the lesser channel.
ch0= ch0(any(ch0,2),:);
ch1= ch1(any(ch1,2),:);
ch_sum=zeros(length(ch0)+length(ch1),1);
ch_sum(1:length(ch0))=ch0;
ch_sum(length(ch0)+1:length(ch0)+length(ch1))=ch1;
ch_sum=sortrows(ch_sum);
% filename='fluctuation'
% ch0=load([filename,'1','.txt']);
% ch1=load([filename,'2','.txt']);


%%Remember the photon time data need to be in time increasing order 
%%for the algorithm to work.

%%
% h=histogram(ch1(ch_sum<500),10000);  %%%%the number of photons per 50 ms(500/10000)
% lambda=mean(h.Values)
% max_v=max(h.Values)
% h1=histogram(h.Values,'Normalization','pdf'); %the distribution of 50 ms
% h1.BinEdges=[-0.5:1:max_v-0.5];
% legend('data')
% hold on
% h2=histogram(poissrnd(lambda,10000),'Normalization','pdf','BinWidth',1,'FaceColor','green');
% h2.BinEdges=[-0.5:1:max_v-0.5];
% hold off
%% calculate
M=50;%%%%%%number of bins for g2
time=100*10^-9;%%%%% lag time range from -time to +time
deltatau=time/M*2;
nij=zeros(1,M);
n0=zeros(1,M);

n1=zeros(1,M);
taumin=zeros(1,M);
taumax=zeros(1,M);
indexl=ones(1,M);
indexh=ones(1,M);
for i=1:M
    taumin(i)=(i-M/2)*deltatau;
    taumax(i)=(i-M/2+1)*deltatau;
end   
T=max([ch0',ch1']);
tauave=(taumin+taumax)/2;
% calculate denominator
% for ti=1:length(ch0)
%     for k=1:M
%         if ch0(ti)<=T-tauave(k)
%             n0(k)=n0(k)+1;
%         end
%     end
% end
% 
% for ti=1:length(ch1)
%     for k=1:M
%         if ch1(ti)>=tauave(k)
%             n1(k)=n1(k)+1;
%         end
%     end
% end

%%%%% below is the improved method to calculate denominator
for k=1:M
    for ti=length(ch0):-1:1
        if ch0(ti)<=T-tauave(k)
        n0(k)=ti;
        break
        end
    end
end
for k=1:M
    for ti=1:length(ch1)
        if ch1(ti)>=tauave(k)
            n1(k)=length(ch1)-ti+1;
            break
        end
    end
end


% calculate coincidences
for ti=1:length(ch0)
    for k=1:M
        while(indexl(k)<length(ch1)&& ch1(indexl(k))<ch0(ti)+taumin(k))
            indexl(k)=indexl(k)+1;    
        end
        while(indexh(k)<length(ch1)&& ch1(indexh(k))<=ch0(ti)+taumax(k))
            indexh(k)=indexh(k)+1;
        end
        nij(k)=nij(k)+indexh(k)-indexl(k);
    end
end
g2=nij;%%%%%%%use this to see coincident count
%g2=nij.*(T-tauave)./n0./n1./deltatau; %%this is real g2
lagtime=(1:M)*deltatau-time;
%%% plot
%figure(1)
%g2=g2/max(g2);
plot(lagtime*10^9,g2,'c*')
%plot(lagtime*10^9,g2)
ylim([0,max(g2)*1.1])
%ylim([0,9.1])
xlabel('Time(ns)')
ylabel('g2')
%legend('pi transition')
%ylabel('coincident counts')
ax = gca;
ax.FontSize = 17;
ax.FontWeight='bold';
ax.LineWidth = 1;
toc