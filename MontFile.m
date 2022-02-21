% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
function [ mcbias ] = MontFile( rss )

num1=21; % 21 for 0, 5, 10, ..., 100 percentiles; 11 for 0, 10, 20, ..., 100 percentiles
num2=18; % number of parameters
yr=[2025:5:2090];
num3=num1*num2;
elasmu=1.45; % elasticity of marginal utility of consumption

for constrains=1:2

t=[5:5:25]; 
for ri=1:3
    r=rss(ri);
    duds=zeros(num3+100000,10,10);
    for mc=1:(num3+100000)
        load(strcat('D:\monte carlo',num2str(constrains),'\output_util-mc',num2str(mc),'.dat'),'-mat');
        for s=1:10
            for p=1:10
                A=zeros(1,5);
                A(1:5)=output_util(r,s:(s+4),p);
                [rs,ms,bs] = regression(t,A);
                duds(mc,s,p)=ms/(1-r*0.001)^(s*5); % discounting to the mitigation year
            end
        end
    end
    if r<100
        save(strcat('output\montecarlo',num2str(constrains),'_duds_rou0',num2str(r),'.dat'),'duds');
    end
end

t=[5:5:25];
duds2d=zeros(num3+100000,10,50);
duds2d_dk=zeros(num3+100000,10,5); % 6 for dk=2%, ... 6%
% gt=[0.00938	0.00994	0.01017	0.01004	0.00975	0.00924	0.00848	0.00769	0.00694	0.00609];
% gt=[0.01634	0.01606	0.01555	0.01476	0.01391	0.01289	0.01169	0.01051	0.00942	0.00827];

for mc=1:(num3+100000)
%     display(mc);
    load(strcat('D:\monte carlo',num2str(constrains),'\output_util-mc',num2str(mc),'.dat'),'-mat');
    for s=1:10
        for r=1:50
            A=zeros(1,5);
            A(1:5)=output_util(r,s:(s+4),10);
            [rs,ms,bs] = regression(t,A);
            duds2d(mc,s,r)=ms/(1-r*0.001)^(s*5); % discounting to the mitigation year
        end
    end
    for s=1:10
        for dk=1:5
            rou=0.01+dk*0.01-gt(s)*1.45;
            duds2d_dk(mc,s,dk)=duds2d(mc,s,max(1,min(50,round(rou*100))));
        end
    end
end
save(strcat('output\montecarlo',num2str(constrains),'_duds_rou_s.dat'),'duds2d');

mcbias=zeros(num3+100000,30);
obserror=[1	0.1	2	0.4	0.45	0.1	0.1	0.1	0.01	1	1	1	0.04	0.01	0.04];
for mc=1:(num3+100000)
    load(strcat('D:\monte carlo',num2str(constrains),'\output_obsm-mc',num2str(mc),'.dat'),'-mat');
    output_obsm(1:45,7)=smooth(output_obsm(1:45,7),5); % LandC obs
    output_obsm(1:45,17)=smooth(output_obsm(1:45,17),5); % Energy price obs
    A=output_obsm(1:28,25:30); output_obsm(:,25:30)=0;
    output_obsm(1:11,25)=A(18:28,1); % Energy model
    output_obsm(1:11,26)=A(18:28,2)./mean(A(5:16,2),1); % Energy model
    output_obsm(1:11,27)=A(18:28,3); % Energy model
    output_obsm(1:11,28)=A(18:28,4)./mean(A(5:16,4),1); % GDP model
    output_obsm(1:11,29)=A(18:28,5); % Energy price model
    output_obsm(1:11,30)=A(18:28,6)./mean(A(5:16,6),1); % Energy price model
    %
    devmc=zeros(15,1);
    devmc(1)=sqrt(sum((output_obsm(:,1)-output_obsm(:,2)).*(output_obsm(:,1)-output_obsm(:,2)),1)); %  Damage    
    devmc(2)=max(abs(output_obsm(:,3)-output_obsm(:,4)),[],1); %  temperature
    devmc(3)=max(abs(output_obsm(:,5)-output_obsm(:,6)),[],1); %  Atmco2
    devmc(4)=sqrt(mean((output_obsm(10:43,7)-output_obsm(10:43,8)).*(output_obsm(10:43,7)-output_obsm(10:43,8)),1)); %  LandC
    devmc(5)=sqrt(mean((output_obsm(5:43,9)-output_obsm(5:43,10)).*(output_obsm(5:43,9)-output_obsm(5:43,10)),1)); %  OceanC
    devmc(6)=max(abs(output_obsm(:,12)./output_obsm(:,11)-1),[],1); %  Capital 1971-2015
    devmc(7)=max(abs(output_obsm(:,14)./output_obsm(:,13)-1),[],1); %  GDP 1971-2015
    devmc(8)=max(abs(output_obsm(:,16)./output_obsm(:,15)-1),[],1); %  Energy 1971-2015
    devmc(9)=sqrt(mean((output_obsm(1:43,17)-output_obsm(1:43,18)).*(output_obsm(1:43,17)-output_obsm(1:43,18)),1)); %  Energy price 1971-2015    
    devmc(10)=sqrt(mean((output_obsm(1:26,19)-output_obsm(1:26,22)).*(output_obsm(1:26,19)-output_obsm(1:26,22)),1)); %  EUE rate 1971-2015
    devmc(11)=sqrt(mean((output_obsm(1:26,20)-output_obsm(1:26,23)).*(output_obsm(1:26,20)-output_obsm(1:26,23)),1)); %  EPE rate 1971-2015
    devmc(12)=sqrt(mean((output_obsm(1:26,21)-output_obsm(1:26,24)).*(output_obsm(1:26,21)-output_obsm(1:26,24)),1)); %  ENE rate 1971-2015      
    devmc(13)=sqrt(mean((output_obsm(1:11,25)-output_obsm(1:11,26)).*(output_obsm(1:11,25)-output_obsm(1:11,26)),1)); %  E in COVID9
    devmc(14)=sqrt(mean((output_obsm(1:11,27)-output_obsm(1:11,28)).*(output_obsm(1:11,27)-output_obsm(1:11,28)),1)); %  Y in COVID9
    devmc(15)=sqrt(mean((output_obsm(2:5,29)-output_obsm(2:5,30)).*(output_obsm(2:5,29)-output_obsm(2:5,30)),1)); %  P in COVID9
    mcbias(mc,1:15)=devmc(1:15,1);    
    for obs=1:15
        if devmc(obs)>obserror(obs)
            mcbias(mc,obs+15)=mcbias(mc,obs+15)+1;
        end
    end
end
save(strcat('output\montecarlo',num2str(constrains),'_mcbias.dat'),'mcbias');

end

end

