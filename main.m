
clear all;
close all;
clc ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%
info.maxit = 5000;
info.model=1;
%                  = 0 pooled model without fixed effects (default, x may contain an intercept)
%                  = 1 spatial fixed effects (x may not contain an intercept)
%                  = 2 time period fixed effects (x may not contain an intercept)
%                  = 3 spatial and time period fixed effects (x may not contain an intercept)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxv=30;
kk=14;
Co=1; %W for Co2 and others
sep=1; % if co2 is added seperatly 0 or 1
print=1;
A=xlsread('cigarette.xls');
CHOSE=[1:kk];
CHOSE2=[1:kk];
CHOSE3= 3;  %%% 1-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FREM=cell(2*maxv+3,13,kk);
FREMM=cell(2*maxv+3,3,kk);
FREMMF=cell(2*(2*maxv+3),2,kk);
ZA=cell(11,3,kk);
ZAA=cell(11,3,size(CHOSE,2));
EFFTO=cell(maxv,6+1,kk);

W2=xlsread('nan.xls',1);
NA=size(W2,1);

mohsen1=1; % sesnsitivity anaylize in different sheets for each variables 1 to number of variables
for mohsen=1:mohsen1
    FREM=cell(2*maxv+3,13,kk);
    FREMM=cell(2*maxv+3,3,kk);
    FREMMF=cell(2*(2*maxv+3),2,kk);
    ZA=cell(11,3,kk);
    ZAA=cell(11,3,size(CHOSE,2));
    EFFTO=cell(maxv,6+1,kk);
    for ii=1:kk
        T=50; % number of time periods
        N=648; % number of regions
        
        filnam= {'1200.txt';'2000.txt';'2500.txt';'3000.txt';'4000.txt';'5000.txt';'6000.txt';'7000.txt';'8000.txt';'9000.txt';'10000.txt';'12000.txt';'14000.txt';'19000.txt'};
        W1=textread(char(filnam(1)));
        Wb=zeros(N,N);
        for p=1:size(W1,1)
            aa=W1(p,1);
            bb=W1(p,2);
            Wb(aa,bb)=1/W1(p,3);
        end
        Wb(:,W2)=[];
        Wb(W2,:)=[];
        W=normw(Wb); % function of LeSage
        
        filnam= {'1200.txt';'2000.txt';'2500.txt';'3000.txt';'4000.txt';'5000.txt';'6000.txt';'7000.txt';'8000.txt';'9000.txt';'10000.txt';'12000.txt';'14000.txt';'19000.txt'};
        W1c=textread(char(filnam(9)));
        Wcb=zeros(N,N);
        for p=1:size(W1c,1)
            aa=W1c(p,1);
            bb=W1c(p,2);
            Wcb(aa,bb)=1/W1c(p,3);
        end
        Wcb(:,W2)=[];
        Wcb(W2,:)=[];
        Wc=normw(Wcb); % function of LeSage
        
        
        if sep==1
            filnam= {'1200.txt';'2000.txt';'2500.txt';'3000.txt';'4000.txt';'5000.txt';'6000.txt';'7000.txt';'8000.txt';'9000.txt';'10000.txt';'12000.txt';'14000.txt';'19000.txt'};
            W1s=textread(char(filnam(13)));
            Wsb=zeros(N,N);
            for p=1:size(W1s,1)
                aa=W1s(p,1);
                bb=W1s(p,2);
                Wsb(aa,bb)=1;
            end
            Wsb(:,W2)=[];
            Wsb(W2,:)=[];
            Ws=normw(Wsb); % function of LeSage
        end
        N=N-NA;
        % dimensions of the problem
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % row-normalize W
        ms=size(A,2);
        AA=A(2:end,6:ms);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %namee = strvcat('tourist','constant','gdpp','pric','ope','govsta','SociCon','InvPro','IntCon','ExtCon','Corr','MilPol','RelTen','LawOr','EthTen','DemAcc','BurQua');
        y = A(2:end,5);
        
        ESSM=cell(1,2);
        LAT =    strvcat('in 50–60°S','in 40–50°S','in 30–40°S','in 20–30°S','in 10–20°S','in 0–10°S','in 0–10°N','in 10–20°N','in 20–30°N','in 30–40°N','in 40–50°N','in 50–60°N','in 60–70°N','in 70–80°N',' (130–30°W)',' (10–60°E)',' (60–160°E)');
        if mohsen==1
            ESSM{1,2}=['Deforestion ' LAT(ii,:)];
            ESSM{1,1}=['Total Other Deforestions '];
        end
        
        %namee =    strvcat('eco','constant','CO2','CO2-2','CO2-3','CO2-31','GDP','GDP2','UR','ENP','ENG','TR','POP','Built-up Land','Carbon','Cropland','Fishing Grounds','Forest Products','Grazing Land','co2g','WT','WT2','LOD','Sunspot','Thermosteric','Barystatic','GMSL');
        namee =    strvcat('tem','constant','in 50–60°S','in 40–50°S','in 30–40°S','in 20–30°S','in 10–20°S','in 0–10°S','in 0–10°N','in 10–20°N','in 20–30°N','in 30–40°N','in 40–50°N','in 50–60°N','in 60–70°N','in 70–80°N',strvcat(ESSM),'GWP','PRE','WINUn+','WINUn-','WINVn+','WINVn-','WINU+','WINU-','WINV+','WINV-','HCL','ML-CL','Sunspot','LOD','GMSL');
        CODE=zeros(1,size(namee,1));
        CAM=[1, 2];
        COL=[3];
        numexo=0;
        global numexo
        global print
        CODE(1, 1)=1;
        CODE(1, 2)=1;
        CODE(1, 17)=1;
        CODE(1, 18)=1;
        
        for jj=19:30
            CODE(1, jj)=1;
        end
        
        NV=sum(CODE)-2;
        CODE7=CODE;
        CODE7(1,1)=0;
        CODE7(1,2)=0;
        selectt=logical (CODE7);
        MSELmohsen=CODE(1, 3:size(CODE,2));
        selectxmohsen=logical (MSELmohsen);
        msmo=NV-numexo;
        for imm=1:size(MSELmohsen,2)
            if sum(MSELmohsen(1:imm))>msmo
                MSELmohsen(imm)=0;
            end
        end
        selectxmohsen2=logical (MSELmohsen);
        x=AA(:,selectxmohsen);
        x2=AA(:,selectxmohsen2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x(:,2)=x(:,1).*AA(:,ii);
        x2(:,2)=x2(:,1).*AA(:,ii);
        
        x(:,1)=x(:,1)-x(:,2);
        x2(:,1)=x2(:,1)-x2(:,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        namefi=namee(selectt,:);
        WM=repmat('W*',sum(CODE7)-numexo,1);
        vnamess=strvcat(namefi,[WM namefi(1:(NV-numexo),:)]);
        FREM{1,1,ii}='constant';
        if CHOSE3==1
            for i=1:NV
                FREM{i+1,1,ii}=vnamess(i,:);
            end
            FREM{NV+2,1,ii}='W*dep.var.';
        elseif CHOSE3==2
            for i=1:NV
                FREM{i+1,1,ii}=vnamess(i,:);
            end
            FREM{NV+2,1,ii}='spat.aut.';
        elseif CHOSE3==3
            for i=1:2*NV-numexo
                FREM{i+1,1,ii}=vnamess(i,:);
            end
            FREM{2*NV+2-numexo,1,ii}='W*dep.var.';
        elseif CHOSE3==4
            for i=1:NV
                FREM{i+1,1,ii}=vnamess(i,:);
            end
            FREM{NV+2,1,ii}='W*dep.var.';
            FREM{NV+3,1,ii}='teta';
        elseif CHOSE3==5
            for i=1:NV
                FREM{i+1,1,ii}=vnamess(i,:);
            end
            FREM{NV+2,1,ii}='spat.aut.';
            FREM{NV+3,1,ii}='teta';
        elseif CHOSE3==6
            for i=1:2*NV-numexo
                FREM{i+1,1,ii}=vnamess(i,:);
            end
            FREM{2*NV+2-numexo,1,ii}='W*dep.var.';
            FREM{2*NV+3-numexo,1,ii}='teta';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        wx=zeros(N*T,NV-numexo);
        wxc=zeros(N*T,NV-numexo);
        wxs=zeros(N*T,NV-numexo);
        
        for t=1:T
            t1=(t-1)*N+1;t2=t*N;
            wx(t1:t2,:)=(W*x2(t1:t2,:));
        end
        if sep==1
            for t=1:T
                t1=(t-1)*N+1;t2=t*N;
                wxc(t1:t2,:)=(Wc*x2(t1:t2,:));
            end
            wx(:,CAM)=wxc(:,CAM);
            for t=1:T
                t1=(t-1)*N+1;t2=t*N;
                wxs(t1:t2,:)=(Ws*x2(t1:t2,:));
            end
            wx(:,COL)=(wxs(:,COL));
        else
            for t=1:T
                t1=(t-1)*N+1;t2=t*N;
                wxc(t1:t2,:)=(Wc*x2(t1:t2,:));
            end
            wx(:,CAM)=wxc(:,CAM);
        end
        
        
        
        
        
        for t=1:T
            t1=(t-1)*N+1;t2=t*N;
            wym(t1:t2,:)=(W*y(t1:t2,:));
        end
        xconstant=ones(N*T,1);
        [nobs K]=size(x);
        
        % ----------------------------------------------------------------------------------------
        % Spatial and time period fixed effects + spatially lagged dependent variable
        info.lflag=0; % required for exact results
        info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
        % New routines to calculate effects estimates
        if info.model==0
            resultssa=sar_panel_FE(y,[xconstant x],W,T,info);
            select=logical (CODE);
            vnames=strvcat(namee(select,:));
            FREM(1:NV+2,2,ii)=num2cell(resultssa.parm(1:end-1));
            FREM(1:NV+2,3,ii)=num2cell(norm_prb(resultssa.tstat(1:end)));
        else
            resultssa=sar_panel_FE(y,x,W,T,info);
            CODE1=CODE;
            CODE1(1,2)=0;
            select=logical (CODE1);
            vnames=strvcat(namee(select,:));
            FREM(2:NV+2,2,ii)=num2cell(resultssa.parm(1:end-1));
            FREM(2:NV+2,3,ii)=num2cell(norm_prb(resultssa.tstat(1:end)));
        end
        
        % Print out coefficient estimates
        if print==1
            prt_sp(resultssa,vnames,1);
        end
        % Print out effects estimates
        spat_model=0;
        EFFECTS1=panel_effects_sar(resultssa,vnames,W);
        
        % needed for Hausman test later on
        logliklagsa=resultssa.lik;
        blagfesa=resultssa.parm(1:end-1);
        covblagfesa=resultssa.cov(1:end-1,1:end-1);
        % ----------------------------------------------------------------------------------------
        % Spatial and time period fixed effects + spatially spatial error model
        info.lflag=0; % required for exact results
        info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
        % New routines to calculate effects estimates
        if info.model==0
            resultsse=sem_panel_FE(y,[xconstant x],W,T,info);
            select=logical (CODE);
            vnames=strvcat(namee(select,:));
            FREM(1:NV+2,4,ii)=num2cell(resultsse.parm(1:end-1));
            FREM(1:NV+2,5,ii)=num2cell(norm_prb(resultsse.tstat(1:end)));
        else
            resultsse=sem_panel_FE(y,x,W,T,info);
            CODE2=CODE;
            CODE2(1,2)=0;
            select=logical (CODE2);
            vnames=strvcat(namee(select,:));
            FREM(2:NV+2,4,ii)=num2cell(resultsse.parm(1:end-1));
            FREM(2:NV+2,5,ii)=num2cell(norm_prb(resultsse.tstat(1:end)));
        end
        % Print out coefficient estimates
        if print==1
            prt_sp(resultsse,vnames,1);
            
        end
        % Print out effects estimates
        % needed for Hausman test later on
        logliklagse=resultsse.lik;
        blagfese=resultsse.parm(1:end-1);
        covblagfese=resultsse.cov(1:end-1,1:end-1);
        % ----------------------------------------------------------------------------------------
        % Spatial and time period fixed effects + spatially lagged dependent variable + spatially
        % independent variables
        % No bias correction
        info.lflag=0; % required for exact results
        info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
        % New routines to calculate effects estimates
        if info.model==0
            resultsd=sar_panel_FE(y,[xconstant x wx],W,T,info);
            CODE3=CODE;
            CODE3(1,1)=0;
            CODE3(1,2)=0;
            select=logical (CODE3);
            vnames=strvcat('Y','constant',namee(select,:),namee(select,:));
            FREM(1:2*NV+2,6,ii)=num2cell(resultsd.parm(1:end-1));
            FREM(1:2*NV+2,7,ii)=num2cell(norm_prb(resultsd.tstat(1:end)));
            FREM(2*NV+3,6,ii)=num2cell(resultsd.rsqr);
            FREM(2*NV+4,7,ii)=num2cell(resultsd.lik);
        else
            resultsd=sar_panel_FE(y,[x wx],W,T,info);
            CODE4=CODE;
            CODE4(1,1)=0;
            CODE4(1,2)=0;
            select=logical (CODE4);
            WM=repmat('W*',sum(CODE4)-numexo,1);
            namefi=namee(select,:);
            vnames=strvcat('Y',namee(select,:),[WM namefi(1:(NV-numexo),:)]);
            FREM(2:2*NV+2-numexo,6,ii)=num2cell(resultsd.parm(1:end-1));
            FREM(2:2*NV+2-numexo,7,ii)=num2cell(norm_prb(resultsd.tstat(1:end)));
            FREM(2*NV+3-numexo,6,ii)=num2cell(resultsd.rsqr);
            FREM(2*NV+4-numexo,7,ii)=num2cell(resultsd.lik);
        end
        
        % Print out coefficient estimates
        if print==1
            
            prt_sp(resultsd,vnames,1);
        end
        % Print out effects estimates
        spat_model=1;
        EFFECTS2=panel_effects_sdm(resultsd,vnames,W);
        
        % needed for Hausman test later on
        logliklag=resultsd.lik;
        blagfe=resultsd.parm(1:end-1);
        covblagfe=resultsd.cov(1:end-1,1:end-1);
        
        % Wald test for spatial Durbin model against spatial lag model
        fprintf(1,'Wald test for spatial Durbin model against spatial lag model');
        
        btemp=resultsd.parm;
        varcov=resultsd.cov;
        Rafg=zeros(K-numexo,2*K+2-numexo);
        for k=1:K-numexo
            Rafg(k,K+k)=1; % R(1,3)=0 and R(2,4)=0;
        end
        Wald_spatial1=(Rafg*btemp)'*inv(Rafg*varcov*Rafg')*Rafg*btemp
        prob_spatial1=1-chis_cdf (Wald_spatial1, K-numexo) % probability greater than 0.05 points to insignificance
        % LR test spatial Durbin model against spatial lag model (requires
        % estimation results of the spatial lag model to be available)
        
        % Wald test spatial Durbin model against spatial error model
        fprintf(1,'Wald test for spatial Durbin model against spatial error model');
        
        R=zeros(K-numexo,1);
        for k=1:K-numexo
            R(k)=btemp(2*K+1-numexo)*btemp(k)+btemp(K+k); % k changed in 1, 7/12/2010
            %   R(1)=btemp(5)*btemp(1)+btemp(3);
            %   R(2)=btemp(5)*btemp(2)+btemp(4);
        end
        Rafg=zeros(K-numexo,2*K+2-numexo);
        for k=1:K-numexo
            Rafg(k,k)    =btemp(2*K+1-numexo); % k changed in 1, 7/12/2010
            Rafg(k,K+k)  =1;
            Rafg(k,2*K+1-numexo)=btemp(k);
            %   Rafg(1,1)=btemp(5);Rafg(1,3)=1;Rafg(1,5)=btemp(1);
            %   Rafg(2,2)=btemp(5);Rafg(2,4)=1;Rafg(2,5)=btemp(2);
        end
        Wald_spatial2=R'*inv(Rafg*varcov*Rafg')*R
        prob_spatial2=1-chis_cdf (Wald_spatial2,K-numexo) % probability greater than 0.05 points to insignificance
        % LR test spatial Durbin model against spatial error model (requires
        % estimation results of the spatial error model to be available
        fprintf(1,'LR test spatial Durbin model against spatial lag model');
        LR_spatial1=-2*(resultssa.lik-resultsd.lik)
        prob_spatialLR1=1-chis_cdf (LR_spatial1,K) % probability greater than 0.05 points to insignificance
        
        fprintf(1,'LR test spatial Durbin model against spatial error model');
        
        LR_spatial2=-2*(resultsse.lik-resultsd.lik)
        prob_spatialLR2=1-chis_cdf (LR_spatial2,K) % probability greater than 0.05 points to insignificance
        
        ZA(3,2,ii)=num2cell(Wald_spatial1);
        ZA(4,2,ii)=num2cell(prob_spatial1);
        ZA(5,2,ii)=num2cell(Wald_spatial2);
        ZA(6,2,ii)=num2cell(prob_spatial2);
        ZA(7,2,ii)=num2cell(LR_spatial1);
        ZA(8,2,ii)=num2cell(prob_spatialLR1);
        ZA(9,2,ii)=num2cell(LR_spatial2);
        ZA(10,2,ii)=num2cell(prob_spatialLR2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% random effects estimator by ML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Spatial random effects and time period fixed effects + spatially lagged dependent variable
        % independent variables
        if  info.model==3
            [ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,2); % 2=time dummies
            resultssar=sar_panel_RE(ywith,xwith,W,T,info);
        else
            resultssar=sar_panel_RE(y,x,W,T,info);
        end
        CODE5=CODE;
        CODE5(1,2)=0;
        select=logical (CODE5);
        vnames=strvcat(namee(select,:));
        if print==1
            
            prt_sp(resultssar,vnames,1);
        end
        FREM(2:NV+3,8,ii)=num2cell(resultssar.parm(1:end-1));
        FREM(2:NV+3,9,ii)=num2cell(norm_prb(resultssar.tstat(1:end)));
        if info.model>0
            % ----------------------------------------------------------------------------------------
            % needed for Hausman test later on
            logliklagresa=resultssar.lik;
            blagresa=resultssar.parm(1:end-2);
            covblagresa=resultssar.cov(1:end-2,1:end-2);
            % Hausman test FE versus RE
            hausman=(blagfesa-blagresa)'*inv(covblagresa-covblagfesa)*(blagfesa-blagresa);
            dof=length(blagfesa);
            probability=1-chis_prb(abs(hausman),dof);
            % Note: probability < 0.025 implies rejection of random effects model in favor of fixed effects model
            % Use 0.025, since it is a one-sided test
            fprintf(1,'Hausman test-statistic, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',abs(hausman),dof,probability);
            ZA{1,1,ii}='Hausman test-statistic';
            ZA(1,2,ii)=num2cell(abs(hausman));
            ZA(1,3,ii)=num2cell(probability);
            % ----------------------------------------------------------------------------------------
        end
        % Print out effects estimates
        spat_model=0;
        EFFECTS3=panel_effects_sar(resultssar,vnames,W);
        
        % Spatial random effects and time period fixed effects + spatially spatial error model
        % independent variables
        if  info.model==3
            [ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,2); % 2=time dummies
            resultsser=sem_panel_RE(ywith,xwith,W,T,info);
        else
            resultsser=sem_panel_RE(y,x,W,T,info);
        end
        if print==1
            prt_sp(resultsser,vnames,1);
        end
        FREM(2:NV+3,10,ii)=num2cell(resultsser.parm(1:end-1));
        FREM(2:NV+3,11,ii)=num2cell(norm_prb(resultsser.tstat(1:end)));
        % Print out effects estimates
        
        % Spatial random effects and time period fixed effects + spatially lagged dependent variable + spatially
        % independent variables
        if  info.model==3
            [ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x wx],N,T,2); % 2=time dummies
            resultsdr=sar_panel_RE(ywith,xwith,W,T,info);
        else
            resultsdr=sar_panel_RE(y,[x wx],W,T,info);
        end
        CODE6=CODE;
        CODE6(1,1)=0;
        CODE6(1,2)=0;
        select=logical (CODE6);
        WM=repmat('W*',sum(CODE4)-numexo,1);
        namefi=namee(select,:);
        vnames=strvcat('Y',namee(select,:),[WM namefi(1:(NV-numexo),:)]);
        if print==1
            prt_sp(resultsdr,vnames,1);
        end
        FREM(2:2*NV+3-numexo,12,ii)=num2cell(resultsdr.parm(1:end-1));
        FREM(2:2*NV+3-numexo,13,ii)=num2cell(norm_prb(resultsdr.tstat(1:end)));
        if info.model>0
            
            % ----------------------------------------------------------------------------------------
            % needed for Hausman test later on
            logliklagre=resultsdr.lik;
            blagre=resultsdr.parm(1:end-2);
            covblagre=resultsdr.cov(1:end-2,1:end-2);
            % Hausman test FE versus RE
            hausman=(blagfe-blagre)'*inv(covblagre-covblagfe)*(blagfe-blagre);
            dof=length(blagfe);
            probability=1-chis_prb(abs(hausman),dof);
            % Note: probability < 0.025 implies rejection of random effects model in favor of fixed effects model
            % Use 0.025, since it is a one-sided test
            fprintf(1,'Hausman test-statistic, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',abs(hausman),dof,probability);
            ZA{2,1,ii}='Hausman test-statistic';
            ZA(2,2,ii)=num2cell(abs(hausman));
            ZA(2,3,ii)=num2cell(probability);
            % ----------------------------------------------------------------------------------------
        end
        % Print out effects estimates
        spat_model=1;
        EFFECTS4=panel_effects_sdm(resultsdr,vnames,W);
        % Wald test for spatial lag model
        fprintf(1,'Wald test for spatial Durbin model against spatial lag model');
        
        btemp=resultsdr.parm(1:2*K+2-numexo);
        varcov=resultsdr.cov(1:2*K+2-numexo,1:2*K+2-numexo);
        Rafg=zeros(K-numexo,2*K+2-numexo);
        for k=1:K-numexo
            Rafg(k,K+k)=1; % R(1,3)=0 and R(2,4)=0;
        end
        Wald_spatial3=(Rafg*btemp)'*inv(Rafg*varcov*Rafg')*Rafg*btemp
        prob_spatial3= 1-chis_cdf (Wald_spatial3, K-numexo) % probability greater than 0.05 points to insignificance
        
        % Wald test for spatial error model
        fprintf(1,'Wald test for spatial Durbin model against spatial error model');
        
        R=zeros(K-numexo,1);
        for k=1:K-numexo
            R(k)=btemp(2*K+1-numexo)*btemp(k)+btemp(K+k); % k changed in 1, 7/12/2010
            %   R(1)=btemp(5)*btemp(1)+btemp(3);
            %   R(2)=btemp(5)*btemp(2)+btemp(4);
        end
        Rafg=zeros(K-numexo,2*K+2-numexo);
        for k=1:K-numexo
            Rafg(k,k)    =btemp(2*K+1-numexo); % k changed in 1, 7/12/2010
            Rafg(k,K+k)  =1;
            Rafg(k,2*K+1-numexo)=btemp(k);
            %   Rafg(1,1)=btemp(5);Rafg(1,3)=1;Rafg(1,5)=btemp(1);
            %   Rafg(2,2)=btemp(5);Rafg(2,4)=1;Rafg(2,5)=btemp(2);
        end
        Wald_spatial4=R'*inv(Rafg*varcov*Rafg')*R
        prob_spatial4= 1-chis_cdf (Wald_spatial4,K-numexo) % probability greater than 0.05 points to insignificance
        % LR test spatial Durbin model against spatial lag model (requires
        % estimation results of the spatial lag model to be available)
        fprintf(1,'LR test spatial Durbin model against spatial lag model');
        LR_spatial3=-2*(resultssar.lik-resultsdr.lik)
        prob_spatialLR3=1-chis_cdf (LR_spatial3,K) % probability greater than 0.05 points to insignificance
        
        % LR test spatial Durbin model against spatial error model (requires
        % estimation results of the spatial error model to be available
        fprintf(1,'LR test spatial Durbin model against spatial error model');
        LR_spatial4=-2*(resultsser.lik-resultsdr.lik)
        
        ZA{3,1,ii}='Wald test for spatial Durbin model against spatial lag model';
        ZA{5,1,ii}='Wald test for spatial Durbin model against spatial error model';
        ZA{7,1,ii}='LR test spatial Durbin model against spatial lag model';
        ZA{9,1,ii}='LR test spatial Durbin model against spatial error model';
        ZA(3,3,ii)=num2cell(Wald_spatial3);
        ZA(4,3,ii)=num2cell(prob_spatial3);
        ZA(5,3,ii)=num2cell(Wald_spatial4);
        ZA(6,3,ii)=num2cell(prob_spatial4);
        ZA(7,3,ii)=num2cell(LR_spatial3);
        ZA(8,3,ii)=num2cell(prob_spatialLR3);
        ZA(9,3,ii)=num2cell(LR_spatial4);
        ZA{11,2,ii}='fixed effects estimator';
        ZA{11,3,ii}='random effects estimator';
        
        
        
        if CHOSE3==1
            EFFECTS2mo = reshape(EFFECTS1(1:3*size(x,2),[1 3]), [size(x,2), 6]);
        elseif CHOSE3==3
            EFFECTS2mo = reshape(EFFECTS2(1:3*size(x,2),[1 3]), [size(x,2), 6]);
        elseif CHOSE3==4
            EFFECTS2mo = reshape(EFFECTS3(1:3*size(x,2),[1 3]), [size(x,2), 6]);
        elseif CHOSE3==6
            EFFECTS2mo = reshape(EFFECTS4(1:3*size(x,2),[1 3]), [size(x,2), 6]);
        end
        EFFECTS2mo = EFFECTS2mo(:, [3,6,1,4,2,5]);
        EFFTO(1:size(x,2),2:7,ii)=num2cell(EFFECTS2mo);
        for ik=1:size(x,2)
            EFFTO{ik,1,ii}=strvcat(namefi(ik,:));
        end
        
        
    end
    
    for i=1:kk
        FREMM(:,1,i)=FREM(:,1,i);
        FREMM(:,2:3,i)=FREM(:,2*(CHOSE3-1)+2:2*(CHOSE3-1)+3,i);
    end
    
    for i=1:kk
        for ii=1:2*(2*maxv+3)
            C = mod(ii,2);
            if C==1
                FREMMF (ii,2,i)=FREMM(fix(ii/2)+1,2,i);
                FREMMF (ii,1,i)=FREMM(fix(ii/2)+1,1,i);
                
            else
                FREMMF (ii,2,i)=FREMM(ii/2,3,i);
                FREMMF {ii,1,i}=[];
            end
        end
    end
    for i=1:size(CHOSE, 2)
        ZAA(:,:,i)=ZA(:,:,CHOSE(i));
    end
    
    FREM1=[];
    FREM2=[];
    ZA1=[];
    EFFTO1=[];
    
    for n=1:kk
        FREM1 = cat(2,FREM1,FREMM(:,:,n));
    end
    
    for n=1:kk
        FREM2 = cat(2,FREM2,FREMMF(:,:,n));
    end
    
    for n=1:size(CHOSE, 2)
        ZA1 = cat(2,ZA1,ZAA(:,:,n));
    end
    for n=1:size(CHOSE2, 2)
        EFFTO1 = cat(1,EFFTO1,EFFTO(:,:,n));
    end
    
    
    xlswrite('Coefficient',FREM2,mohsen);
    xlswrite('Tests',ZA1,mohsen);
    xlswrite('Direct-Indirect Coefficient',EFFTO1,mohsen);
end
beep