clear all;
close all;
clc ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxv=12;
kk=14;
T=50; % number of time periods
N=648; % number of regions
A=xlsread('cigarette.xls');
    
FREMMMK=cell(4,9,8);
CHOSE=[1:kk];
pr=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FREMM=cell(maxv+9,9,kk);
FREMMpr=cell(maxv+3,3,kk);
FREMMFpr=cell(2*(maxv+1)+2,2,kk);
FREMMM=cell(4,9,kk);


for ii=1:kk
    % dimensions of the problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W1=xlsread('matrix.xls',ii);
W=zeros(N,N);
for p=1:size(W1,1)
    aa=W1(p,1);
    bb=W1(p,2);
    W(aa,bb)=1/W1(p,3);
end

W=normw(W); % function of LeSage
    ms=size(A,2);
    AA=A(:,6:ms);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %namee = strvcat('tourist','constant','gdpp','pric','ope','govsta','SociCon','InvPro','IntCon','ExtCon','Corr','MilPol','RelTen','LawOr','EthTen','DemAcc','BurQua');
    y = A(:,5);

    %namee =    strvcat('eco','constant','CO2','CO2-2','CO2-3','CO2-31','GDP','GDP2','UR','ENP','ENG','TR','POP','Built-up Land','Carbon','Cropland','Fishing Grounds','Forest Products','Grazing Land','co2g','WT','WT2','LOD','Sunspot','Thermosteric','Barystatic','GMSL');
    namee =    strvcat('tem','constant','Urban','Cropland','Forest','Forest (Unknown/Other)','Forest (Evergreen, needle leaf) ','Forest ( Evergreen, broad leaf) ','Forest (Deciduous, needle leaf) ','Forest (Deciduous, broad leaf) ','Forest (Mixed) ','CO2E','CO2F','PRE','WINUn+','WINUn-','WINVn+','WINVn-','WINU+','WINU-','WINV+','WINV-','HCL','MCL','LCL','ML-CL','TCL','Thermosteric','Barystatic','Sunspot','LOD','GMSL','CHG','GHC-FOSSIL','GHC-FOSSIL_AGRI','CHG2-Fossil','CHG2-LULUCF','CHG2-TOT');
    CODE=zeros(1,size(namee,1));
    CODE(1, 1)=1;
    CODE(1, 2)=1;
    CODE(1, 3)=1;
    CODE(1, 4)=1;
    CODE(1, 5)=1;
    CODE(1, 6)=0;
    CODE(1, 7)=0;
    CODE(1, 8)=0;
    CODE(1, 9)=0;
    CODE(1, 10)=0;
    CODE(1, 11)=0;
    CODE(1, 12)=1;
    CODE(1, 13)=0;
    CODE(1, 14)=1;
    CODE(1, 15)=1;
    CODE(1, 16)=1;
    CODE(1, 17)=1;
    CODE(1, 18)=1;
    CODE(1, 19)=0;
    CODE(1, 20)=0;
    CODE(1, 21)=1;
    CODE(1, 22)=1;
    CODE(1, 23)=1;
    CODE(1, 24)=0;
    CODE(1, 25)=0;
    CODE(1, 25)=1;
    CODE(1, 30)=0;
    CODE(1, 31)=1;
    CODE(1, 32)=0;
    CODE(1, 36)=0;
        CODE(1, 37)=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%must change%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    select=logical (CODE);
    vnames=namee(select,:);

    NV=sum(CODE)-2;
    MSEL=CODE(1, 3:size(CODE,2));
    selectx=logical (MSEL);
    x=AA(:,selectx);
    xconstant=ones(N*T,1);
    [nobs K]=size(x);
    % written by: J.Paul Elhorst summer 2010
    % University of Groningen
    % Department of Economics
    % 9700AV Groningen
    % the Netherlands
    % j.p.elhorst@rug.nl
    %
    % REFERENCES:
    % Elhorst JP (2010) Matlab Software for Spatial Panels. Under review.
    %
    % Elhorst JP (2010) Spatial Panel Data Models. In Fischer MM, Getis A (Eds.)
    % Handbook of Applied Spatial Analysis, Ch. C.2. Springer: Berlin Heidelberg New York.
    %
    % ----------------------------------------------------------------------------------------
    FREMM{1,2,ii}='Pooled OLS';
    FREMM{1,4,ii}='Spatial fixed effects';
    FREMM{1,6,ii}='Time-period fixed effects';
    FREMM{1,8,ii}='Spatial and time-period fixed effects';
    
    % ols estimation
    results=ols(y,[xconstant x]);
    prt_reg(results,vnames,1);
    sige=results.sige*((nobs-K)/nobs);
    loglikols=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
    
    for i=1:NV+1
        FREMM{i+1,1,ii}=vnames(i+1,:);
    end
    FREMM(2:NV+2,2,ii)=num2cell(results.beta(1:end));
    FREMM(2:NV+2,3,ii)=num2cell(norm_prb(results.tstat(1:end)));
    % The (robust)LM tests developed by Elhorst
    FREMM{NV+7,1,ii}='LogL';
    FREMM(NV+7,2,ii)=num2cell(loglikols);
    FREMM{NV+8,1,ii}='R^2';
    FREMM(NV+8,2,ii)=num2cell(results.rsqr);
    LMsarsem_panel(results,W,y,[xconstant x]); % (Robust) LM tests
    
    % The lm tests developed by Donald Lacombe
    % see http://www.rri.wvu.edu/lacombe/~lacombe.htm
    
    lm1=lmlag_panel(y,[xconstant x],W);
    prt_tests(lm1);
    FREMM(NV+3,2,ii)=num2cell(lm1.lm);
    FREMM(NV+3,3,ii)=num2cell(lm1.prob);
    FREMM{NV+3,1,ii}='LM spatial lag';
    
    lm2=lmerror_panel(y,[xconstant x],W);
    prt_tests(lm2);
    FREMM(NV+4,2,ii)=num2cell(lm2.lm);
    FREMM(NV+4,3,ii)=num2cell(lm2.prob);
    FREMM{NV+4,1,ii}='LM spatial error';
    
    lm3=lmlag_robust_panel(y,[xconstant x],W);
    prt_tests(lm3);
    FREMM(NV+5,2,ii)=num2cell(lm3.lm);
    FREMM(NV+5,3,ii)=num2cell(lm3.prob);
    FREMM{NV+5,1,ii}='robust LM spatial lag';
    
    lm4=lmerror_robust_panel(y,[xconstant x],W);
    prt_tests(lm4);
    FREMM(NV+6,2,ii)=num2cell(lm4.lm);
    FREMM(NV+6,3,ii)=num2cell(lm4.prob);
    FREMM{NV+6,1,ii}='robust LM spatial error';
    
    % ----------------------------------------------------------------------------------------
    % spatial fixed effects + (robust) LM tests for spatial lag and spatial error model
    
    % fixed effects, within estimator
    % demeaning of the y and x variables
    model=1;
    [ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
    results=ols(ywith,xwith);
    FREMM(3:NV+2,4,ii)=num2cell(results.beta(1:end));
    FREMM(3:NV+2,5,ii)=num2cell(norm_prb(results.tstat(1:end)));
    CODE1=CODE;
    CODE1(1,2)=0;
    select=logical (CODE1);
    vnames=namee(select,:);
    
    prt_reg(results,vnames);
    sfe=meanny-meannx*results.beta; % including the constant term
    yme = y - mean(y);
    et=ones(T,1);
    error=y-kron(et,sfe)-x*results.beta;
    rsqr1 = error'*error;
    rsqr2 = yme'*yme;
    FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
    sige=results.sige*((nobs-K)/nobs);
    logliksfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
    fprintf(1,'spatial fixed effects + (robust) LM tests for spatial lag and spatial error model');
    FREMM(NV+7,4,ii)=num2cell(logliksfe);
    FREMM(NV+8,4,ii)=num2cell(FE_rsqr2);
    
    LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests
    
    lm1=lmlag_panel(ywith,xwith,W);
    prt_tests(lm1);
    FREMM(NV+3,4,ii)=num2cell(lm1.lm);
    FREMM(NV+3,5,ii)=num2cell(lm1.prob);
    
    lm2=lmerror_panel(ywith,xwith,W);
    prt_tests(lm2);
    FREMM(NV+4,4,ii)=num2cell(lm2.lm);
    FREMM(NV+4,5,ii)=num2cell(lm2.prob);
    
    lm3=lmlag_robust_panel(ywith,xwith,W);
    prt_tests(lm3);
    FREMM(NV+5,4,ii)=num2cell(lm3.lm);
    FREMM(NV+5,5,ii)=num2cell(lm3.prob);
    
    lm4=lmerror_robust_panel(ywith,xwith,W);
    prt_tests(lm4);
    FREMM(NV+6,4,ii)=num2cell(lm4.lm);
    FREMM(NV+6,5,ii)=num2cell(lm4.prob);
    if pr==1
        FREMMpr(1:NV+1, 1:3, ii) =FREMM(2:NV+2, [1 4 5], ii);
        FREMMpr(NV+2:NV+3, 1:3, ii) =FREMM([NV+7 NV+8], [1 4 5], ii);
    end
    % ----------------------------------------------------------------------------------------
    % time-period fixed effects + (robust) LM tests for spatial lag and spatial error model
    
    % fixed effects, within estimator
    % demeaning of the y and x variables
    model=2;
    [ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
    results=ols(ywith,xwith);
    FREMM(3:NV+2,6,ii)=num2cell(results.beta(1:end));
    FREMM(3:NV+2,7,ii)=num2cell(norm_prb(results.tstat(1:end)));
    CODE1=CODE;
    CODE1(1,2)=0;
    select=logical (CODE1);
    vnames=namee(select,:);
    prt_reg(results,vnames);
    tfe=meanty-meantx*results.beta; % including the constant term
    yme = y - mean(y);
    en=ones(N,1);
    error=y-kron(tfe,en)-x*results.beta;
    rsqr1 = error'*error;
    rsqr2 = yme'*yme;
    FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
    sige=results.sige*((nobs-K)/nobs);
    logliktfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
    fprintf(1,'time-period fixed effects + (robust) LM tests for spatial lag and spatial error model');
    FREMM(NV+7,6,ii)=num2cell(logliktfe);
    FREMM(NV+8,6,ii)=num2cell(FE_rsqr2);
    
    LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests
    
    lm1=lmlag_panel(ywith,xwith,W);
    prt_tests(lm1);
    FREMM(NV+3,6,ii)=num2cell(lm1.lm);
    FREMM(NV+3,7,ii)=num2cell(lm1.prob);
    
    lm2=lmerror_panel(ywith,xwith,W);
    prt_tests(lm2);
    FREMM(NV+4,6,ii)=num2cell(lm2.lm);
    FREMM(NV+4,7,ii)=num2cell(lm2.prob);
    
    lm3=lmlag_robust_panel(ywith,xwith,W);
    prt_tests(lm3);
    FREMM(NV+5,6,ii)=num2cell(lm3.lm);
    FREMM(NV+5,7,ii)=num2cell(lm3.prob);
    
    lm4=lmerror_robust_panel(ywith,xwith,W);
    prt_tests(lm4);
    FREMM(NV+6,6,ii)=num2cell(lm4.lm);
    FREMM(NV+6,7,ii)=num2cell(lm4.prob);
    if pr==2
        FREMMpr(1:NV+1, 1:3, ii) =FREMM(2:NV+2, [1 6 7], ii);
        FREMMpr(NV+2:NV+3, 1:3, ii) =FREMM([NV+7 NV+8], [1 6 7], ii);
    end
    % ----------------------------------------------------------------------------------------
    % spatial and time period fixed effects + (robust) LM tests for spatial lag and spatial error model
    
    % fixed effects, within estimator
    % demeaning of the y and x variables
    model=3;
    [ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
    results=ols(ywith,xwith);
    FREMM(3:NV+2,8,ii)=num2cell(results.beta(1:end));
    FREMM(3:NV+2,9,ii)=num2cell(norm_prb(results.tstat(1:end)));
    CODE1=CODE;
    CODE1(1,2)=0;
    select=logical (CODE1);
    vnames=namee(select,:);
    prt_reg(results,vnames);
    intercept=mean(y)-mean(x)*results.beta;
    sfe=meanny-meannx*results.beta-kron(en,intercept);
    tfe=meanty-meantx*results.beta-kron(et,intercept);
    yme = y - mean(y);
    ent=ones(N*T,1);
    error=y-kron(tfe,en)-kron(et,sfe)-x*results.beta-kron(ent,intercept);
    rsqr1 = error'*error;
    rsqr2 = yme'*yme;
    FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
    sige=results.sige*((nobs-K)/nobs);
    loglikstfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
    FREMM(NV+7,8,ii)=num2cell(loglikstfe);
    FREMM(NV+8,8,ii)=num2cell(FE_rsqr2);
    
    fprintf(1,'spatial and time period fixed effects + (robust) LM tests for spatial lag and spatial error model');
    
    LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests
    
    lm1=lmlag_panel(ywith,xwith,W);
    prt_tests(lm1);
    FREMM(NV+3,8,ii)=num2cell(lm1.lm);
    FREMM(NV+3,9,ii)=num2cell(lm1.prob);
    
    lm2=lmerror_panel(ywith,xwith,W);
    prt_tests(lm2);
    FREMM(NV+4,8,ii)=num2cell(lm2.lm);
    FREMM(NV+4,9,ii)=num2cell(lm2.prob);
    
    lm3=lmlag_robust_panel(ywith,xwith,W);
    prt_tests(lm3);
    FREMM(NV+5,8,ii)=num2cell(lm3.lm);
    FREMM(NV+5,9,ii)=num2cell(lm3.prob);
    
    lm4=lmerror_robust_panel(ywith,xwith,W);
    prt_tests(lm4);
    FREMM(NV+6,8,ii)=num2cell(lm4.lm);
    FREMM(NV+6,9,ii)=num2cell(lm4.prob);
    if pr==3
        FREMMpr(1:NV+1, 1:3, ii) =FREMM(2:NV+2, [1 8 9], ii);
        FREMMpr(NV+2:NV+3, 1:3, ii) =FREMM([NV+7 NV+8], [1 8 9], ii);
    end
    % ----------------------------------------------------------------------------------------
    % Tests for the joint significance of spatial and/or time-period fixed effects
    LR=-2*(logliktfe-loglikstfe);
    dof=N;
    probability=1-chis_prb(LR,dof);
    % Note: probability > 0.05 implies rejection of spatial fixed effects
    fprintf(1,'LR-test joint significance spatial fixed effects, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',LR,dof,probability);
    FREMM{NV+9,1,ii}='LR-test';
    FREMM(NV+9,6,ii)=num2cell(LR);
    FREMM(NV+9,7,ii)=num2cell(probability);
    
    LR=-2*(logliksfe-loglikstfe);
    dof=T;
    probability=1-chis_prb(LR,dof);
    % Note: probability > 0.05 implies rejection of spatial fixed effects
    fprintf(1,'LR-test joint significance time-periode fixed effects, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',LR,dof,probability);
    fprintf(1,'Note: probability > 0.05 implies rejection of spatial fixed effects');
    
    FREMM(NV+9,4,ii)=num2cell(LR);
    FREMM(NV+9,5,ii)=num2cell(probability);
    FREMMM(1:4,:,ii)=FREMM([1 NV+3 NV+4 NV+9],:,ii);
    
    for jj=1:2*(NV+1)
        C = mod(jj,2);
        if C==1
            FREMMFpr (jj,2,ii)=FREMMpr(fix(jj/2)+1,2,ii);
            FREMMFpr (jj,1,ii)=FREMMpr(fix(jj/2)+1,1,ii);
            
        else
            FREMMFpr (jj,2,ii)=FREMMpr(jj/2,3,ii);
            FREMMFpr {jj,1,ii}=[];
        end
    end
    FREMMFpr(2*NV+3:2*NV+4,1:2, ii)=FREMMpr(NV+2:NV+3,1:2, ii);
    
end
for i=1:size(CHOSE, 2)
    FREMMMK(:,:,i)=FREMMM(:,:,CHOSE(i));
end

FREMM1=[];
FREMMMK1=[];

for n=1:kk
    FREMM1 = cat(1,FREMM1,FREMM(:,:,n));
end

xlswrite('MOGHAYESEH',FREMM1);

for n=1:size(CHOSE, 2)
    FREMMMK1 = cat(1,FREMMMK1,FREMMMK(:,:,n));
end

xlswrite('KHOLASEMOGHAYESEH',FREMMMK1);

FREMMFpr1=[];
for n=1:kk
    FREMMFpr1 = cat(2,FREMMFpr1,FREMMFpr(:,:,n));
end
xlswrite('MOGHAYESEH1',FREMMFpr1);
