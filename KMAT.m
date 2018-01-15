function [zModel,ii] = KMAT(fdinput,pfdata,Kc,K,Cs,rShoreFor,nsigma,keepall)
    % Profile model written by Kristen Splinter (October 2014) and implemented
    % in this code by Joshua Simmons 14/01/2015
    % Input:
    % fdinput = fielddata input - structure containing H,T,WL,Ho,Hb,diss,datesWLWave
    % pfdata = profile data as a structure .x,.z
    % Kc = parameter
    % K = parameter
    % Cs = parmaeter
    % keepall - boolean - 1 = keep z at all timesteps and 0 = only final z result kept.
    %
    % model description: 
    % first want to add in methods from Larson, Kubota and Erikson (04)
    % second lets add in PH12 model
    % OLD - define equilibrium beach slope based on Sunamura (1984) as Beq = 0.12*Hb./(g^0.5*T)*D^(0.25)
    % define runup based on Stockdon et al. (2006)
    % z is elevation above MSL, positive up
    % x positive cross-shore
    % h positive down (depth)
    % dh/dx is beach slope - need limits to define this by? time varying? +/-
    % swash?
    
    %constants - shouldnt be changed so can stayin in function
    g= 9.81;
    thetam = 30; %angle of repose, degrees
    mu = 0.7; %sediment packing factor

    %get all the values out f the structure for neatness
    H = fdinput.H; T = fdinput.T;
    WL = fdinput.WL; Ho = fdinput.Ho;
    Hb = fdinput.Hb; diss = fdinput.diss;
    datesWLWave = fdinput.datesWLWave;
    x1 = pfdata.x;
        
    %as a first stab, let's calculate beta +/- swash around WL+setup 
    setup = 0.17*Ho;    %from Davidson paper,  Ref Guza Thorton, 1981
    swash = 0.7*Ho;     %from Davidson paper, Ref Guza Thorton, 1982
    
    MWL = WL+setup;
    
    Lo = 1.56.*T.^2;
    diss(diss<0)=rShoreFor.*diss(diss<0);
    
    zModel = nan(length(MWL),length(x1));
    zModel(1,:) = pfdata.z;
    dQdx = zeros(size(zModel));
    zDuneToe=nan(size(zModel,1),1);
    dVResidual=zeros(size(zDuneToe));
    
    dt = abs(datesWLWave(2)-datesWLWave(1)).*24*3600;
    dx = abs(x1(2)-x1(1)); % XX this doesnt account for irregular grid
    
    %start the model
    ii=1;
    countDry=0;
    while ii<length(datesWLWave)
       %CALC sediment transport over 'active zone' between R2 and WL using
       %Stockdon formulatin for etabar, sigma_s, R2
        beta(ii) = abs(calcBeachSlopeSigmaSwash(x1,zModel(ii,:),min(WL(ii)-swash(ii),0),max(WL(ii)+swash(ii),0)));
        etabar(ii) = abs(0.35.*beta(ii).*sqrt(Ho(ii).*Lo(ii)));  %from Stockdon 06
        sigma_s(ii) = sqrt(Ho(ii).*Lo(ii).*(0.563.*(beta(ii).^2)+0.004))./2;%from Stockdon 06

        R2(ii) = calcR2(abs(beta(ii)),Ho(ii),T(ii)); %Stockdon eqn, includes setup but not tide (WL)
        magQ(ii,:) = Kc*2*sqrt(2*g).*(R2(ii)).^(3./2).*(1-(zModel(ii,:)-(WL(ii)))./(R2(ii))).^2;    %LKE04 use start at SWL

        con = find(zModel(ii,:)<=(WL(ii)),1,'first');
        magQ(ii,con:end)=magQ(ii,con);
        keep=find(zModel(ii,:)>(R2(ii)+WL(ii)),1,'last');

        magQ(ii,1:keep)=0;
        dissQ(ii)=-1.*diss(ii);
        dhdx=diff(zModel(ii,:))./dx;    %first gradient in slope
        avTerm(ii,:)=repmat(tand(thetam),size(dhdx))./(repmat(tand(thetam),size(dhdx)).^2-dhdx.^2); %avalanching
        avTerm(ii,1:keep)=0;    %still want to not allow this to turn on if the dune erosion model is on.%KS commented this out 24/03/15


        p(ii,:) =1-cdf('norm',zModel(ii,:),WL(ii)+etabar(ii),sigma_s(ii));
        Nc(ii,:) = p(ii,:).*(dt./T(ii)); 
        Qnet(ii,:)=magQ(ii,:).*dissQ(ii); %over single swash cycle, only where water is

        %START OF DUNE EROSION CHECK
        stScarp=find(zModel(ii,:)>=WL(ii)+R2(ii),1,'last');
        en=min(length((zModel(ii,:)))-2,find(zModel(ii,:)<=WL(ii),1,'first'));
        en2=min(length((zModel(ii,:)))-2,find(zModel(ii,:)<=WL(ii)-swash(ii),1,'first'));
        stDune=find(zModel(ii,:)>max(MWL+swash),1,'last');

        if isempty(stScarp)
            stScarp=1;
        end
        if isempty(en)
            en=length(zModel(ii,:))-2;
        end
        if isempty(en2)
            en2=length(zModel(ii,:))-2;
        end
        %dhdx=diff(zModel(ii,:))./dx;    %first gradient in slope
        dh2dx2=diff(zModel(ii,:),2)./(dx.^2);   %second gradient in slope
        [maxdh2,ind]=max(dh2dx2(stScarp:en2));   %max indicates an up turn to define a scarp toe
        [maxdh2Dune,indDune]=find(dh2dx2(stDune:en2)>0.2,1,'last');   %max indicates an up turn to define a scarp toe

        [~,ind2]=min(dh2dx2(stScarp:stScarp+ind));   %min if negative indicates a down turn and the top
        %    of the scarp hopefully withing the swash are
        [~,ind2Dune]=min(dh2dx2(stDune:stDune+indDune));   %min if negative indicates a down turn and the top
        if (dh2dx2(stDune+indDune)<maxdh2Dune & dh2dx2(stDune+indDune-2)<maxdh2Dune) %we've found a curvature that defines our dune toe (this is the one we will act on to turn on erosion)
            zDuneToe(ii)=zModel(ii,stDune+indDune);
            xDuneToe(ii) = x1(stDune+indDune);
            avTerm(ii,1:stDune+indDune)=0;
        else
            zDuneToe(ii)=NaN;
            xDuneToe(ii)=NaN;
        end

        zDuneCrest(ii)=max(zModel(ii,:));
        xDuneCrest(ii)=x1(zModel(ii,:)==max(zModel(ii,:)));
        
        if (maxdh2>0.3 & dh2dx2(stScarp+ind)<maxdh2 & dh2dx2(stScarp+ind-2)<maxdh2) %we've found a curvature that defines our scarp toe (this is the one we will act on to turn on erosion)
            avTerm(ii,1:stScarp+ind)=0;
            %%%%let's turn the dune erosion model on     
            zScarpToe(ii)=zModel(ii,stScarp+ind);
            xScarpToe(ii) = x1(stScarp+ind);
            if xScarpToe(ii)==xDuneToe(ii)      %this is if the scarp and the dune are the same thi
                zScarpCrest(ii)=zDuneCrest(ii);
                xScarpCrest(ii)=xDuneCrest(ii);
            else
                zScarpCrest(ii)=zModel(ii,stScarp+ind2);
                xScarpCrest(ii)=x1(stScarp+ind2);
            end             
            pOverwashScarp(ii) =1-cdf('norm',zScarpCrest(ii),WL(ii)+etabar(ii),sigma_s(ii));
            if pOverwashScarp(ii)<=0.50;   %if there is greater than a 50% chance that the scarp is overwashed, don't treat it as a dune erosion
                Bo(ii) = abs((zScarpToe(ii)-zModel(ii,end))./(xScarpToe(ii)-x1(end)));   %base slope between dune toe and WL which will be the projection inland of the dune retreat
                Bt=-Bo(ii);  %trajectory slope
                zbT = [Bt.*(x1(1:stScarp+ind)) + zScarpToe(ii)-Bt.*x1(stScarp+ind)];  %trajectory that dune toe receeds.
                V(ii) = sum(abs(diff(x1(1:2))).*(zModel(ii,1:stScarp+ind)));   %volume from back to dune toe
                clear Vc
                Vc = cumsum(abs(diff(x1(1:2))).*(zModel(ii,1:stScarp+ind)-zbT(1:stScarp+ind)));  %cumulative volume above the dune trajectory
                Vc = Vc(end)-Vc;  %as we are going from Shore -> sea in this model need to reverse everything
                etabarScarp(ii) = 0.35.*Bo(ii).*sqrt(Ho(ii).*Lo(ii));
                sigma_sScarp(ii) = sqrt(Ho(ii).*Lo(ii).*(0.563.*(Bo(ii).^2)+0.004))./2.*nsigma./2;
                zR(ii) = 1.1.*(etabarScarp(ii)+ sigma_sScarp(ii));
                zTotal(ii) = zR(ii) + WL(ii);
                pScarp(ii) =1-cdf('norm',zScarpToe(ii),etabarScarp(ii)+WL(ii),sigma_sScarp(ii));
                NcScarp(ii) = pScarp(ii).*(dt./T(ii));
                if ii>1
                    dV(ii) = 4.*Cs.*(max(zTotal(ii)-zScarpToe(ii),0)).^2.*NcScarp(ii);
                    dVT(ii) = dV(ii) - dVResidual(ii-1);
                else
                    dVT(ii) = 4.*Cs.*(max(zTotal(ii)-zScarpToe(ii),0)).^2.*NcScarp(ii);
                end
                %these cause a bit of a problem when the dune is quite large compared to
                %the runup... We'd still expect undercutting and failure

                if dVT(ii)<0
                    kk=stScarp+ind;
                else
                    [val kk] = (min((Vc-dVT(ii)).^2)); %find grid point where dune toe is
                end
                dxScarp(ii) = x1(kk)-xScarpToe(ii);
                dVResidual(ii) = Vc(kk)-dVT(ii);
                zModel(ii,:) = [zModel(ii,1:kk-1) zbT(kk:stScarp+ind) zModel(ii,stScarp+ind+1:end)];
            end    %end turn on dune model
        end
        % setup for 4thorder accuracy [tmp(1) tmp(1) Qnet(ii,:) tmp(end) tmp(end)];
        Qnet(ii,:)=Qnet(ii,:).*[0 avTerm(ii,:)];  
        tmp=Qnet(ii,:);
        clear tmp2
        tmp2=[Qnet(ii,:) tmp(end) tmp(end) 0 0; ...
             tmp(1) Qnet(ii,:) tmp(end) tmp(end) 0; ...
              0 tmp(1) tmp(1) Qnet(ii,:) tmp(end); ...
             0 0 tmp(1) tmp(1) Qnet(ii,:)];
        dQdxtmp = nansum([-1.*tmp2(1,:); 8.*tmp2(2,:); -8.*tmp2(3,:); tmp2(4,:)])./(12*dx);
        dQdx(ii,:)=dQdxtmp(3:end-2);
        filt = hanning(5);    %comment out KS 2015/03/24
        tmp=conv([dQdx(ii,1) dQdx(ii,1) dQdx(ii,:) dQdx(ii,end) dQdx(ii,end)],filt./sum(filt),'same');
        dQdx(ii,:)=tmp(3:end-2); clear tmp

        if ii==1
            zModel(ii+1,:)=zModel(ii,:)- 1./mu.*(dQdx(ii,:)).*Nc(ii,:);
        else
            zModel(ii+1,:)=zModel(ii,:)- 1./mu.*(3*dQdx(ii,:)-dQdx(ii-1,:))./2.*Nc(ii,:);
        end
        filt=hann(5); filt=filt./sum(filt);
        zModel(ii+1,stScarp+ind+2:end-2)=conv(zModel(ii+1,stScarp+ind:end),filt,'valid');

        %%calc 2nd derivative in h for diffusion term for stability
        tmp=zModel(ii+1,:);
        clear tmp2
        tm1=tmp(1)-(tmp(2)-tmp(1));
        tm2= tmp(1)-(tmp(3)-tmp(1));
        tp2=tmp(end)+(tmp(end)-tmp(end-2));
        tp1=tmp(end)+(tmp(end)-tmp(end-1));
        tmp2=[tmp tp1 tp2 0 0; ...
             tm1 tmp tp1 tp2 0; ...
             tm2 tm1 tmp tp1 tp2;...
              0 tm2 tm1 tmp tp1; ...
             0 0 tm2 tm1 tmp];

        dh2dx2tmp = nansum([-1.*tmp2(1,:); 16.*tmp2(2,:);-30.*tmp2(3,:); +16.*tmp2(4,:); -1.*tmp2(5,:)])./(12*dx.^2);
        dh2dx2Model(ii,:)=dh2dx2tmp(3:end-2);

        filt = hann(11);  %KS comment out 2015/03/24
        tmp=conv([repmat(dh2dx2Model(ii,1),1,5) dh2dx2Model(ii,:) repmat(dh2dx2Model(ii,end),1,5)],filt./sum(filt),'same');
        dh2dx2Model(ii,:)=tmp(6:end-5); clear tmp

        zModel(ii+1,:)=zModel(ii+1,:)+K.*dh2dx2Model(ii,:).*Nc(ii,:); %[0 diff(zModel(ii,:),1)].* assumes probability of impact

        filt=hann(5); filt=filt./sum(filt);
        zModel(ii+1,stScarp+ind+2:end-2)=conv(zModel(ii+1,stScarp+ind:end),filt,'valid');
        %this is the avalanching part
        pDuneToeWet(ii) =1-cdf('norm',zDuneToe(ii),WL(ii)+etabar(ii),sigma_s(ii));
        if pDuneToeWet(ii)<0.5
            countDry=countDry+1;
        else
            countDry=0;
        end
        if countDry>24; % & abs(zDuneCrest(ii)-zDuneToe(ii))>2 & abs((zDuneCrest(ii)-zDuneToe(ii))./(xDuneCrest(ii)-xDuneToe(ii)))>tand(thetam);
            dhdxtmp=diff(zModel(ii+1,:))./dx;  
            slumpId=find(abs(dhdxtmp)>tand(thetam));
            if(isempty(slumpId)==0)
                dhSlump=zModel(ii+1, slumpId(1))-zModel(ii+1,slumpId(length(slumpId))+1);
                if dhSlump>2
                    dxSlump=dhSlump./tand(thetam);
                    [~, p1x] = min((x1 - (xDuneToe(ii)-dxSlump/2)).^2);
                    [~, p2x] = min((x1 - (xDuneToe(ii)+dxSlump/2)).^2);
                    dhSlump2=zModel(ii+1,p1x)-zModel(ii+1,p2x);
                    dxsNew=x1(p2x)-x1(p1x);

                    while (dhSlump2./dxsNew)>tand(thetam)
                        p1x=p1x-1;
                        p2x=p2x+1;
                        dhSlump2=zModel(ii+1,p1x)-zModel(ii+1,p2x);
                         dxsNew=x1(p2x)-x1(p1x);
                    end
                    zSlump=-dhSlump2*(x1(p1x:p2x)-x1(p2x))./(x1(p2x)-x1(p1x));
                    zModel(ii+1,p1x:p2x)=zSlump+zModel(ii+1,p2x);
                    countDry=0;
                end
            end
        end

        ii=ii+1;
    end
        if ~keepall
            tmpz = zModel;
            clear zModel
            zModel(1,:) = tmpz(end,:);
        end
end
