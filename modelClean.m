%model clean
function [SkillRMSE,SkillR2,zModelFinal]=modelClean(MWL,x1,xp,zp,WL, swash,Ho,Lo,T,thetam,zpI,stdz,dateProfile,datesWLWave,Kc,Cs,K,diss,dx,dt,mu,nsigma)
    g=9.81
    zModel = nan(length(MWL),length(x1));
    zModel(1,:)=interp1(xp,zp(1,:),x1);
    for ss=1:size(zp,1)
        zpI(ss,:)=interp1(xp,zp(ss,:),x1);
        zstdI(ss,:)=interp1(xp,stdz(ss,:),x1);
    end
    zTrue = nan(length(MWL),length(x1));
    zTrueStd=zTrue;
    [cc,ia,ib]=intersect(dateProfile,datesWLWave);
    zTrue(ib,:)=zpI(ia,:); %is the measured profile I think
    zTrueStd(ib,:)=zstdI(ia,:);
    dQdx = zeros(size(zModel));
    clear beta R2 magQ dissQ Qnet tmp dQdxtmp dQdx p Nc
    zDuneToe=nan(size(zModel,1),1);
    dVResidual=zeros(size(zDuneToe));

    ii=1;
    countDry=0;
    while ii<ib(end) %:length(MWL)
        %CALC sediment transport over 'active zone' between R2 and WL
        beta(ii) = abs(calcBeachSlopeSigmaSwash(x1,zModel(ii,:),min(WL(ii)-swash(ii),0),max(WL(ii)+swash(ii),0)));
        etabar(ii) = abs(0.35.*beta(ii).*sqrt(Ho(ii).*Lo(ii)));  %from Stockdon 06
        sigma_s(ii) = sqrt(Ho(ii).*Lo(ii).*(0.563.*(beta(ii).^2)+0.004))./2;%from Stockdon 06

        R2(ii) = calcR2(abs(beta(ii)),Ho(ii),T(ii)); %Stockdon eqn, includes setup but not tide (WL)
        magQ(ii,:) = Kc*2*sqrt(2*g).*(R2(ii)).^(3./2).*(1-(zModel(ii,:)-(WL(ii)))./(R2(ii))).^2;    %LKE04 use start at SWL


        con = find(zModel(ii,:)<=(WL(ii)),1,'first');
        magQ(ii,con:end)=magQ(ii,con);
        keep=find(zModel(ii,:)>(R2(ii)+WL(ii)),1,'last');

        magQ(ii,1:keep)=0;
        %dissQ(ii) = tand(thetam)./(tand(thetam).^2-beta(ii).^2).*(beta(ii)-betaEq(ii)); %+ values mean offshore and should be flattening
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

        [mindh2,ind2]=min(dh2dx2(stScarp:stScarp+ind));   %min if negative indicates a down turn and the top
    %    of the scarp hopefully withing the swash are
        [mindh2Dune,ind2Dune]=min(dh2dx2(stDune:stDune+indDune));   %min if negative indicates a down turn and the top
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
%        
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
        %Bt=-betaEq(ii);
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

        % dz/dt=-1/mu*dQx/dx.% x pos offshore, Qpos offshore, z (or h) pos up
        %add in effect of gradients in Q
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
    end % %%%
%     %     
    tmp2=zModel(:,find(x1==xpStart):find(x1==xpEnd));
    tmp=zTrue(:,find(x1==xpStart):find(x1==xpEnd));
%    SkillRMSE(jj,mm,nn)=rmse(tmp(~isnan(tmp)),tmp2(~isnan(tmp)));
%    SkillR2(jj,mm,nn)=(corr(tmp(~isnan(tmp)),tmp2(~isnan(tmp)))).^2;
   SkillRMSE=rmse(tmp(~isnan(tmp)),tmp2(~isnan(tmp)));
   SkillR2=(corr(tmp(~isnan(tmp)),tmp2(~isnan(tmp)))).^2;
   zModelFinal=zModel(ii-1,:);
