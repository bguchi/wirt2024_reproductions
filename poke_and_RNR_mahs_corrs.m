% for a=1:12
%     for b=1:3
%         iFRhld=vertcat(alliFR_ports{b,goodsess(a)},alliFR_portsR{b,goodsess(a)},...
%             alliFR_portsNR{b,goodsess(a)});
%         [allcoe_RNR{b,a},allscorsRNR{b,a}]=pca(zscore(iFRhld));
%     end
%     iFRhld_real=vertcat(alliFR_ports_real{1,goodsess(a)}(1:73,:),...
%         alliFR_ports_real{2,goodsess(a)}(1:73,:),alliFR_ports_real{3,goodsess(a)}(1:73,:));
%     [allcoe_real{1,a},allscors_real{1,a}]=pca(zscore(iFRhld_real));
% end
% 
% 
% for b=1:3
%     x=0;
%     for a=1:12
%         for i=1:size(allcoe_RNR{b,a},1)
%             x=x+1;
%             coehld_all{b}(x,1:3)=allcoe_RNR{b,a}(i,1:3);
%             coehld_all{b}(x,4:6)=allcoe_real{1,a}(i,1:3);
%         end
%     end
% end
% 
% 
% % calculate single cell MahD RNR curves for each port

for a=1:12
    for b=1:3
        clear iFRhld coe sco idx cells
        for i=1:size(alliFR_ports_real{b,goodsess(a)},2)
            iFRhld=vertcat(alliFR_ports_real{b,goodsess(a)},alliFR_portsR{b,goodsess(a)},...
                alliFR_portsNR{b,goodsess(a)});
            [coe,sco]=pca(zscore(iFRhld));
            % sco=tsne(zscore(iFRhld));
            numtri=(length(sco)/2);numR=length(alliFR_portsR{b,goodsess(a)}); numNR=size(alliFR_portsNR{b,goodsess(a)},1);
            for c=1:(length(sco)/2)-1
                idx{1}=c:c+2;
                idx{2}=numtri+1:numtri+10;
                idx{3}=numtri+numR+1: numtri+numR+10;
                res=MahDis_James_accel(zscore(iFRhld(:,i)),idx,.05);
                cellMahs_RNR{b,a}{i}(c,1)=res.Mah(1,2);
                cellMahs_RNR{b,a}{i}(c,2)=res.Mah(1,3);
            end
        end
    end
end

for b=1:3
    x=0;
    for a=1:12
        for i=1:size(alliFR_ports_real{b,goodsess(a)},2)
            x=x+1;
            hld=(cellMahs_RNR{b,a}{i}(:,1))-(cellMahs_RNR{b,a}{i}(:,2));
            cellMahs_RNRall{b}(:,x)=hld(1:70,1);
            mod=ones(70,1);
            mod(:,2)=TrialTimes{1,goodsess(a)}(3:72,1)';
            [betasRNR_time_real,~,~,~,statslin]=regress(movmean(hld(1:70,1),mv),mod);
            rho_cells_RNR_time_real(x,b)=statslin(1,1);
            mod=ones(70,1);
            
            mod(:,2)=allmodelmovs{1,goodsess(a)}(3:72,4)';
            [betasRNR_beh_real,~,~,~,statslin]=regress(movmean(hld(1:70,1),mv),mod);
            rho_cells_RNR_beh_real(x,b)=statslin(1,1);
            
        end
    end
end

rho_cells_RNR_time_real(isnan(rho_cells_RNR_time_real))=0;
rho_cells_RNR_beh_real(isnan(rho_cells_RNR_beh_real))=0;


for a=1:548
rho_cells_RNR_time_real(a,4)=max(rho_cells_RNR_time_real(a,1:3));
rho_cells_RNR_beh_real(a,4)=max(rho_cells_RNR_beh_real(a,1:3));
end
figure; plot(histc((rho_cells_RNR_beh_real(:,4)-rho_cells_beh_real(:,1)),[-1:0.1:1])); hold on;
plot(histc((rho_cells_RNR_time_real(:,4)-rho_cells_time_real(:,1)),[-1:0.1:1]));

% based on difference curves for 25_75 - 75_25
cellMahs_RNR_lohiDiff=(cellMahs_RNRall{1})-(cellMahs_RNRall{3});
 x=0;
    for a=1:12
        for i=1:size(alliFR_ports_real{b,goodsess(a)},2)
            x=x+1;
            hld=cellMahs_RNR_lohiDiff(:,x);
            cellMahs_RNRall{b}(:,x)=hld(1:70,1);
            mod=ones(70,1);
            mod(:,2)=TrialTimes{1,goodsess(a)}(3:72,1)';
            [betasRNR_time_real,~,~,~,statslin]=regress(movmean(hld(1:70,1),1),mod);
            rho_cells_RNRdiff_time_real(x,1)=statslin(1,1);
            mod=ones(70,1);
            
            mod(:,2)=allmodelmovs{1,goodsess(a)}(3:72,4)';
            [betasRNR_beh_real,~,~,~,statslin]=regress(movmean(hld(1:70,1),1),mod);
            rho_cells_RNRdiff_beh_real(x,1)=statslin(1,1);
            
        end
        idx=find(~isnan(Mah_sme_align_diff(:,a)));
        hld=Mah_sme_align_diff(idx(3:72),a);

            mod=ones(70,1);
            mod(:,2)=TrialTimes{1,goodsess(a)}(3:72,1)';
            [betasRNR_time_real,~,~,~,statslin]=regress(movmean(hld(1:70,1),1),mod);
            rho_ens_RNRdiff_time_real(a,1)=statslin(1,1);
            mod=ones(70,1);
            
            mod(:,2)=allmodelmovs{1,goodsess(a)}(3:72,4)';
            [betasRNR_beh_real,~,~,~,statslin]=regress(movmean(hld(1:70,1),1),mod);
            rho_ens_RNRdiff_beh_real(a,1)=statslin(1,1);
    end
    
    
 for b=1:3
    for a=1:12   
        for c=1:3
        hld=allscorsRNR{b,a}(3:72,c);
                mod=ones(70,1);
            mod(:,2)=allmodelmovs{1,goodsess(a)}(3:72,4)';
            [betasRNR_beh_real,~,~,~,statslin]=regress(movmean(hld(1:70,1),1),mod);
            rho_ens_RNRPCs_beh_real{b}(a,c)=statslin(1,1);
        end
    end
 end
 
 
 % regress time off of ens dMah RNR -25-75 diff curves
 
ens_Mah_RNRdiff_real_resids=NaN(200,12);
for s=1:12
    idx=find(~isnan(Mah_sme_align_diff(:,s)));
    hldr=Mah_sme_align_diff(idx(3:72),s);
    hldr(:,2)=TrialTimes{goodsess(s)}(3:length(hldr)+2,1);
    [~,~,ens_Mah_RNRdiff_real_resids(1:70,s),~,stats]=regress(hldr(:,1), [ ones(length(hldr),1) hldr(:,2)]);
end
     
    for a=1:12   
        hld=ens_Mah_RNRdiff_real_resids(3:72,a);
                mod=ones(70,1);
            mod(:,2)=allmodelmovs{1,goodsess(a)}(3:72,4)';
            [betasRNR_beh_resids,~,~,~,statslin]=regress(movmean(hld(1:70,1),1),mod);
            rho_ens_RNRPCs_beh_resids(a,1)=statslin(1,1);
        end
