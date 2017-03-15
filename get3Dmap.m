function [model model2 model3 model4 model5]    = get3Dmap(es_train, es_test, variable, nGrid, flag_sp, bins1, bins2, bins3)
        X_train(:,1) = [eval(['es_train.' variable{1}])];
        X_train(:,2) = [eval(['es_train.' variable{2}])];
        X_train(:,3) = [eval(['es_train.' variable{3}])];
        X_test(:,1) = [eval(['es_test.' variable{1}])];
        X_test(:,2) = [eval(['es_test.' variable{2}])];
        X_test(:,3) = [eval(['es_test.' variable{3}])];
        
        Y_train = X_train(:,1) + nGrid(1)*(X_train(:,2)-1) + nGrid(1)*nGrid(2)*(X_train(:,3)-1);
        Y_test  = X_test(:,1) + nGrid(1)*(X_test(:,2)-1) + nGrid(1)*nGrid(2)*(X_test(:,3)-1);
        
        pos2d1 = nGrid(1)*(X_test(:,2)-1) + X_test(:,1);
        pos2d2 = nGrid(1)*(X_test(:,3)-1) + X_test(:,1);
        pos2d3 = nGrid(2)*(X_test(:,3)-1) + X_test(:,2);
        
        
        if nGrid(1)>100
            n1 = 10;
            grids1 = 1:2:n1+2;
        else
            n1 = 5;
            grids1 = 1:n1;
        end
        if nGrid(2)>100
            n2 = 10;
            grids2 = 1:2:n2+2;
        else
            n2 = 5;
            grids2 = 1:n2;
        end
        if nGrid(3)>100
            n3 = 10;
            grids3 = 1:2:n3+2;
        else
            n3 = 5;
            grids3 = 1:n3;
        end
        
        for icell = 1:numCells
            display(['Processing cell: ' num2str(icell)]);
            
            zua_train = conv(es_train.spikeTrain(:,icell), sfilt, 'same');
            zua_test = conv(es_test.spikeTrain(:,icell), sfilt, 'same');
            test = es_test.spikeTrain(:,icell);
            
            % 3D: PVR
            scMap = full(sparse(Y_train, 1, (zua_train), nGrid(1)*nGrid(2)*nGrid(3),1));
            scMap = reshape(scMap,nGrid(1),nGrid(2),nGrid(3));
            occMap = full(sparse(Y_train, 1, 1, nGrid(1)*nGrid(2)*nGrid(3),1));
            occMap = reshape(occMap,nGrid(1),nGrid(2),nGrid(3));
            % 1D marginals:
            % P
            scMap1d1 = full(sparse(X_train(:,1), 1, (zua_train), nGrid(1), 1));
            occMap1d1 = full(sparse(X_train(:,1), 1, 1, nGrid(1), 1));
            % V
            scMap1d2 = full(sparse(X_train(:,2), 1, (zua_train), nGrid(2), 1));
            occMap1d2 = full(sparse(X_train(:,2), 1, 1, nGrid(2), 1));
            % R
            scMap1d3 = full(sparse(X_train(:,3), 1, (zua_train), nGrid(3), 1));
            occMap1d3 = full(sparse(X_train(:,3), 1, 1, nGrid(3), 1));
            % 2D Marginals
            % PV
            scMap2d1 = full(sparse(X_train(:,1), X_train(:,2), (zua_train), nGrid(1), nGrid(2)));
            occMap2d1 = full(sparse(X_train(:,1), X_train(:,2), 1, nGrid(1), nGrid(2)));
            % PR
            scMap2d2 = full(sparse(X_train(:,1), X_train(:,3), (zua_train), nGrid(1), nGrid(3)));
            occMap2d2 = full(sparse(X_train(:,1), X_train(:,3), 1, nGrid(1), nGrid(3)));
            % VR
            scMap2d3 = full(sparse(X_train(:,2), X_train(:,3), (zua_train), nGrid(2), nGrid(3)));
            occMap2d3 = full(sparse(X_train(:,2), X_train(:,3), 1, nGrid(2), nGrid(3)));
            
            %get the FR map by smoothing
            %             grids1 = [1 30 60 90 120 150];
            %             grids2 = [1 15 30];
            %             grids3 = [1 15 30];
            
            EV  = zeros(length(grids1),length(grids2),length(grids3));
            EV2 = zeros(length(grids1),length(grids2),length(grids3));
            EV3 = zeros(length(grids1),length(grids2),length(grids3));
            EV4 = zeros(length(grids1),length(grids2),length(grids3));
            EV5 = zeros(length(grids1),length(grids2),length(grids3));
            
            for win1=1:length(grids1)
                fprintf('\n P %d ',win1);
                for win2=1:length(grids2)
                    fprintf(' \n ');
                    for win3=1:length(grids3)
                        fprintf('. ');
                        %                         display(['Cell ' num2str(icell) ' ' num2str(win1) ' ' num2str(win2) ' ' num2str(win3) ' ']);
                        %3D
                        FRMap = special_smooth_3d(scMap, [1./grids1(win1) 1./grids2(win2) 1./grids3(win3)], flag_sp, bins1, bins2, bins3, [n1 n2 n3])...
                            ./special_smooth_3d(occMap, [1./grids1(win1) 1./grids2(win2) 1./grids3(win3)], flag_sp, bins1, bins2,bins3, [n1 n2 n3]);
                        
                        fullFRvector    = reshape(FRMap, [],1);
                        pred = fullFRvector(Y_test);
                        spred  = conv(pred,sfilt, 'same');
                        EV(win1,win2,win3) = calculateError(zua_train, zua_test, spred, test, pred);
                        
                        %2D
                        FRMap2d1 = special_smooth_2d(scMap2d1, [1./grids1(win1) 1./grids2(win2)], [0 0], bins1, bins2, [n1 n2])...
                            ./special_smooth_2d(occMap2d1, [1./grids1(win1) 1./grids2(win2)], [0 0], bins1, bins2, [n1 n2]);
                        FRvector2d1    = reshape(FRMap2d1, [],1);
                        FRMap2d2 = special_smooth_2d(scMap2d2, [1./grids1(win1) 1./grids3(win3)], [0 0], bins1, bins3, [n1 n3])...
                            ./special_smooth_2d(occMap2d2, [1./grids1(win1) 1./grids3(win3)], [0 0], bins1, bins3, [n1 n3]);
                        FRvector2d2    = reshape(FRMap2d2, [],1);
                        FRMap2d3 = special_smooth_2d(scMap2d3, [1./grids2(win2) 1./grids2(win3)], [0 0], bins2, bins3, [n2 n3])...
                            ./special_smooth_2d(occMap2d3, [1./grids2(win2) 1./grids3(win3)], [0 0], bins2, bins3, [n2 n3]);
                        FRvector2d3    = reshape(FRMap2d3, [],1);
                        %1D
                        FRMap1d1 = special_smooth_1d(scMap1d1, 1./grids1(win1), flag_sp(1), bins1, n1)...
                            ./special_smooth_1d(occMap1d1, 1./grids1(win1), flag_sp(1), bins1, n1);
                        FRMap1d2 = special_smooth_1d(scMap1d2, 1./grids2(win2), flag_sp(2), bins2, n2)...
                            ./special_smooth_1d(occMap1d2, 1./grids2(win2), flag_sp(2), bins2, n2);
                        FRMap1d3 = special_smooth_1d(scMap1d2, 1./grids2(win3), flag_sp(3), bins3, n3)...
                            ./special_smooth_1d(occMap1d3, 1./grids2(win3), flag_sp(3), bins3, n3);
                        
                        %model2     P X V X R PmVmR
                        pred2   = FRMap1d1(X_test(:,1)).*FRMap1d2(X_test(:,2)).*FRMap1d3(X_test(:,3));
                        spred2   = conv(pred2,sfilt, 'same');
                        [EV2(win1,win2,win3) corr2(win1,win2,win3) L2(win1,win2,win3) Q2(win1,win2,win3)] = calculateError(zua_train, zua_test, spred2, test, pred2);
                        %model3     P X (V,R) PmVR
                        temp_pred = FRvector2d3(pos2d3);
                        pred3     = FRMap1d1(X_test(:,1)).*temp_pred;
                        spred3   = conv(pred3,sfilt, 'same');
                        [EV3(win1,win2,win3) corr3(win1,win2,win3) L3(win1,win2,win3) Q3(win1,win2,win3)] = calculateError(zua_train, zua_test, spred3, test, pred3);
                        %model4     (P,V) X R PVmR
                        temp_pred = FRvector2d1(pos2d1);
                        pred4     = FRMap1d3(X_test(:,3)).*temp_pred;
                        spred4   = conv(pred4,sfilt, 'same');
                        [EV4(win1,win2,win3) corr4(win1,win2,win3) L4(win1,win2,win3) Q4(win1,win2,win3)] = calculateError(zua_train, zua_test, spred4, test, pred4);
                        %model5     (P,R) X V PRmV
                        temp_pred = FRvector2d2(pos2d2);
                        pred5     = FRMap1d2(X_test(:,2)).*temp_pred;
                        spred5   = conv(pred5,sfilt, 'same');
                        [EV5(win1,win2,win3) corr5(win1,win2,win3) L5(win1,win2,win3) Q5(win1,win2,win3)] = calculateError(zua_train, zua_test, spred5, test, pred5);
                        
                    end
                end
            end
            EV  = round(100*EV);
            EV2 = round(100*EV2);
            EV3 = round(100*EV3);
            EV4 = round(100*EV4);
            EV5 = round(100*EV5);
            
            idx  = find(reshape(EV ,[],1) == max(reshape(EV ,[],1)),1,'first');
            idx2 = find(reshape(EV2,[],1) == max(reshape(EV2,[],1)),1,'first');
            idx3 = find(reshape(EV3,[],1) == max(reshape(EV3,[],1)),1,'first');
            idx4 = find(reshape(EV4,[],1) == max(reshape(EV4,[],1)),1,'first');
            idx5 = find(reshape(EV5,[],1) == max(reshape(EV5,[],1)),1,'first');
            
            
            %model      P,V,R     PVR
            model(icell).swin = idx;
            %model2     P X V X R PmVmR
            model2(icell).swin = idx2;
            if ~isempty(idx2)
                model2(icell).EV   = 0.01*EV2(idx2);
                model2(icell).corr   = corr2(idx2);
                model2(icell).L   = L2(idx2);
                model2(icell).Q   = Q2(idx2);
                model2(icell).bins{1} = bins1;
                model2(icell).bins{2} = bins2;
                model2(icell).bins{3} = bins3;
            else
                model2(icell).swin   = NaN;
                model2(icell).EV   = NaN;
                model2(icell).corr = NaN;
                model2(icell).L   = NaN;
                model2(icell).Q   = NaN;
                model2(icell).bins{1} = bins1;
                model2(icell).bins{2} = bins2;
                model2(icell).bins{3} = bins3;
            end
            %model3     P X (V,R) PmVR
            model3(icell).swin = idx3;
            if ~isempty(idx3)
                model3(icell).EV   = 0.01*EV3(idx3);
                model3(icell).corr   = corr3(idx3);
                model3(icell).L   = L3(idx3);
                model3(icell).Q   = Q3(idx3);
                model3(icell).bins{1} = bins1;
                model3(icell).bins{2} = bins2;
                model3(icell).bins{3} = bins3;
            else
                model3(icell).swin   = NaN;
                model3(icell).EV   = NaN;
                model3(icell).corr = NaN;
                model3(icell).L   = NaN;
                model3(icell).Q   = NaN;
                model3(icell).bins{1} = bins1;
                model3(icell).bins{2} = bins2;
                model3(icell).bins{3} = bins3;
            end%model4     (P,V) X R PVmR
            model4(icell).swin = idx4;
            if ~isempty(idx4)
                model4(icell).EV   = 0.01*EV4(idx4);
                model4(icell).corr   = corr4(idx4);
                model4(icell).L   = L4(idx4);
                model4(icell).Q   = Q4(idx4);
                model4(icell).bins{1} = bins1;
                model4(icell).bins{2} = bins2;
                model4(icell).bins{3} = bins3;
            else
                model4(icell).swin   = NaN;
                model4(icell).EV   = NaN;
                model4(icell).corr = NaN;
                model4(icell).L   = NaN;
                model4(icell).Q   = NaN;
                model4(icell).bins{1} = bins1;
                model4(icell).bins{2} = bins2;
                model4(icell).bins{3} = bins3;
            end
            %model5     (P,R) X V PRmV
            model5(icell).swin = idx5;
            if ~isempty(idx5)
                model5(icell).EV   = 0.01*EV5(idx5);
                model5(icell).corr = corr5(idx5);
                model5(icell).L   = L5(idx5);
                model5(icell).Q   = Q5(idx5);
                model5(icell).bins{1} = bins1;
                model5(icell).bins{2} = bins2;
                model5(icell).bins{3} = bins3;
            else
                model5(icell).swin   = NaN;
                model5(icell).EV   = NaN;
                model5(icell).corr = NaN;
                model5(icell).L   = NaN;
                model5(icell).Q   = NaN;
                model5(icell).bins{1} = bins1;
                model5(icell).bins{2} = bins2;
                model5(icell).bins{3} = bins3;
            end
            XY = rem(idx,length(grids1)*length(grids2));
            idxX = rem(XY, length(grids1));
            if idxX==0; idxX = length(grids1); end
            idxY = ((XY-idxX)./length(grids1)) + 1;
            if idxY==0; idxY = length(grids2); end
            idxZ = ((idx-idxX-(idxY-1)*length(grids1))./(length(grids1)*length(grids2))) + 1;
            if idxZ==0; idxZ = length(grids3); end
            
            win1 = (idxX);
            win2 = (idxY);
            win3 = (idxZ);
            model(icell).swin = [win1 win2 win3];
            
            
            
            display([num2str(idx) '_' num2str(idxX) '_' num2str(idxY) '_' num2str(idxZ)]);
            
            if ~(isempty(win1) & isnan(win1) & isempty(win2) & isnan(win2) & isempty(win3) & isnan(win3))
                t1 = tic;
                FRMap = special_smooth_3d(scMap, [1./grids1(win1) 1./grids2(win2) 1./grids3(win3)], flag_sp, bins1, bins2, bins3, [n1 n2 n3])...
                    ./special_smooth_3d(occMap, [1./grids1(win1) 1./grids2(win2) 1./grids3(win3)], flag_sp, bins1, bins2,bins3, [n1 n2 n3]);
                toc(t1)
                fullFRvector    = reshape(FRMap, [],1);
                pred = fullFRvector(Y_test);
                spred  = conv(pred,sfilt, 'same');
                
                model(icell).tuning = FRMap;
                [model(icell).EV model(icell).corr model(icell).L model(icell).Q model(icell).train_mean]= calculateError(zua_train, zua_test, spred, test, pred);
                model(icell).bins{1} = bins1;
                model(icell).bins{2} = bins2;
                model(icell).bins{3} = bins3;
            else
                model(icell).swin   = NaN;
                model(icell).EV = NaN;
                model(icell).corr = NaN;
                model(icell).L = NaN;
                model(icell).Q = NaN;
                model(icell).bins{1} = bins1;
                model(icell).bins{2} = bins2;
                model(icell).bins{3} = bins3;
            end
            %
        end
    end