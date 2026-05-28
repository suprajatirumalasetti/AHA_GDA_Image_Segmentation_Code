function [BestX,BestF,HisBestFit,VisitTable]=AHA_GDA(FunIndex,MaxIt,nPop,I, Dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FunIndex: The index of function.                      %
% MaxIt: The maximum number of iterations.              %
% nPop: The size of hummingbird population.             %
% Dim: The dimensionality of problem.                    %
% BestX: The best solution found so far.                %
% BestF: The best fitness corresponding to BestX.       %
% HisBestFit: History best fitness over iterations.     %
% VisitTable: The visit table.                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Low,Up] = FunRange(FunIndex);
    nPop4 = floor(nPop / 4);
    pop1 = LatinHypercubeInitialization(nPop4, Dim, Low, Up);
    pop2 = SobolInitialization(nPop4, Dim, Low, Up);
    pop3 = HaltonInitialization(nPop4, Dim, Low, Up);
    pop4 = SierpinskyInitialization(nPop4, Dim, Low, Up);
    PopPos = [pop1; pop2; pop3; pop4];
    PopPos = sort(round(PopPos), 2);
    remPop = nPop - 4*nPop4;
    if remPop > 0
        extraPop = rand(remPop, Dim) .* (Up - Low) + Low;
        extraPop = sort(round(extraPop));
        PopPos = [PopPos; extraPop];
    end
    PopFit = zeros(1, nPop);
    for i = 1:nPop
        PopFit(i) = BenFunctions(PopPos(i,:), FunIndex, Dim, I);
    end
    [BestF, idxBest] = min(PopFit);
    BestX = PopPos(idxBest,:);
    HisBestFit = zeros(MaxIt,1);
    VisitTable = zeros(nPop);
    VisitTable(logical(eye(nPop))) = NaN;
    W = BestF * 1.1; 
    decayRate = (W - BestF) / MaxIt;
    for It = 1:MaxIt
        DirectVector = zeros(nPop, Dim);
        for i = 1:nPop
            r = rand;
            if r < 1/3
                RandDim = randperm(Dim);
                if Dim >= 3
                    RandNum = ceil(rand*(Dim-2)+1);
                else
                    RandNum = ceil(rand*(Dim-1)+1);
                end
                DirectVector(i, RandDim(1:RandNum)) = 1;
            elseif r > 2/3
                DirectVector(i,:) = 1;
            else
                RandNum = ceil(rand*Dim);
                DirectVector(i,RandNum) = 1;
            end
            if rand < 0.5 
                [MaxUnvisitedTime, TargetFoodIndex] = max(VisitTable(i,:));
                MUT_Index = find(VisitTable(i,:) == MaxUnvisitedTime);
                if length(MUT_Index) > 1
                    [~, Ind] = min(PopFit(MUT_Index));
                    TargetFoodIndex = MUT_Index(Ind);
                end
                newPopPos = PopPos(TargetFoodIndex,:) + randn * DirectVector(i,:) .* ...
                    (PopPos(i,:) - PopPos(TargetFoodIndex,:));
                newPopPos = sort(round(SpaceBound(newPopPos, Up, Low)));
                newPopFit = BenFunctions(newPopPos, FunIndex, Dim, I);
                if (newPopFit < PopFit(i)) || (newPopFit <= W)
                    PopFit(i) = newPopFit;
                    PopPos(i,:) = newPopPos;
                    VisitTable(i,:) = VisitTable(i,:) + 1;
                    VisitTable(i, TargetFoodIndex) = 0;
                    VisitTable(:,i) = max(VisitTable,[],2) + 1;
                    VisitTable(i,i) = NaN;
                else
                    VisitTable(i,:) = VisitTable(i,:) + 1;
                    VisitTable(i, TargetFoodIndex) = 0;
                end

            else 
                newPopPos = PopPos(i,:) + randn * DirectVector(i,:) .* PopPos(i,:);
                newPopPos = sort(round(SpaceBound(newPopPos, Up, Low)));
                newPopFit = BenFunctions(newPopPos, FunIndex, Dim, I);
                if (newPopFit < PopFit(i)) || (newPopFit <= W)
                    PopFit(i) = newPopFit;
                    PopPos(i,:) = newPopPos;
                    VisitTable(i,:) = VisitTable(i,:) + 1;
                    VisitTable(:,i) = max(VisitTable,[],2) + 1;
                    VisitTable(i,i) = NaN;
                else
                    VisitTable(i,:) = VisitTable(i,:) + 1;
                end
            end
        end
        if mod(It, 2 * nPop) == 0
            [~, MigrationIndex] = max(PopFit);
            PopPos(MigrationIndex,:) = sort(round(rand(1,Dim).*(Up-Low) + Low));
            PopFit(MigrationIndex) = BenFunctions(PopPos(MigrationIndex,:), FunIndex, Dim, I);
            VisitTable(MigrationIndex,:) = VisitTable(MigrationIndex,:) + 1;
            VisitTable(:,MigrationIndex) = max(VisitTable,[],2) + 1;
            VisitTable(MigrationIndex, MigrationIndex) = NaN;
        end
        [currBestF, currBestIdx] = min(PopFit);
        if currBestF < BestF
            BestF = currBestF;
            BestX = PopPos(currBestIdx,:);
        end
        HisBestFit(It) = BestF;
        W = W - decayRate;
        if W < BestF
            W = BestF;
        end
    end
end


