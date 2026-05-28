
function Fit=BenFunctions(X,FunIndex,Dim,I)
      switch FunIndex
          case 1
                Fit= M_CEM(I, X);
          otherwise
            error('Invalid FunIndex. For the proposed image segmentation method, use FunIndex = 1.');
    end

end
