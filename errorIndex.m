function CalculatedIndex = errorIndex(ObservedMatrix, PredictedMatrix, Index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description: This function calculates the prediction error based  on  the                    %
%% observed data, predicted values, and favorite error index.                                   %
%%                                                                                              %
%% Interface:                                                                                   %
%%                                                                                              %
%% Inputs:                                                                                      %
%%                                                                                              %
%%                                                                                              %
%% ObservedMatrix:  N * m , N is the number of observations, m is the prediction horizon        %
%% PredictedMatrix: N * m , N is the number of predictions , m is the prediction horizon        %
%% Index: string,  this input indicates the index type. The possible choices for this input is: %
%%                                                                                              %
%% MSE     ==> Mean Squared Error                         = mean(e^2_t)                         %
%% RMSE    ==> Root Mean Squared Error                    = sqrt(MSE)                           %
%% MAE     ==> Mean Absolute Error                        = mean(abs(e_t))                      %
%% MDAE    ==> Median Absolute Error                      = median(abs(e_t))                    %
%% MAPE    ==> Mean Absolute Percentage Error             = mean((p_t))                         %
%% MDAPE   ==> Median Absolute Percentage Error           = median((p_t))                       %
%% SMAPE   ==> Symmetric Mean Absolute Percentage Error   = mean(2abs(Y_t - F_t)/(Y_t + F_t))   %
%% SMDAPE  ==> Symmetric Median Absolute Percentage Error = median(2abs(Y_t - F_t)/(Y_t + F_t)) %
%% MRAE    ==> Mean Relative Absolute Error               = mean(abs(r_t))                      %
%% MDRAE   ==> Median Relative Absolute Error             = median(abs(r_t))                    %
%% GMRAE   ==> Geometric Mean Relative Absolute Error     = gmean(abs(r_t))                     %
%% RELMAE  ==> Relative Mean Absolute Error               = MAE / MAE_b                         %
%% RELMSE  ==> Relative Mean Squared Error                = MSE / MSE_b                         %
%% RELRMSE ==> Relative Root Mean Squared Error           = RMSE / RMSE_b                       %
%% LMR     ==> Log Mean Squared Error Ratio               = log(RelRMSE)                        %
%% PB      ==> Percentage Better                          = 100 * mean(I{r_t < 1})              %
%% PBMAE   ==> Percentage Better Mean Absolute Error      = 100 * mean(I{MAE < MAE_b})          %
%% PBMSE   ==> Percentage Better Mean Squared Error       = 100 * mean(I{MSE < MSE_b})          %
%% PBRMSE  ==> Percentage Better Root Mean Squared Error  = 100 * sqrt(mean(I{MSE < MSE_b}))    %
%%                                                                                              %
%% Output:                                                                                      %
%% CalculatedIndex: m * 1, vector doubles, m is the prediction horizon                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 0)
    disp('Dude, I need at least two inputs to calucualte the error index!');
    CalculatedIndex = 0;
    Index = 'MSE';
elseif (nargin == 1)
    disp('Dude, I need at least two inputs to calucualte the error index!');
    CalculatedIndex = 0;
    Index = 'MSE';
elseif (nargin == 2)
    Index = 'MSE';
end

Index = upper(Index);

BaseError = ObservedMatrix(2:end) - ObservedMatrix(1:end - 1);

ErrorB = sum(abs(BaseError) , 1);
MAEB = ErrorB / size(BaseError , 1);

ErrB = BaseError' * BaseError;
ErrorB = diag (ErrB);
clear ErrB;

MSEB  = ErrorB / size(BaseError , 1);
RMSEB = sqrt(ErrorB / size(BaseError , 1));

switch Index
    
    case 'MSE'
        Err = (ObservedMatrix - PredictedMatrix).^2;
        Error = sum(Err , 1);
        clear Err;
        CalculatedIndex = Error / size(ObservedMatrix , 1);
        
    case 'RMSE'
        Err = (ObservedMatrix - PredictedMatrix).^2;
        Error = sum(Err , 1);
        clear Err;
        CalculatedIndex = sqrt(Error / size(ObservedMatrix , 1));
        
    case 'MAE'
        Err = abs(ObservedMatrix - PredictedMatrix);
        Error = sum(Err , 1);
        clear Err;
        CalculatedIndex = Error / size(ObservedMatrix , 1);
        
    case 'MDAE'
        Err = abs(ObservedMatrix - PredictedMatrix);
        Error = median(Err , 1);
        clear Err;
        CalculatedIndex = Error;
        
    case 'MAPE'
        PercentErr = 100 * (ObservedMatrix - PredictedMatrix) ./ ObservedMatrix;
        Error = sum(abs(PercentErr) , 1);
        clear PercentErr;
        CalculatedIndex = Error / size(ObservedMatrix , 1);
        
    case 'MDAPE'
        PercentErr = abs(100 * (ObservedMatrix - PredictedMatrix) ./ ObservedMatrix);
        Error = median(PercentErr , 1);
        clear PercentErr;
        CalculatedIndex = Error;
        
    case 'SMAPE'
        SymmetricErr = 2 * abs(ObservedMatrix - PredictedMatrix) ./ (ObservedMatrix + PredictedMatrix);
        Error = sum(SymmetricErr , 1);
        clear SymmetricErr;
        CalculatedIndex = Error / size(ObservedMatrix , 1);
        
    case 'SMDAPE'
        SymmetricErr = 2 * abs(ObservedMatrix - PredictedMatrix) ./ (ObservedMatrix + PredictedMatrix);
        Error = median(SymmetricErr , 1);
        clear SymmetricErr;
        CalculatedIndex = Error;
        
    case 'MRAE'
        RelErr = ObservedMatrix(2 : end) - PredictedMatrix(2 : end);
        Error = sum(abs(RelErr./BaseError),1);
        clear RelErr;
        CalculatedIndex = Error / (size(ObservedMatrix , 1) - 1);
        
    case 'MDRAE'
        RelErr = ObservedMatrix(2 : end) - PredictedMatrix(2 : end);
        Error = median(abs(RelErr./BaseError) , 1);
        clear RelErr;
        CalculatedIndex = Error;
        
    case 'GMRAE'
        RelErr = ObservedMatrix(2 : end) - PredictedMatrix(2 : end);
        Error  = abs(RelErr./BaseError);
        clear RelErr;
        CalculatedIndex = (prod(Error,1)) ^ (1/(size(ObservedMatrix , 1) - 1));
        
    case 'RELMAE'
        Err = abs(ObservedMatrix - PredictedMatrix);
        Error = sum(Err , 1);
        clear Err;
        MAE = Error / size(ObservedMatrix , 1);
        CalculatedIndex = MAE ./ MAEB;
        
    case 'RELRMSE'
        Err = (ObservedMatrix - PredictedMatrix).^2;
        Error = sum(Err , 1);
        clear Err;
        RMSE = sqrt(Error / size(ObservedMatrix , 1));
        CalculatedIndex = RMSE ./ RMSEB;
        
    case 'RELMSE'
        Err = (ObservedMatrix - PredictedMatrix).^2;
        Error = sum(Err , 1);
        clear Err;
        MSE = Error / size(ObservedMatrix , 1);
        CalculatedIndex = MSE ./ MSEB;
    
    case 'LMR'
        Err = (ObservedMatrix - PredictedMatrix)' * (ObservedMatrix - PredictedMatrix);
        Error = diag (Err);
        clear Err;
        MSE = Error / size(ObservedMatrix , 1);
        CalculatedIndex = log(MSE / MSEB);
    
    case 'PB'
        RelErr = ObservedMatrix(2 : end) - PredictedMatrix(2 : end);
        Error = abs(RelErr./BaseError);
        clear RelErr;
        CalculatedIndex = 100 * sum((Error < 1) , 1) / (size(ObservedMatrix , 1) - 1);
        
    case 'PBMAE'
        Err = abs(ObservedMatrix(2:end) - PredictedMatrix(2:end));
        Error = (Err < abs(BaseError));
        clear Err;
        CalculatedIndex = 100 * sum((Error < 1) , 1) / size(BaseError , 1);
    
    case 'PBMSE'
        Err  = (ObservedMatrix(2:end) - PredictedMatrix(2:end)).^2;
        ErrB = (BaseError).^2;
        Error =  (Err < ErrB);
        clear Err;
        CalculatedIndex = sum(Error , 1) / size(ObservedMatrix , 1);
        
    case 'PBRMSE'
        Err  = (ObservedMatrix(2:end) - PredictedMatrix(2:end)).^2;
        ErrB = (BaseError).^2;
        Error =  (Err < ErrB);
        clear Err;
        CalculatedIndex = sqrt(sum(Error , 1) / size(ObservedMatrix , 1));
        
    otherwise
        warning('Unknown Error Index!');
        CalculatedIndex = 0;
end


end