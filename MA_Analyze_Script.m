%% Plot moment arm data for all trials
root = 'G:\My Drive\UW NMBL\Tendon\Tendon Mechanics Database\TMD_TDchildren\TMD Subject Data\';
subjectID = 'TMD_TD01';
data_path = [root '\' subjectID '\Data Analysis\MA\Data'];
plot_path = [root '\' subjectID '\Data Analysis\MA\Plots'];
kinemID = 'Ex';
createCmbdMAvAngPlot(data_path,plot_path,subjectID,kinemID);

%% Fit second-order polynomial to good trials
fits_path = [root '\' subjectID '\Data Analysis\MA\Fits'];
cond = 'K'; % A for Ankle, K for Knee
trials = [1 2];
db_path = root;
db_name = 'TMD';
getMAFit(data_path,fits_path,subjectID,cond,kinemID,trials,db_path,db_name);