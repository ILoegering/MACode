%% Plot moment arm data for all trials
data_path = '';
plot_path = '';
subjectID = '';
kinemID = '';
createCmbdMAvAngPlot();

%% Fit second-order polynomial to good trials
fits_path = '';
cond = '';
trials = '';
db_path = '';
db_name = '';
getMAFit();