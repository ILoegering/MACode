function [output,idx] = selectWhere(db,select,tFields,tValues,operators)
%selectWhere - extract records that fulfill specified conditions
%Extract the fields indicated in select of the records from database db
%for which the specified conditions are met. Conditions are defined by
%tFields, tValues, and operators. For example, if the first entries in
%tFields, tValues, and operators are 'age', 50, and '<', respectively, then
%only records for which the 'age' field is less than 50 will be extracted.
%This function was modeled after the SELECT statement with a WHERE clause
%in the database language SQL.
%
% Syntax:  [out,idx] = selectWhere(s,outFields,tFields,tValues,varargin)
%
% Inputs:
%    db (required) - struct array (1 x numRecords)
%           Database
%    select (required) - string cell array (1 x numOutFields)
%           Fields to extract from db
%    tFields (required) - cell array (1 x numConditions)
%           Fields that must take on values in tValues
%    tValues (required) - cell array (1 x numConditions)
%           Values on which fields in tFields must take or handles of
%           functions, such as @isnan, that test a condition
%    operators (required) - cell array (1 x numConditions)
%           Logical operators (e.g., '==','~=','>','<','>=','<=','~'); if
%           the value in tValues is non-numeric (e.g., char array or
%           function handle), use either an empty char array ('') or '=='
%           to test for equality or either negation ('~') or '~=' to
%           test for inequality
%
% Outputs:
%    out - struct array (1 x numRecordsSelected)
%           Subset of input database db 
%    idx - double array (1 x numRecordsSelected)
%           Indices of selected records
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Isaac Loegering
% UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: isaacloegering@gmail.com
% February 2019; Last revision: 14-May-2019
%------------- BEGIN CODE --------------
% Find indices of rows that satisfy the specified conditions
testArrays = nan(numel(tFields),numel(db));
ind = 1;
for i = 1:numel(tValues)
    try
        if (isnumeric(tValues{i}))
            eval(['testArrays(ind,:) = ([db.(tFields{i})] ' operators{i} ' tValues{i});'])
        elseif(ischar(tValues{i}))
            if (strcmp(operators{i},'=='))
                operators{i} = '';
            elseif (strcmp(operators{i},'~='))
                operators{i} = '~';
            end
            eval(['testArrays(ind,:) = ' operators{i} 'strcmp({db.(tFields{i})},tValues{i});'])
        elseif(isa(tValues{i},'function_handle'))
            if (strcmp(operators{i},'=='))
                operators{i} = '';
            elseif (strcmp(operators{i},'~='))
                operators{i} = '~';
            end
            eval(['testArrays(ind,:) = ' operators{i} 'cellfun(tValues{i},{db.(tFields{i})});'])
        else
            %TODO: This error may be caught by the try-catch clause,
            %resulting in two error messages being displayed. Are they both
            %applicable?
            error(['This function was not designed to handle the data '...
                'type of tValues{' num2str(i) '}.']);
        end
    catch
        error(['Condition ' num2str(i) ' is incorrectly defined. The '...
            'combination of arguments specified in tFields, tValues, '...
            'and operators is invalid. Ensure data types match and '...
            'logical operators specified are appropriate for those '...
            'comparisons of those data types.'])
    end
    ind = ind + 1;
end
testDecision = ones(1,size(testArrays,2));
for i = 1:size(testArrays,1)
    testDecision = testDecision & testArrays(i,:);
end
idx = find(testDecision==1);
output = db(idx);
fn = fieldnames(db);
fn = fn(~ismember(fn,select)==1);
output = rmfield(output,fn);
end
%------------- END OF CODE --------------