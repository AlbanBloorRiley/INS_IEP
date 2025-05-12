function Sys_Out = Sys_Output(Vals,Ops,SysFixed)
%SYS_OUTPUT Forms optimal spin structure
%Sys_Out = Sys_Output(Y,Ops,SysFixed) calculates the spin structure defined
%by Vals.*Ops + SysFixed wher Vals is an array of doubles, Ops is a cell 
% structure containing a basis of all the spin structure fields,  each 
% corresponding to the equivalent element of Vals and SysFixed is the spin 
% structure formed of only fixed variables.
%
%If Vals is a structure array then Sys_Out will also be a structure array 
% of equal size where the function is applied element-wise.

if isstruct(SysFixed)
    S=SysFixed.S;
else
    S=SysFixed;
end



f = fieldnames(Ops);
CellOps = struct2cell(Ops);
k=1;
for i=1:length(CellOps)       %iterate through degree of Stev and then J/ee
    if isempty(CellOps{i})     %Skip empty
        %AssignOps{N}{i} = {};
        for  N= 1:size(Vals,2)
            OpsOut{N}{i} = [];
        end
        continue
    end
    for j=1:length(CellOps{i})     % iterate through elements of B2,B4,etc.
        for  N= 1:size(Vals,2)      %iterate through deflations
            AssignOps{N}{i}{j} = CellOps{i}{j}*Vals(k,N);
        end
        k=k+1;
    end
    if length(CellOps{i})>1
        for  N= 1:size(Vals,2)
            OpsOut{N}{i} = zeros(size(AssignOps{N}{i}{1}));
            for j = 1:length(CellOps{i})
                %                 for j = 2:length(CellOps{i})
                %                 OpsOut{N}{i}=OpsOut{N}{i}+plus(AssignOps{N}{i}{j-1},AssignOps{N}{i}{j});
                OpsOut{N}{i}=OpsOut{N}{i}+AssignOps{N}{i}{j};
            end
        end
    else
        for  N= 1:size(Vals,2)
            OpsOut{N}{i}=cell2mat(AssignOps{N}{i});
        end
    end
end


if isstruct(SysFixed)
    for  N= size(Vals,2):-1:1
        Sys_Out(N) = cell2struct(OpsOut{N}',f,1);
    end
    ff = fieldnames(SysFixed);
    for  N= size(Vals,2):-1:1
        for i = 1:length(ff)
            if ~isfield(Sys_Out(N),ff{i})||isempty(Sys_Out(N).(ff{i}))
                Sys_Out(N).(ff{i}) = SysFixed.(ff{i});
            else
            Sys_Out(N).(ff{i})(logical(SysFixed.(ff{i}))) = SysFixed.(ff{i})(logical(SysFixed.(ff{i})));
            end
        end
    end
else
    for  N= size(Vals,2):-1:1
        Sys_Out(N) = cell2struct(OpsOut{N}',f,1);
    end
    for  N= size(Vals,2):-1:1
        Sys_Out(N).S = S;

        Sys_Out = rmfield(Sys_Out, f(structfun(@isempty, Sys_Out(N))));

    end
end


% f = fieldnames(Ops);
% COps = struct2cell(Ops);
%
%     k=1;
%     for i=1:length(COps)       %Stevens Operators
%         if isempty(COps{i})
%             continue
%         end
%         for j=1:length(COps{i})
%
%             COps{i}{j} = COps{i}{j}*Y(k);
%             k=j+1;
%         end
%         if length(COps{i})>1
%             COps{i}=plus(COps{i}{:});
%         else
%             COps{i}=cell2mat(COps{i});
%         end
%     end
%     Sys_Out{N} = cell2struct(COps,f,1);
