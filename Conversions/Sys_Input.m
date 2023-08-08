function [A,A0,xint,Ops,SysFixed] = Sys_Input(Sys0,Vary)
xint=[];
A={};
Ops={};
N_electrons = length(Sys0.S);
n = prod(2*Sys0.S+1);
Varynames = fieldnames(Vary);
Sys0names = fieldnames(Sys0);
A0=sparse(n,n);

if ~any(contains(Sys0names,'J'))&&~any(contains(Sys0names,'ee'))&&N_electrons>1
    error('Please provide an electron-electron exchange coupling. Use either the J or ee field')
end
stevk="";
for i = 1:length(Varynames)
    if Varynames{i}(1) =='B'
        if ~any(contains(Sys0names,Varynames{i})) %validates that Sy0 and Vary have same fields.
            error('Please input the same Stevens opertors in the Sys0 and Vary structures')
        end
        if (str2double(Varynames{i}(2:end)))>12
            error('easyspin only accepts Stevens operators up to degree 12')
        end
        %         stevk(end+1) = str2num(Varynames{i}(2:end));    %lists stevens operators used
        stevk(length(stevk)+1) = (Varynames{i});    %lists stevens operators used
    end
end
stevk=stevk(2:end); %removes first, empty string
for Bk = stevk
    SysFixed.(Bk) = zeros(size(Sys0.(Bk)));
    SysFixed.(Bk)(~logical((Vary.(Bk)))) = Sys0.(Bk)(~logical((Vary.(Bk))));
    Vary.(Bk) = logical(Vary.(Bk)); %will use uncertainty data later
    [TempOps,TempInt,TempA] = StevOps(Sys0,Vary,Bk);
    Ops.(Bk) = TempOps.(Bk);
    if size(TempInt,2)>1
        xint =[xint;TempInt'];
    else
        xint =[xint;TempInt];
    end
    A = [A,TempA];
end




if  any(contains(Sys0names,'J'))
    if any(contains(Sys0names,'ee'))
        error('Please only use one of the J or ee fields to input the exchange interaction')
    end
    SysFixed.J = zeros(size(Sys0.J));
    SysFixed.J(~logical((Vary.J))) = Sys0.J(~logical((Vary.J)));
    [TempOps, TempInt, TempA] = JOps(Sys0,Vary);
    Ops.J = TempOps.J;
    if size(TempInt,2)>1
        xint =[xint;TempInt'];
    else
        xint =[xint;TempInt];
    end
    A = [A,TempA];
elseif any(contains(Sys0names,'ee'))    %Used Sys0.ee
    SysFixed.ee = zeros(size(Sys0.ee));
    SysFixed.ee(~logical((Vary.ee))) = Sys0.ee(~logical((Vary.ee)));

    [TempOps, TempInt, TempA] = eeOps(Sys0,Vary);
    Ops.ee = TempOps.ee;
    if size(TempInt,2)>1
        xint =[xint;TempInt'];
    else
        xint =[xint;TempInt];
    end
    A = [A,TempA];
end


SysFixed.S = Sys0.S;
A0=ham(SysFixed,[0 0 0],'sparse');
end


function [Ops, Int, A] = eeOps(Sys0,Vary)
Ops.ee = []; A = []; Int = [];
Sys.S = Sys0.S;
N_electrons = length(Sys.S);
% if all(size(Sys0.ee) == [N_electrons,3])||all(size(Sys0.ee) == [3*N_electrons,3])
[initials] = unique(abs(Sys0.ee),'stable');
initials = initials';
if all(size(Sys0.ee) == [N_electrons,3])

    OpDim= [3,(N_electrons-1)*(N_electrons)/2];
elseif all(size(Sys0.ee) == [3*(N_electrons-1)*(N_electrons)/2,3])
    OpDim=[3*(N_electrons-1)*(N_electrons)/2,3];
else
    error('Please input an ee field that is Nx3 or 3Nx3')
end
%          idx=reshape(idx,size(Sys0.ee));
for j = 1:length(initials)  %iterate along row
    %             pidx = ((idx==j).*Vary.ee)>0;
    %             nidx = ((idx==j).*Vary.ee)<0;

    idx=zeros(OpDim);
    idx(logical(((Sys0.ee== -initials(j))).*Vary.ee)) = -1;
    idx(logical(((Sys0.ee== initials(j))).*Vary.ee)) = 1;    %the pinned entries
    if any(any(idx))
        if initials(j) == 0
            warning("Any ee parameters initialised to 0 will not be varied.")
            continue
%         elseif  length(unique(abs(Sys0.ee(logical(idx)))))>1
%             error("")
        end
        
        [Ops.ee{end+1}, A{end+1}] = CalculateOp(Sys,'ee',"ham_ee",OpDim,idx);
        Int(end+1) = initials(j);
    end
end
end


function [Op, A] = CalculateOp(Sys,FN,OpType,OpDim,idx,col)
if nargin>6
    error("can only parse a column or row")
end
pidx = idx>0;
nidx = idx<0;
if nargin>5
    Op = zeros(OpDim);
    Op(pidx,col) = 1;
    Op(nidx,col) = -1;
    Sys.(FN) = Op;
    A= feval(OpType,Sys,1:length(Sys.S),'sparse');
else
    Op = zeros(OpDim);
    Op(pidx) = 1;
    Op(nidx) = -1;
    Sys.(FN) = Op;
    A= feval(OpType,Sys,1:length(Sys.S),'sparse');
end
end


function [Ops, Int, A] = StevOps(Sys0,Vary,Bk)
A=[];   Ops.(Bk)=[];    Int=[];
Sys.S = Sys0.S;
N_electrons = length(Sys.S);
Sys.J = zeros(1,(N_electrons-1)*(N_electrons)/2);
for i = 1:size(Sys0.(Bk),2)    %Iterate through columns
    [initials,~, idx] = unique(Sys0.(Bk)(:,i),'stable');
    initials = initials';
    for j = 1:length(initials)
        lidx=((idx==j).*Vary.(Bk)(:,i));
        if any(logical(lidx))
            [Ops.(Bk){end+1}, A{end+1}] = CalculateOp(Sys,Bk,"ham_zf",size(Sys0.(Bk)),lidx,i);
            Int(end+1) = initials(j);
        end
    end
end
end
function [Ops, Int, A] = JOps(Sys0,Vary)
Ops.J = []; A = []; Int = [];
Sys.S = Sys0.S;
N_electrons = length(Sys.S);
[initials,~, idx] = unique(Sys0.J,'stable');
initials = initials';
if size(idx,1)>1
    idx = idx';
end
for j = 1:length(initials)  %iterate along row
    lidx = ((idx==j).*Vary.J);
    if any(logical(lidx))
        [Ops.J{end+1}, A{end+1}] = CalculateOp(Sys,'J',"ham_ee",[1,(N_electrons-1)*(N_electrons)/2],lidx);
        Int(end+1) = initials(j);

    end
end
end

%
% Sys = Cr_Spin_Sys_3(3);
% Sys.B2(1) = 1.7*1e3;
% Sys.B2(1,3) = -3.2*1e3;
% Sys.B2(3) = 0;
% Sys.ee=[ones(1,3)*3.53*1e5;0 0 0; ones(1,3)*-3.53*1e5];
% Sys.ee=[eye(3)*3.53*1e5;zeros(3); eye(3)*-3.53*1e5];
%
% % Sys = rmfield(Sys, 'ee'); Sys.J = [3.53*1e5, 0, 3.53*1e5]; Vary.J = [1,1,1];
% Vary = Sys;
% Vary.B2(:,[1,3])=1;
%

