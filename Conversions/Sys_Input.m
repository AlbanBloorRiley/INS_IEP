function [A,A0,xint,Ops,SysFixed] = Sys_Input(Sys0,Vary)
%   SYS_INPUT calculates the basis of matrices(A) and initial values (xint) to be
%   varied, as well as the matrix of fixed values (A0) as well as the 
% corresponding spin structure (SysFixed).
%
%   Inputs:
%   Sys0 - easyspin structure containing the current guess of the spin system.  
%   Vary - logical structure defining which variables in Sys0 to vary.
%   
%   Outputs:
%   A - Cell structure of basis matrices
%   A0 = Matrix determined by fixed variables
%   xint - Initial values as defined by Sys0
%   Ops - Cell structure containing basis of all the spin structure fields,
% each corresponding to the equivalent element of xint.
%   SysFixed - the spin structure formed of only fixed variables.

xint=[];
A={};
Ops={};
N_electrons = length(Sys0.S);
Varynames = fieldnames(Vary);
Sys0names = fieldnames(Sys0);

if ~all(contains(Sys0names,[Varynames;"S"]))
    error("Please use all the same fields in 'Sys0' and 'Vary' inputs.")
end

if ~any(contains(Sys0names,'J'))&&~any(contains(Sys0names,'ee'))&&N_electrons>1
    error('Please provide an electron-electron exchange coupling. Use either the J or ee field')
end
stevk="";   %adds empty string
for i = 1:length(Varynames)
    if Varynames{i}(1) =='B'
        if ~any(contains(Sys0names,Varynames{i})) %validates that Sy0 and Vary have same fields.
            error('Please input the same Stevens opertors in the Sys0 and Vary structures')
        end
        if (str2double(Varynames{i}(2:end)))>12
            error('easyspin only accepts Stevens operators up to degree 12')
        end
        stevk(length(stevk)+1) = (Varynames{i});    %lists stevens operators used
    end
end
stevk=stevk(2:end); %removes first, empty string

for Bk = stevk    %Iterates through all the Stevens operators used 
    SysFixed.(Bk) = zeros(size(Sys0.(Bk)));     %Initialises fixed values to 0                                                  
    SysFixed.(Bk)(~logical((Vary.(Bk)))) = Sys0.(Bk)(~logical((Vary.(Bk))));    %saves the fixed values
    %----------------------------------------%
    %Could use the values to define a range like in easyspin.
    Vary.(Bk) = logical(Vary.(Bk)); 
    %----------------------------------------%
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
    Vary.J = logical(Vary.J);
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
    Vary.ee = logical(Vary.ee);   %Does this break anisotropy?
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
if all(size(Sys0.ee) == [(N_electrons-1)*(N_electrons)/2,3])

    OpDim= [(N_electrons-1)*(N_electrons)/2,3];
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
        if initials(j) == 0
            initials(j) = 1e-18;
        end
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
            if initials(j) == 0
                initials(j) = 1e-18;
            end
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
        if initials(j) == 0
            initials(j) = 1e-18;
        end
        Int(end+1) = initials(j);

    end
end
end


