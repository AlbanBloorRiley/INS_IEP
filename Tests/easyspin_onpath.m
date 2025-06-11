function ok = easyspin_onpath()
EasySpinPath = fileparts(which('easyspin'));
ok = true;
if isempty(EasySpinPath)
    ok = false;
end