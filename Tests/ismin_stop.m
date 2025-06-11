function ok = ismin_stop()
f = 1e-18;
x = rand(6,1);
p = rand(6,1)*1e-6;
gradf = rand(6,1)*1e-6;
NIter = 1001;

ConvergenceParams.FunctionTolerance = 1e-8;
ConvergenceParams.StepTolerance = 1e-4;
ConvergenceParams.GradientTolerance = 1e-4;
ConvergenceParams.MaximumIterations = 1000;
ConvergenceParams.RelativeStepTolerance = 1e-4;


ok(1) = ~ismin(1e10,1e10,1e10,1e10,0,ConvergenceParams);
ok(2) = ismin(1e10,1e10,1e10,1e10,0,ConvergenceParams, [Inf,NaN]);
ok(3) = ismin(1e10,1e10,p,1e10,0,ConvergenceParams);
ok(4) = ismin(inf,1e10,1e10,1e10,0,ConvergenceParams);
ok(5) = ismin(f,1e10,1e10,1e10,0,ConvergenceParams);
ok(6) = ismin(1e10,1e10,1e10,gradf,0,ConvergenceParams);
ok(7) = ismin(1e10,1e10,1e10,1e10,NIter,ConvergenceParams);
ok(8) = ismin(1e10,x,p,1e10,0,ConvergenceParams);


