//+
f = 4;
//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Line(5) = {1, 3};
//+
MeshSize {1, 3} = 0.01 / f;
//+
MeshSize {2, 4} = 0.1 / f;
//+
Field[1] = Distance;
//+
Field[1].CurvesList = {5};
//+
Field[1].Sampling = 100;
//+
Field[2] = Threshold;
//+
Field[2].DistMax = 0.15;
//+
Field[2].DistMin = 0.02;
//+
Field[2].InField = 1;
//+
Field[2].SizeMax = 0.1 / f;
//+
Field[2].SizeMin = 0.01 / f;
//+
Background Field = 2;
