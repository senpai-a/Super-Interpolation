(* ::Package:: *)

(*\:6e05\:6d17\:6837\:672c,\:5c06\:5bbd\:9ad8\:4e3a\:5947\:6570\:7684\:56fe\:50cf\:88c1\:5207\:4e00\:884c\:6216\:4e00\:5217\:4f7f\:5176\:6210\:4e3a\:5076\:6570\:5bbd\:9ad8,\:5c06\:591a\:901a\:9053\:56fe\:50cf\:5206\:79bb*)
regularImg[img_/;ImageQ[img]]:=Module[{w,h,sized},
w=ImageDimensions[img][[1]];
h=ImageDimensions[img][[2]];
If[w>=6&&h>=6,
sized=ImageCrop[img,{w-Mod[w,2],h-Mod[h,2]},{Left,Top}];
If[ImageChannels[sized]>1,
Image[#,"Byte"]&/@ImageData[sized,"Byte",Interleaving->False],
sized
]
]
];

(*\:4eceHR\:56fe\:50cf\:5f97\:5230LR\:56fe\:50cf\:ff0c\:5bbd\:9ad8\:51cf\:534a*)
getLR[img_/;ImageQ[img]]:=Module[{w,h},
w=ImageDimensions[img][[1]];
h=ImageDimensions[img][[2]];
ImageResize[img,{w/2,h/2},Resampling->"Cubic"]
];

(*\:8f93\:5165HR img\:ff0c\:8f93\:51fa{2x2 HR patches, 3x3 LR patches}*)
getPatches[hrImg_/;ImageQ[hrImg]]:=Module[{w,h,lrImg,hrImg1,hrPatches,lrPatches},
w=ImageDimensions[hrImg][[1]];
h=ImageDimensions[hrImg][[2]];
lrImg=getLR[hrImg];
hrImg1=ImageCrop[hrImg,{w-4,h-4},Center];
hrPatches=Flatten[ImagePartition[hrImg1,{2,2}]];
lrPatches=Flatten[ImagePartition[lrImg,{3,3},{1,1}]];
MapThread[List,{hrPatches,lrPatches}]
];

EO=Compile[{{p,_Integer,2}},
Module[{h,v,gh,gv,len,dir},
h=({
 {1, 1},
 {-1, -1}
})p;
v=({
 {1, -1},
 {1, -1}
})p;
gh=Total[Flatten[h]];
gv=Total[Flatten[v]];
len=N[Sqrt[gh*gh+gv*gv]];
dir=If[gv!=0,
N[ArcTan[gh/gv]],
If[gh>0,
\[Pi]/2,
-\[Pi]/2
]
];
If[len<15,
0,
Ceiling[Mod[4./\[Pi] (dir+\[Pi]/8.),4.]]
]
],
CompilationTarget->"C",
Parallelization->True
];

(*\:8f93\:51653x3 LR Patch\:ff0c\:8ba1\:7b97EO Index = 0..624*)
EOIndex[lrPatch_/;ImageQ[lrPatch]]:=Module[{parts},
parts=ImagePartition[lrPatch,{2,2},{1,1}];
parts=Flatten[parts];
parts=ParallelMap[ImageData[#,"Byte"]&,parts];
Total[{125,25,5,1} ParallelMap[EO,parts]]
];
getEOIndex[patchPair_]:=Append[patchPair,EOIndex[patchPair[[2]]]];

solveMap=Compile[{{Xc,_Integer,2},{Yc,_Integer,2},{\[Lambda],_Real}},
Xc.Yc\[Transpose].Inverse[Yc.Yc\[Transpose]+\[Lambda] IdentityMatrix[9]]
];
(*\:8f93\:5165\:4e00\:4e2aEOClass\:7684{HR-LR Patches}\:ff0c\:4f7f\:7528L2\:6b63\:5219\:5316\:7cfb\:6570\[Lambda]\:8ba1\:7b97MSE\:610f\:4e49\:6700\:4f73\:7ebf\:6027\:6620\:5c04*)
findLinearMapping[patches_/;ListQ[patches],\[Lambda]_/;NumberQ[\[Lambda]]]:=Module[{Xc,Yc},
If[Length[patches]==0,
({
 {.25, .25, 0, .25, .25, 0, 0, 0, 0},
 {0, .25, .25, 0, .25, .25, 0, 0, 0},
 {0, 0, 0, .25, .25, 0, .25, .25, 0},
 {0, 0, 0, 0, .25, .25, 0, .25, .25}
}),
Xc=ParallelMap[#[[1]]&,patches];
Yc=ParallelMap[#[[2]]&,patches];
Xc=ParallelMap[N[ImageData[#,"Byte"]]&,Xc];
Yc=ParallelMap[N[ImageData[#,"Byte"]]&,Yc];
Xc=Transpose[Flatten/@Xc];
Yc=Transpose[Flatten/@Yc];
solveMap[Xc,Yc,\[Lambda]]
]
];

(*\:4ece\:5df2\:6e05\:6d17\:7684HR\:56fe\:50cf\:6837\:672c\:5217\:8868\:8bad\:7ec3\:6620\:5c04\:ff0c\:8f93\:51fa625\:4e2a\:77e9\:9635\:7684\:5217\:8868*)
train[hrImgList_/;ListQ[hrImgList],\[Lambda]_/;NumberQ[\[Lambda]]]:=Module[{patches,classes},
patches=Flatten[ParallelMap[getPatches,hrImgList],1];
patches=getEOIndex/@patches;
classes=ParallelTable[Select[patches,#[[3]]==i&],{i,0,624}];
Map[findLinearMapping[#,\[Lambda]]&,classes]
];

cdot=Compile[{{M,_Real,2},{y,_Integer,1}},y[[1]] M[[1,1]]+y[[2]] M[[1,2]]+y[[3]] M[[1,3]]+y[[4]] M[[1,4]]+y[[5]] M[[1,5]]+y[[6]] M[[1,6]]+y[[7]] M[[1,7]]+y[[8]] M[[1,8]]+y[[9]] M[[1,9]]+y[[1]] M[[2,1]]+y[[2]] M[[2,2]]+y[[3]] M[[2,3]]+y[[4]] M[[2,4]]+y[[5]] M[[2,5]]+y[[6]] M[[2,6]]+y[[7]] M[[2,7]]+y[[8]] M[[2,8]]+y[[9]] M[[2,9]]+y[[1]] M[[3,1]]+y[[2]] M[[3,2]]+y[[3]] M[[3,3]]+y[[4]] M[[3,4]]+y[[5]] M[[3,5]]+y[[6]] M[[3,6]]+y[[7]] M[[3,7]]+y[[8]] M[[3,8]]+y[[9]] M[[3,9]]+y[[1]] M[[4,1]]+y[[2]] M[[4,2]]+y[[3]] M[[4,3]]+y[[4]] M[[4,4]]+y[[5]] M[[4,5]]+y[[6]] M[[4,6]]+y[[7]] M[[4,7]]+y[[8]] M[[4,8]]+y[[9]] M[[4,9]],CompilationTarget->"C"];
(*\:5bf9\:5355\:901a\:9053lrImg\:5e94\:75282xSR*)
SR[mappings_/;ListQ[mappings]&&AllTrue[mappings,MatrixQ[#]&&Dimensions[#]=={4,9}&],
lrImg_/;ImageQ[lrImg]]:=Module[{lrPatches,map,index,M,y,x},
lrPatches=ImagePartition[lrImg,{3,3},{1,1}];

map=Function[{lrPatch},
index=EOIndex[lrPatch];
M=mappings[[index+1]];
y=Flatten[ImageData[lrPatch,"Byte"]];
x=cdot[M,y];
Image[Partition[x,2],"Byte"]
];
ImageAssemble[Map[map,lrPatches,{2}]]
];

(*\:5bf9\:56fe\:50cfimg\:5e94\:75282xSR*)
imgSR[mappings_/;ListQ[mappings]&&AllTrue[mappings,MatrixQ[#]&&Dimensions[#]=={4,9}&],
img_/;ImageQ[img]]:=Module[{rgb},
If[ImageChannels[img]==1,

SR[mappings,img],

rgb=Image[#,"Byte"]&/@ImageData[img,"Byte",Interleaving->False];
rgb=SR[mappings,#]&/@rgb;
rgb=ImageData[#,"Byte"]&/@rgb;
Image[rgb,"Byte",ColorSpace->"RGB",Interleaving->False]
]
];

(*\:4ece\:6307\:5b9a\:76ee\:5f55\:52a0\:8f7d\:6240\:6709\:56fe\:50cf\:ff0c\:5e76\:6e05\:6d17\:6570\:636e\:5f97\:5230\:8bad\:7ec3\:6837\:672c\:56fe\:50cf\:7684\:5217\:8868*)
getSamples[dir_/;DirectoryQ[dir]||DirectoryQ[NotebookDirectory[]<>dir]]:=Module[{images},
SetDirectory[If[DirectoryQ[dir],dir,NotebookDirectory[]<>dir]];
images=Select[Import[#]&/@FileNames[],ImageQ];
Flatten[regularImg/@images]
];
