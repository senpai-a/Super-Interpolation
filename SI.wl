(*清洗样本,将宽高为奇数的图像裁切一行或一列使其成为偶数宽高,将多通道图像分离*)
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

(*从HR图像得到LR图像，宽高减半*)
getLR[img_/;ImageQ[img]]:=Module[{w,h},
w=ImageDimensions[img][[1]];
h=ImageDimensions[img][[2]];
ImageResize[img,{w/2,h/2},Resampling->"Cubic"]
];

(*输入HR img，输出{2x2 HR patches, 3x3 LR patches}*)
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
h=(1	1
-1	-1

)p;
v=(1	-1
1	-1

)p;
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

(*输入3x3 LR Patch，计算EO Index = 0..624*)
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
(*输入一个EOClass的{HR-LR Patches}，使用L2正则化系数\[Lambda]计算MSE意义最佳线性映射*)
findLinearMapping[patches_/;ListQ[patches],\[Lambda]_/;NumberQ[\[Lambda]]]:=Module[{Xc,Yc},
If[Length[patches]==0,
(.25	.25	0	.25	.25	0	0	0	0
0	.25	.25	0	.25	.25	0	0	0
0	0	0	.25	.25	0	.25	.25	0
0	0	0	0	.25	.25	0	.25	.25

),
Xc=ParallelMap[#[[1]]&,patches];
Yc=ParallelMap[#[[2]]&,patches];
Xc=ParallelMap[N[ImageData[#,"Byte"]]&,Xc];
Yc=ParallelMap[N[ImageData[#,"Byte"]]&,Yc];
Xc=Transpose[Flatten/@Xc];
Yc=Transpose[Flatten/@Yc];
solveMap[Xc,Yc,\[Lambda]]
]
];

(*从已清洗的HR图像样本列表训练映射，输出625个矩阵的列表*)
train[hrImgList_/;ListQ[hrImgList],\[Lambda]_/;NumberQ[\[Lambda]]]:=Module[{patches,classes},
patches=Flatten[ParallelMap[getPatches,hrImgList],1];
patches=getEOIndex/@patches;
classes=ParallelTable[Select[patches,#[[3]]==i&],{i,0,624}];
Map[findLinearMapping[#,\[Lambda]]&,classes]
];

cdot=Compile[{{M,_Real,2},{y,_Integer,1}},y[[1]] M[[1,1]]+y[[2]] M[[1,2]]+y[[3]] M[[1,3]]+y[[4]] M[[1,4]]+y[[5]] M[[1,5]]+y[[6]] M[[1,6]]+y[[7]] M[[1,7]]+y[[8]] M[[1,8]]+y[[9]] M[[1,9]]+y[[1]] M[[2,1]]+y[[2]] M[[2,2]]+y[[3]] M[[2,3]]+y[[4]] M[[2,4]]+y[[5]] M[[2,5]]+y[[6]] M[[2,6]]+y[[7]] M[[2,7]]+y[[8]] M[[2,8]]+y[[9]] M[[2,9]]+y[[1]] M[[3,1]]+y[[2]] M[[3,2]]+y[[3]] M[[3,3]]+y[[4]] M[[3,4]]+y[[5]] M[[3,5]]+y[[6]] M[[3,6]]+y[[7]] M[[3,7]]+y[[8]] M[[3,8]]+y[[9]] M[[3,9]]+y[[1]] M[[4,1]]+y[[2]] M[[4,2]]+y[[3]] M[[4,3]]+y[[4]] M[[4,4]]+y[[5]] M[[4,5]]+y[[6]] M[[4,6]]+y[[7]] M[[4,7]]+y[[8]] M[[4,8]]+y[[9]] M[[4,9]],CompilationTarget->"C"];
(*对单通道lrImg应用2xSR*)
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

(*对图像img应用2xSR*)
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

(*从指定目录加载所有图像，并清洗数据得到训练样本图像的列表*)
getSamples[dir_/;DirectoryQ[dir]||DirectoryQ[NotebookDirectory[]<>dir]]:=Module[{images},
SetDirectory[If[DirectoryQ[dir],dir,NotebookDirectory[]<>dir]];
images=Select[Import[#]&/@FileNames[],ImageQ];
Flatten[regularImg/@images]
];