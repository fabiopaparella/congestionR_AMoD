clear all
clc
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultaxesticklabelinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
symmetricFlag = 1;

N = 357;
global RoadGraph NodesLoc;
RoadGraph=cell(N,1);

RoadGraph{1} = [2 3];
RoadGraph{2} = [1 4 7];
RoadGraph{3} = [1 4 6];
RoadGraph{4} = [2 3 5];
RoadGraph{5} = [6 4 18];
RoadGraph{6} = [3 24 5];
RoadGraph{7} = [2 10];
RoadGraph{8} = [10 9];
RoadGraph{9} = [8 17];
RoadGraph{10} = [8 7 11];
RoadGraph{11} = [12 9 13];
RoadGraph{12} = [7];
RoadGraph{13} = [14 15];
RoadGraph{14} = [12 18];
RoadGraph{15} = [17 21];
RoadGraph{16} = [14 15];
RoadGraph{17} = [9 23];
RoadGraph{18} = [16 19 5];
RoadGraph{19} = [18 20 28];
RoadGraph{20} = [19 16 21];
RoadGraph{21} = [20 22 31];
RoadGraph{22} = [23 37 21];
RoadGraph{23} = [17 38 22];
RoadGraph{24} = [6 59 25];
RoadGraph{25} = [24 26 56 58];
RoadGraph{26} = [27 25];
RoadGraph{27} = [5 53 26 28];
RoadGraph{28} = [27 29 50];
RoadGraph{29} = [28 30 19];
RoadGraph{30} = [20 29 31];
RoadGraph{31} = [30 33 32];
RoadGraph{32} = [31 22];
RoadGraph{33} = [35 34];
RoadGraph{34} = [46 45];
RoadGraph{35} = [44 34];
RoadGraph{36} = [35 32 35];
RoadGraph{37} = [36 39];
RoadGraph{38} = [23 40];
RoadGraph{39} = [38 41];
RoadGraph{40} = [38 41 91];
RoadGraph{41} = [40 42 87];
RoadGraph{42} = [41 44 43];
RoadGraph{43} = [36 44];
RoadGraph{44} = [42 45 83];
RoadGraph{45} = [44 47];
RoadGraph{46} = [48 30];
RoadGraph{47} = [45 49 46];
RoadGraph{48} = [50 29];
RoadGraph{49} = [48 47 51];
RoadGraph{50} = [51 52];
RoadGraph{51} = [49 52 76];
RoadGraph{52} = [51 53];
RoadGraph{53} = [52 64 54 27];
RoadGraph{54} = [58 55 26];
RoadGraph{55} = [54 56 63];
RoadGraph{56} = [57 67 25];
RoadGraph{57} = [60 58 68];
RoadGraph{58} = [57 25 61];
RoadGraph{59} = [24 61 62];
RoadGraph{60} = [62 61 57 71];
RoadGraph{61} = [60 59 58];
RoadGraph{62} = [60 59 69];
RoadGraph{63} = [66 64 55];
RoadGraph{64} = [76 65 53];
RoadGraph{65} = [75 77 66 64];
RoadGraph{66} = [65 67 63 74];
RoadGraph{67} = [66 68 73 56];
RoadGraph{68} = [72 57 67 69];
RoadGraph{69} = [62 70 68];
RoadGraph{70} = [69 92 71];
RoadGraph{71} = [72 70 93];
RoadGraph{72} = [68 71 73];
RoadGraph{73} = [67 94 72 74];
RoadGraph{74} = [66 75 73];
RoadGraph{75} = [65 97 90 74];
RoadGraph{76} = [80 77];
RoadGraph{77} = [78 79 65];
RoadGraph{78} = [79 90];
RoadGraph{79} = [80 77];
RoadGraph{80} = [81 49];
RoadGraph{81} = [47 82];
RoadGraph{82} = [83 84];
RoadGraph{83} = [84 85];
RoadGraph{84} = [85 88];
RoadGraph{85} = [42];
RoadGraph{86} = [87 85];
RoadGraph{87} = [109 91];
RoadGraph{88} = [105 86 89];
RoadGraph{89} = [88 90 81];
RoadGraph{90} = [99 89 75];
RoadGraph{91} = [40 113 87];
RoadGraph{92} = [70 115];
RoadGraph{93} = [71 94 115];
RoadGraph{94} = [117 93 95];
RoadGraph{95} = [74 94 96];
RoadGraph{96} = [98 121 95 102];
RoadGraph{97} = [75 98 99];
RoadGraph{98} = [96 97];
RoadGraph{99} = [98 102];
RoadGraph{100} = [99 89];
RoadGraph{101} = [100];
RoadGraph{102} = [96 103 123];
RoadGraph{103} = [102 104 100];
RoadGraph{104} = [101 103 106];
RoadGraph{105} = [101 106 107];
RoadGraph{106} = [126 104 108];
RoadGraph{107} = [105 110 86];
RoadGraph{108} = [106 111 107];
RoadGraph{109} = [110 113 107];
RoadGraph{110} = [111 107];
RoadGraph{111} = [112 128 108];
RoadGraph{112} = [109 114 111];
RoadGraph{113} = [91 114 109];
RoadGraph{114} = [113 131 112 130];
RoadGraph{115} = [93 92 116 118];
RoadGraph{116} = [115 118];
RoadGraph{117} = [116 119];
RoadGraph{118} = [116 148 119];
RoadGraph{119} = [120 147 118];
RoadGraph{120} = [95 121 119];
RoadGraph{121} = [145 122 120];
RoadGraph{122} = [123 121];
RoadGraph{123} = [103 124 122 143];
RoadGraph{124} = [142 123 125];
RoadGraph{125} = [104 126 124];
RoadGraph{126} = [140 125 127];
RoadGraph{127} = [108 128 126];
RoadGraph{128} = [138 127 129];
RoadGraph{129} = [112 128 130];
RoadGraph{130} = [136 129 132];
RoadGraph{131} = [132 133 114];
RoadGraph{132} = [133 131 130];
RoadGraph{133} = [131 132 134];
RoadGraph{134} = [135 133 167];
RoadGraph{135} = [132 134 136];
RoadGraph{136} = [164 135 137];
RoadGraph{137} = [129 136 138];
RoadGraph{138} = [159 137 139];
RoadGraph{139} = [127 138 140];
RoadGraph{140} = [157 139 141];
RoadGraph{141} = [125 142 140];
RoadGraph{142} = [141 143 155];
RoadGraph{143} = [123 154 142 144];
RoadGraph{144} = [143 145 122];
RoadGraph{145} = [152 144 146];
RoadGraph{146} = [120 145 147];
RoadGraph{147} = [150 146 148];
RoadGraph{148} = [147 118 149];
RoadGraph{149} = [148 150 184];
RoadGraph{150} = [149 151 183];
RoadGraph{151} = [146 150 152];
RoadGraph{152} = [179 151 153];
RoadGraph{153} = [144 152 154];
RoadGraph{154} = [177 143 155 153];
RoadGraph{155} = [176 156 154];
RoadGraph{156} = [141 157 155];
RoadGraph{157} = [174 158 156];
RoadGraph{158} = [139 159 157];
RoadGraph{159} = [172 160 158];
RoadGraph{160} = [161 163 159];
RoadGraph{161} = [137 162];
RoadGraph{162} = [164 163];
RoadGraph{163} = [160 165];
RoadGraph{164} = [165];
RoadGraph{165} = [170 163 166];
RoadGraph{166} = [135 165 167];
RoadGraph{167} = [134 166 168];
RoadGraph{168} = [167 169 237];
RoadGraph{169} = [168 166 170];
RoadGraph{170} = [205 169 171];
RoadGraph{171} = [170 160 172];
RoadGraph{172} = [203 171 173];
RoadGraph{173} = [158 172 174];
RoadGraph{174} = [201 173 175];
RoadGraph{175} = [156 174 176];
RoadGraph{176} = [197 175 177];
RoadGraph{177} = [154 195 176 178];
RoadGraph{178} = [153 177 179];
RoadGraph{179} = [191 178 180];
RoadGraph{180} = [151 179 181];
RoadGraph{181} = [180 187];
RoadGraph{182} = [181];
RoadGraph{183} = [182 184];
RoadGraph{184} = [149 183 185];
RoadGraph{185} = [182 184 186];
RoadGraph{186} = [227 187 185];
RoadGraph{187} = [186 188 226];
RoadGraph{188} = [190 187];
RoadGraph{189} = [188 180];
RoadGraph{190} = [189 192];
RoadGraph{191} = [189 192];
RoadGraph{192} = [220 194];
RoadGraph{193} = [191 178];
RoadGraph{194} = [193 196];
RoadGraph{195} = [177 196 197];
RoadGraph{196} = [195 216 198];
RoadGraph{197} = [195 198];
RoadGraph{198} = [214 200];
RoadGraph{199} = [197 175];
RoadGraph{200} = [199];
RoadGraph{201} = [199 202];
RoadGraph{202} = [173 201 203];
RoadGraph{203} = [202 204 210];
RoadGraph{204} = [208 210 203 171];
RoadGraph{205} = [204 207];
RoadGraph{206} = [169 205];
RoadGraph{207} = [235];
RoadGraph{208} = [207 204 235];
RoadGraph{209} = [204 208];
RoadGraph{210} = [212 211];
RoadGraph{211} = [209 233 210];
RoadGraph{212} = [200 214];
RoadGraph{213} = [211 212];
RoadGraph{214} = [215 216];
RoadGraph{215} = [231 213];
RoadGraph{216} = [218 196 217];
RoadGraph{217} = [215 216 230];
RoadGraph{218} = [194 220];
RoadGraph{219} = [218 217];
RoadGraph{220} = [221];
RoadGraph{221} = [219 228];
RoadGraph{222} = [225 228 223];
RoadGraph{223} = [222 224 252];
RoadGraph{224} = [253 227 223];
RoadGraph{225} = [190];
RoadGraph{226} = [225 223];
RoadGraph{227} = [226 186 224];
RoadGraph{228} = [250 229 222];
RoadGraph{229} = [219 230 228];
RoadGraph{230} = [217 248 231 229];
RoadGraph{231} = [247 232 230];
RoadGraph{232} = [213 233 231];
RoadGraph{233} = [241 234 211];
RoadGraph{234} = [233 209 235];
RoadGraph{235} = [236 208 234 238 239];
RoadGraph{236} = [206 235];
RoadGraph{237} = [168 268 267];
RoadGraph{238} = [268 265 239 235];
RoadGraph{239} = [245 238 244];
RoadGraph{240} = [234];
RoadGraph{241} = [240 242 233];
RoadGraph{242} = [243 246 262 241];
RoadGraph{243} = [242 244];
RoadGraph{244} = [240 239];
RoadGraph{245} = [243 264];
RoadGraph{246} = [232 247 242];
RoadGraph{247} = [259 246 248];
RoadGraph{248} = [247 249 258 230];
RoadGraph{249} = [229 248 250];
RoadGraph{250} = [256 249 251];
RoadGraph{251} = [222 250 252];
RoadGraph{252} = [254 251 253];
RoadGraph{253} = [224 252 286];
RoadGraph{254} = [285 255];
RoadGraph{255} = [254 256 251];
RoadGraph{256} = [283 255 257];
RoadGraph{257} = [249 256 258];
RoadGraph{258} = [248 280 257 259];
RoadGraph{259} = [278 260 258];
RoadGraph{260} = [261 259];
RoadGraph{261} = [246 262];
RoadGraph{262} = [261 263 275 242];
RoadGraph{263} = [243 262 264];
RoadGraph{264} = [273 265 263];
RoadGraph{265} = [238 264 266 271];
RoadGraph{266} = [268 269 265];
RoadGraph{267} = [237 269];
RoadGraph{268} = [237 266 236 238];
RoadGraph{269} = [266 271 270];
RoadGraph{270} = [267 269 272 349];
RoadGraph{271} = [265 272];
RoadGraph{272} = [270 273 271 304];
RoadGraph{273} = [303 272 274];
RoadGraph{274} = [263 273 275];
RoadGraph{275} = [301 274 276 262];
RoadGraph{276} = [278 260];
RoadGraph{277} = [276 275];
RoadGraph{278} = [276 280 279];
RoadGraph{279} = [277 296];
RoadGraph{280} = [278 282 281 258];
RoadGraph{281} = [280 279 295];
RoadGraph{282} = [280 283 257];
RoadGraph{283} = [293 282 284];
RoadGraph{284} = [255 283 285];
RoadGraph{285} = [288 286 284];
RoadGraph{286} = [287 253 285];
RoadGraph{287} = [286 290 288];
RoadGraph{288} = [291 289];
RoadGraph{289} = [284];
RoadGraph{290} = [287 318 291];
RoadGraph{291} = [317 290 292];
RoadGraph{292} = [289 293 291];
RoadGraph{293} = [315 294 292];
RoadGraph{294} = [282 293 295];
RoadGraph{295} = [294 296 281 313];
RoadGraph{296} = [312 297 295];
RoadGraph{297} = [277 296];
RoadGraph{298} = [297 299];
RoadGraph{299} = [298 310 300];
RoadGraph{300} = [309 299 301];
RoadGraph{301} = [275 302 300 308];
RoadGraph{302} = [274 301 303];
RoadGraph{303} = [306 304 302];
RoadGraph{304} = [272 303 305];
RoadGraph{305} = [304 306 348];
RoadGraph{306} = [347 305 307];
RoadGraph{307} = [339 306 327];
RoadGraph{308} = [301 327 309];
RoadGraph{309} = [300 308 326];
RoadGraph{310} = [299 309 311 325];
RoadGraph{311} = [298 312 310];
RoadGraph{312} = [324 311 313];
RoadGraph{313} = [295 312 323 314];
RoadGraph{314} = [294 313 315];
RoadGraph{315} = [321 314 316];
RoadGraph{316} = [292 317 315];
RoadGraph{317} = [319 316 318];
RoadGraph{318} = [290 319 317];
RoadGraph{319} = [318 328 320];
RoadGraph{320} = [319 328 316 321];
RoadGraph{321} = [322 320];
RoadGraph{322} = [314 323 321];
RoadGraph{323} = [313 322 324];
RoadGraph{324} = [331 325 323];
RoadGraph{325} = [310 335 326 324];
RoadGraph{326} = [309 337 327 325];
RoadGraph{327} = [308 307 338 326];
RoadGraph{328} = [320 319 329];
RoadGraph{329} = [322 328 330];
RoadGraph{330} = [329 332 331];
RoadGraph{331} = [332 333 330];
RoadGraph{332} = [334 330 340];
RoadGraph{333} = [331 335];
RoadGraph{334} = [333 336 332];
RoadGraph{335} = [325 336 337 333];
RoadGraph{336} = [335 341];
RoadGraph{337} = [326 335 338 342];
RoadGraph{338} = [327 339 337 343];
RoadGraph{339} = [307 338 344];
RoadGraph{340} = [332 351 341];
RoadGraph{341} = [336 342 340];
RoadGraph{342} = [337 351 341 343];
RoadGraph{343} = [338 344 342 352];
RoadGraph{344} = [339 343 345 354];
RoadGraph{345} = [346 344 355];
RoadGraph{346} = [348 350 345 356];
RoadGraph{347} = [345 348];
RoadGraph{348} = [305 349 347 346];
RoadGraph{349} = [270 348 350];
RoadGraph{350} = [349 346 357];
RoadGraph{351} = [342 340 352];
RoadGraph{352} = [343 351 353];
RoadGraph{353} = [354 352];
RoadGraph{354} = [344 355 353];
RoadGraph{355} = [354 356];
RoadGraph{356} = [355 357 346];
RoadGraph{357} = [350 356];



% define node locations
NodesLoc=zeros(N,2);
NodesLatLon = zeros(N,3);



NodesLatLon(1,1:2) = [40.701394, -74.012270]';
NodesLatLon(2,1:2) = [40.702435, -74.012849]';
NodesLatLon(3,1:2) = [40.703265, -74.007881]';
NodesLatLon(4,1:2) = [40.704151, -74.008944]';
NodesLatLon(5,1:2) = [40.709462, -74.001712]';
NodesLatLon(6,1:2) = [40.708039, -73.999856]';
NodesLatLon(7,1:2) = [40.704615, -74.014287]';
NodesLatLon(8,1:2) = [40.704859, -74.016904]';
NodesLatLon(9,1:2) = [40.707706, -74.015553]';
NodesLatLon(10,1:2) = [40.704883, -74.015295]';
NodesLatLon(11,1:2) = [40.706982, -74.013546]';
NodesLatLon(12,1:2) = [40.707417, -74.011975]';
NodesLatLon(13,1:2) = [40.712081, -74.010135]';
NodesLatLon(14,1:2) = [40.711386, -74.008627]';
NodesLatLon(15,1:2) = [40.712728, -74.009695]';
NodesLatLon(16,1:2) = [40.711996, -74.008101]';
NodesLatLon(17,1:2) = [40.714427, -74.013412]';
NodesLatLon(18,1:2) = [40.712577, -74.004699]';
NodesLatLon(19,1:2) = [40.713130, -74.004109]';
NodesLatLon(20,1:2) = [40.714187, -74.006308]';
NodesLatLon(21,1:2) = [40.714943, -74.008068]';
NodesLatLon(22,1:2) = [40.715513, -74.009202]';
NodesLatLon(23,1:2) = [40.717123, -74.012802]';
NodesLatLon(24,1:2) = [40.709707, -73.991774]';
NodesLatLon(25,1:2) = [40.713838, -73.992654]';
NodesLatLon(26,1:2) = [40.713708, -73.994241]';
NodesLatLon(27,1:2) = [40.713503, -73.998500]';
NodesLatLon(28,1:2) = [40.715166, -74.002175]';
NodesLatLon(29,1:2) = [40.715553, -74.003043]';
NodesLatLon(30,1:2) = [40.716268, -74.004574]';
NodesLatLon(31,1:2) = [40.717061, -74.006334]';
NodesLatLon(32,1:2) = [40.718260, -74.007151]';
NodesLatLon(33,1:2) = [40.718377, -74.005263]';
NodesLatLon(34,1:2) = [40.719498, -74.004370]';
NodesLatLon(35,1:2) = [40.719846, -74.005212]';
NodesLatLon(36,1:2) = [40.720485, -74.006717]';
NodesLatLon(37,1:2) = [40.720581, -74.008529]';
NodesLatLon(38,1:2) = [40.722322, -74.011719]';
NodesLatLon(39,1:2) = [40.721993, -74.008274]';
NodesLatLon(40,1:2) = [40.725764, -74.011030]';
NodesLatLon(41,1:2) = [40.723714, -74.007968]';
NodesLatLon(42,1:2) = [40.722593, -74.006334]';
NodesLatLon(43,1:2) = [40.721800, -74.006436]';
NodesLatLon(44,1:2) = [40.721761, -74.005263]';
NodesLatLon(45,1:2) = [40.720485, -74.003451]';
NodesLatLon(46,1:2) = [40.718667, -74.002583]';
NodesLatLon(47,1:2) = [40.719344, -74.001945]';
NodesLatLon(48,1:2) = [40.717932, -74.001001]';
NodesLatLon(49,1:2) = [40.718493, -74.000542]';
NodesLatLon(50,1:2) = [40.717603, -74.000261]';
NodesLatLon(51,1:2) = [40.718028, -74.000006]';
NodesLatLon(52,1:2) = [40.717197, -73.998858]';
NodesLatLon(53,1:2) = [40.716249, -73.996153]';
NodesLatLon(54,1:2) = [40.715862, -73.994902]';
NodesLatLon(55,1:2) = [40.718122, -73.993898]';
NodesLatLon(56,1:2) = [40.717368, -73.991295]';
NodesLatLon(57,1:2) = [40.716740, -73.989164]';
NodesLatLon(58,1:2) = [40.714022, -73.990262]';
NodesLatLon(59,1:2) = [40.710734, -73.984724]';
NodesLatLon(60,1:2) = [40.715221, -73.984240]';
NodesLatLon(61,1:2) = [40.714457, -73.984750]';
NodesLatLon(62,1:2) = [40.713277, -73.977580]';
NodesLatLon(63,1:2) = [40.719072, -73.993376]';
NodesLatLon(64,1:2) = [40.719455, -73.994440]';
NodesLatLon(65,1:2) = [40.720467, -73.994025]';
NodesLatLon(66,1:2) = [40.719975, -73.992888]';
NodesLatLon(67,1:2) = [40.719223, -73.990362]';
NodesLatLon(68,1:2) = [40.718478, -73.988258]';
NodesLatLon(69,1:2) = [40.714949, -73.976606]';
NodesLatLon(70,1:2) = [40.718893, -73.975136]';
NodesLatLon(71,1:2) = [40.720731, -73.981305]';
NodesLatLon(72,1:2) = [40.722178, -73.986326]';
NodesLatLon(73,1:2) = [40.722959, -73.988601]';
NodesLatLon(74,1:2) = [40.723674, -73.991090]';
NodesLatLon(75,1:2) = [40.724162, -73.992592]';
NodesLatLon(76,1:2) = [40.720820, -73.997656]';
NodesLatLon(77,1:2) = [40.721463, -73.997366]';
NodesLatLon(78,1:2) = [40.722274, -73.997060]';
NodesLatLon(79,1:2) = [40.721602, -73.997800]';
NodesLatLon(80,1:2) = [40.721075, -73.998233]';
NodesLatLon(81,1:2) = [40.721786, -73.999873]';
NodesLatLon(82,1:2) = [40.723364, -74.002966]';
NodesLatLon(83,1:2) = [40.723625, -74.004803]';
NodesLatLon(84,1:2) = [40.724185, -74.004599]';
NodesLatLon(85,1:2) = [40.723789, -74.006117]';
NodesLatLon(86,1:2) = [40.728614, -74.005339]';
NodesLatLon(87,1:2) = [40.728769, -74.007125]';
NodesLatLon(88,1:2) = [40.728372, -74.002876]';
NodesLatLon(89,1:2) = [40.725413, -73.996816]';
NodesLatLon(90,1:2) = [40.725102, -73.995308]';
NodesLatLon(91,1:2) = [40.729146, -74.010636]';
NodesLatLon(92,1:2) = [40.726780, -73.971921]';
NodesLatLon(93,1:2) = [40.728421, -73.975692]';
NodesLatLon(94,1:2) = [40.731320, -73.982603]';
NodesLatLon(95,1:2) = [40.732318, -73.984949]';
NodesLatLon(96,1:2) = [40.733289, -73.987222]';
NodesLatLon(97,1:2) = [40.728366, -73.990813]';
NodesLatLon(98,1:2) = [40.729569, -73.989893]';
NodesLatLon(99,1:2) = [40.729993, -73.991084]';
NodesLatLon(100,1:2) = [40.730568, -73.992455]';
NodesLatLon(101,1:2) = [40.732277, -73.996406]';
NodesLatLon(102,1:2) = [40.734402, -73.989883]';
NodesLatLon(103,1:2) = [40.734808, -73.990816]';
NodesLatLon(104,1:2) = [40.735967, -73.993673]';
NodesLatLon(105,1:2) = [40.733633, -73.999560]';
NodesLatLon(106,1:2) = [40.737368, -73.996911]';
NodesLatLon(107,1:2) = [40.736540, -74.001160]';
NodesLatLon(108,1:2) = [40.738509, -73.999681]';
NodesLatLon(109,1:2) = [40.736875, -74.005581]';
NodesLatLon(110,1:2) = [40.739364, -74.002838]';
NodesLatLon(111,1:2) = [40.739760, -74.002477]';
NodesLatLon(112,1:2) = [40.740936, -74.005418]';
NodesLatLon(113,1:2) = [40.737045, -74.009910]';
NodesLatLon(114,1:2) = [40.742329, -74.008708]';
NodesLatLon(115,1:2) = [40.731065, -73.973963]';
NodesLatLon(116,1:2) = [40.732708, -73.974478]';
NodesLatLon(117,1:2) = [40.734984, -73.979928]';
NodesLatLon(118,1:2) = [40.735423, -73.975079]';
NodesLatLon(119,1:2) = [40.736837, -73.978533]';
NodesLatLon(120,1:2) = [40.737894, -73.980915]';
NodesLatLon(121,1:2) = [40.738870, -73.983211]';
NodesLatLon(122,1:2) = [40.739488, -73.984777]';
NodesLatLon(123,1:2) = [40.740265, -73.986453]';
NodesLatLon(124,1:2) = [40.740897, -73.988015]';
NodesLatLon(125,1:2) = [40.741518, -73.989609]';
NodesLatLon(126,1:2) = [40.742921, -73.992810]';
NodesLatLon(127,1:2) = [40.744150, -73.995642]';
NodesLatLon(128,1:2) = [40.745253, -73.998472]';
NodesLatLon(129,1:2) = [40.746512, -74.001393]';
NodesLatLon(130,1:2) = [40.747673, -74.004162]';
NodesLatLon(131,1:2) = [40.748098, -74.007674]';
NodesLatLon(132,1:2) = [40.749541, -74.006635]';
NodesLatLon(133,1:2) = [40.750236, -74.008365]';
NodesLatLon(134,1:2) = [40.756983, -74.004837]';
NodesLatLon(135,1:2) = [40.755797, -74.002003]';
NodesLatLon(136,1:2) = [40.754580, -73.999152]';
NodesLatLon(137,1:2) = [40.753378, -73.996319]';
NodesLatLon(138,1:2) = [40.752216, -73.993504]';
NodesLatLon(139,1:2) = [40.750995, -73.990653]';
NodesLatLon(140,1:2) = [40.749787, -73.987761]';
NodesLatLon(141,1:2) = [40.748419, -73.984589]';
NodesLatLon(142,1:2) = [40.747741, -73.982951]';
NodesLatLon(143,1:2) = [40.747097, -73.981358]';
NodesLatLon(144,1:2) = [40.746361, -73.979704]';
NodesLatLon(145,1:2) = [40.745713, -73.978119]';
NodesLatLon(146,1:2) = [40.744754, -73.975898]';
NodesLatLon(147,1:2) = [40.743779, -73.973527]';
NodesLatLon(148,1:2) = [40.743210, -73.972153]';
NodesLatLon(149,1:2) = [40.748046, -73.968173]';
NodesLatLon(150,1:2) = [40.748843, -73.969825]';
NodesLatLon(151,1:2) = [40.749818, -73.972218]';
NodesLatLon(152,1:2) = [40.750769, -73.974460]';
NodesLatLon(153,1:2) = [40.751476, -73.976091]';
NodesLatLon(154,1:2) = [40.752143, -73.977743]';
NodesLatLon(155,1:2) = [40.752817, -73.979320]';
NodesLatLon(156,1:2) = [40.753500, -73.980908]';
NodesLatLon(157,1:2) = [40.754849, -73.984138]';
NodesLatLon(158,1:2) = [40.756028, -73.986959]';
NodesLatLon(159,1:2) = [40.757230, -73.989824]';
NodesLatLon(160,1:2) = [40.758425, -73.992667]';
NodesLatLon(161,1:2) = [40.757872, -73.993075]';
NodesLatLon(162,1:2) = [40.758360, -73.994384]';
NodesLatLon(163,1:2) = [40.758921, -73.993933]';
NodesLatLon(164,1:2) = [40.758913, -73.995993]';
NodesLatLon(165,1:2) = [40.759676, -73.995467]';
NodesLatLon(166,1:2) = [40.760871, -73.998364]';
NodesLatLon(167,1:2) = [40.762001, -74.001111]';
NodesLatLon(168,1:2) = [40.771509, -73.994211]';
NodesLatLon(169,1:2) = [40.770293, -73.991378]';
NodesLatLon(170,1:2) = [40.769076, -73.988599]';
NodesLatLon(171,1:2) = [40.767942, -73.985731]';
NodesLatLon(172,1:2) = [40.766712, -73.982916]';
NodesLatLon(173,1:2) = [40.765537, -73.980011]';
NodesLatLon(174,1:2) = [40.764362, -73.977232]';
NodesLatLon(175,1:2) = [40.762995, -73.973930]';
NodesLatLon(176,1:2) = [40.762298, -73.972414]';
NodesLatLon(177,1:2) = [40.761656, -73.970754]';
NodesLatLon(178,1:2) = [40.760945, -73.969148]';
NodesLatLon(179,1:2) = [40.760248, -73.967597]';
NodesLatLon(180,1:2) = [40.759291, -73.965287]';
NodesLatLon(181,1:2) = [40.758321, -73.962941]';
NodesLatLon(182,1:2) = [40.755779, -73.964818]';
NodesLatLon(183,1:2) = [40.752594, -73.967146]';
NodesLatLon(184,1:2) = [40.751856, -73.965016]';
NodesLatLon(185,1:2) = [40.754836, -73.962220]';
NodesLatLon(186,1:2) = [40.758731, -73.958701]';
NodesLatLon(187,1:2) = [40.760043, -73.961660]';
NodesLatLon(188,1:2) = [40.760822, -73.962833]';
NodesLatLon(189,1:2) = [40.760590, -73.964367]';
NodesLatLon(190,1:2) = [40.761246, -73.963862]';
NodesLatLon(191,1:2) = [40.761574, -73.966622]';
NodesLatLon(192,1:2) = [40.762202, -73.966117]';
NodesLatLon(193,1:2) = [40.762257, -73.968264]';
NodesLatLon(194,1:2) = [40.762886, -73.967651]';
NodesLatLon(195,1:2) = [40.762913, -73.969816]';
NodesLatLon(196,1:2) = [40.763555, -73.969347]';
NodesLatLon(197,1:2) = [40.763569, -73.971458]';
NodesLatLon(198,1:2) = [40.764225, -73.970971]';
NodesLatLon(199,1:2) = [40.764266, -73.972992]';
NodesLatLon(200,1:2) = [40.764908, -73.972559]';
NodesLatLon(201,1:2) = [40.765660, -73.976240]';
NodesLatLon(202,1:2) = [40.766876, -73.979072]';
NodesLatLon(203,1:2) = [40.768215, -73.981436]';
NodesLatLon(204,1:2) = [40.772985, -73.982086]';
NodesLatLon(205,1:2) = [40.774179, -73.984811]';
NodesLatLon(206,1:2) = [40.775363, -73.987737]';
NodesLatLon(207,1:2) = [40.774809, -73.984458]';
NodesLatLon(208,1:2) = [40.773930, -73.982238]';
NodesLatLon(209,1:2) = [40.773568, -73.981582]';
NodesLatLon(210,1:2) = [40.771791, -73.979186]';
NodesLatLon(211,1:2) = [40.772479, -73.978807]';
NodesLatLon(212,1:2) = [40.768008, -73.970307]';
NodesLatLon(213,1:2) = [40.768658, -73.969828]';
NodesLatLon(214,1:2) = [40.767321, -73.968667]';
NodesLatLon(215,1:2) = [40.767970, -73.968264]';
NodesLatLon(216,1:2) = [40.766652, -73.967053]';
NodesLatLon(217,1:2) = [40.767302, -73.966649]';
NodesLatLon(218,1:2) = [40.765964, -73.965438]';
NodesLatLon(219,1:2) = [40.766633, -73.964984]';
NodesLatLon(220,1:2) = [40.765334, -73.963824]';
NodesLatLon(221,1:2) = [40.765945, -73.963446]';
NodesLatLon(222,1:2) = [40.768811, -73.958476]';
NodesLatLon(223,1:2) = [40.767798, -73.956055]';
NodesLatLon(224,1:2) = [40.766537, -73.950909]';
NodesLatLon(225,1:2) = [40.763118, -73.962512]';
NodesLatLon(226,1:2) = [40.762143, -73.960116]';
NodesLatLon(227,1:2) = [40.760729, -73.956736]';
NodesLatLon(228,1:2) = [40.769747, -73.960595]';
NodesLatLon(229,1:2) = [40.770396, -73.962260]';
NodesLatLon(230,1:2) = [40.771141, -73.963849]';
NodesLatLon(231,1:2) = [40.771829, -73.965438]';
NodesLatLon(232,1:2) = [40.772460, -73.967053]';
NodesLatLon(233,1:2) = [40.776223, -73.975932]';
NodesLatLon(234,1:2) = [40.777407, -73.978782]';
NodesLatLon(235,1:2) = [40.778629, -73.981759]';
NodesLatLon(236,1:2) = [40.779776, -73.984508]';
NodesLatLon(237,1:2) = [40.785582, -73.984029]';
NodesLatLon(238,1:2) = [40.783767, -73.979917]';
NodesLatLon(239,1:2) = [40.783118, -73.978328]';
NodesLatLon(240,1:2) = [40.780654, -73.976462]';
NodesLatLon(241,1:2) = [40.779508, -73.973662]';
NodesLatLon(242,1:2) = [40.782048, -73.971694]';
NodesLatLon(243,1:2) = [40.783233, -73.974595]';
NodesLatLon(244,1:2) = [40.781953, -73.975528]';
NodesLatLon(245,1:2) = [40.784512, -73.977370]';
NodesLatLon(246,1:2) = [40.777006, -73.963748]';
NodesLatLon(247,1:2) = [40.776299, -73.962084]';
NodesLatLon(248,1:2) = [40.775669, -73.960494]';
NodesLatLon(249,1:2) = [40.774962, -73.958880]';
NodesLatLon(250,1:2) = [40.774313, -73.957392]';
NodesLatLon(251,1:2) = [40.773377, -73.955021]';
NodesLatLon(252,1:2) = [40.772326, -73.952700]';
NodesLatLon(253,1:2) = [40.770320, -73.947554]';
NodesLatLon(254,1:2) = [40.776910, -73.949396]';
NodesLatLon(255,1:2) = [40.777885, -73.951792]';
NodesLatLon(256,1:2) = [40.778840, -73.953961]';
NodesLatLon(257,1:2) = [40.779489, -73.955601]';
NodesLatLon(258,1:2) = [40.780177, -73.957240]';
NodesLatLon(259,1:2) = [40.780864, -73.958855]';
NodesLatLon(260,1:2) = [40.781514, -73.960394]';
NodesLatLon(261,1:2) = [40.780826, -73.960898]';
NodesLatLon(262,1:2) = [40.785314, -73.969373]';
NodesLatLon(263,1:2) = [40.786537, -73.972274]';
NodesLatLon(264,1:2) = [40.787664, -73.975049]';
NodesLatLon(265,1:2) = [40.788332, -73.976512]';
NodesLatLon(266,1:2) = [40.789975, -73.980321]';
NodesLatLon(267,1:2) = [40.796276, -73.976285]';
NodesLatLon(268,1:2) = [40.785066, -73.982591]';
NodesLatLon(269,1:2) = [40.795474, -73.975831]';
NodesLatLon(270,1:2) = [40.796105, -73.974923]';
NodesLatLon(271,1:2) = [40.793947, -73.972274]';
NodesLatLon(272,1:2) = [40.794653, -73.971770]';
NodesLatLon(273,1:2) = [40.794080, -73.970357]';
NodesLatLon(274,1:2) = [40.792839, -73.967557]';
NodesLatLon(275,1:2) = [40.791697, -73.964734]';
NodesLatLon(276,1:2) = [40.787915, -73.955766]';
NodesLatLon(277,1:2) = [40.788590, -73.955240]';
NodesLatLon(278,1:2) = [40.787280, -73.954179]';
NodesLatLon(279,1:2) = [40.787922, -73.953679]';
NodesLatLon(280,1:2) = [40.786544, -73.952609]';
NodesLatLon(281,1:2) = [40.787260, -73.952083]';
NodesLatLon(282,1:2) = [40.785882, -73.950941]';
NodesLatLon(283,1:2) = [40.785200, -73.949345]';
NodesLatLon(284,1:2) = [40.784316, -73.947151]';
NodesLatLon(285,1:2) = [40.783262, -73.944761]';
NodesLatLon(286,1:2) = [40.782979, -73.944048]';
NodesLatLon(287,1:2) = [40.783668, -73.943513]';
NodesLatLon(288,1:2) = [40.783958, -73.944226]';
NodesLatLon(289,1:2) = [40.785005, -73.946616]';
NodesLatLon(290,1:2) = [40.788717, -73.937871]';
NodesLatLon(291,1:2) = [40.789634, -73.940061]';
NodesLatLon(292,1:2) = [40.790618, -73.942409]';
NodesLatLon(293,1:2) = [40.791555, -73.944700]';
NodesLatLon(294,1:2) = [40.792257, -73.946272]';
NodesLatLon(295,1:2) = [40.792926, -73.947894]';
NodesLatLon(296,1:2) = [40.793619, -73.949497]';
NodesLatLon(297,1:2) = [40.794293, -73.951113]';
NodesLatLon(298,1:2) = [40.796954, -73.949200]';
NodesLatLon(299,1:2) = [40.798258, -73.952463]';
NodesLatLon(300,1:2) = [40.799467, -73.955278]';
NodesLatLon(301,1:2) = [40.800547, -73.958257]';
NodesLatLon(302,1:2) = [40.801826, -73.961089]';
NodesLatLon(303,1:2) = [40.803008, -73.963851]';
NodesLatLon(304,1:2) = [40.804178, -73.966643]';
NodesLatLon(305,1:2) = [40.815813, -73.958402]';
NodesLatLon(306,1:2) = [40.813325, -73.956240]';
NodesLatLon(307,1:2) = [40.810768, -73.952623]';
NodesLatLon(308,1:2) = [40.804440, -73.955415]';
NodesLatLon(309,1:2) = [40.803213, -73.952502]';
NodesLatLon(310,1:2) = [40.802010, -73.949620]';
NodesLatLon(311,1:2) = [40.800657, -73.946460]';
NodesLatLon(312,1:2) = [40.799998, -73.944869]';
NodesLatLon(313,1:2) = [40.799305, -73.943203]';
NodesLatLon(314,1:2) = [40.798657, -73.941626]';
NodesLatLon(315,1:2) = [40.797941, -73.939975]';
NodesLatLon(316,1:2) = [40.796975, -73.937814]';
NodesLatLon(317,1:2) = [40.796021, -73.935442]';
NodesLatLon(318,1:2) = [40.794327, -73.931419]';
NodesLatLon(319,1:2) = [40.801418, -73.930293]';
NodesLatLon(320,1:2) = [40.802770, -73.933626]';
NodesLatLon(321,1:2) = [40.803725, -73.935832]';
NodesLatLon(322,1:2) = [40.804407, -73.937438]';
NodesLatLon(323,1:2) = [40.805020, -73.938924]';
NodesLatLon(324,1:2) = [40.805713, -73.940651]';
NodesLatLon(325,1:2) = [40.807804, -73.945319]';
NodesLatLon(326,1:2) = [40.808963, -73.948291]';
NodesLatLon(327,1:2) = [40.810201, -73.951113]';
NodesLatLon(328,1:2) = [40.805111, -73.931989]';
NodesLatLon(329,1:2) = [40.807043, -73.933776]';
NodesLatLon(330,1:2) = [40.811508, -73.934796]';
NodesLatLon(331,1:2) = [40.812110, -73.935982]';
NodesLatLon(332,1:2) = [40.814064, -73.934811]';
NodesLatLon(333,1:2) = [40.812757, -73.937649]';
NodesLatLon(334,1:2) = [40.814712, -73.936192]';
NodesLatLon(335,1:2) = [40.814098, -73.940861]';
NodesLatLon(336,1:2) = [40.815416, -73.939930]';
NodesLatLon(337,1:2) = [40.815314, -73.943713]';
NodesLatLon(338,1:2) = [40.816552, -73.946535]';
NodesLatLon(339,1:2) = [40.817086, -73.948021]';
NodesLatLon(340,1:2) = [40.819710, -73.934526]';
NodesLatLon(341,1:2) = [40.820437, -73.936192]';
NodesLatLon(342,1:2) = [40.821618, -73.939060]';
NodesLatLon(343,1:2) = [40.822834, -73.941882]';
NodesLatLon(344,1:2) = [40.824123, -73.944858]';
NodesLatLon(345,1:2) = [40.825251, -73.947634]';
NodesLatLon(346,1:2) = [40.826435, -73.950475]';
NodesLatLon(347,1:2) = [40.818946, -73.952242]';
NodesLatLon(348,1:2) = [40.820102, -73.955094]';
NodesLatLon(349,1:2) = [40.821172, -73.957531]';
NodesLatLon(350,1:2) = [40.827094, -73.952015]';
NodesLatLon(351,1:2) = [40.828135, -73.934886]';
NodesLatLon(352,1:2) = [40.829176, -73.937271]';
NodesLatLon(353,1:2) = [40.830399, -73.940200]';
NodesLatLon(354,1:2) = [40.830914, -73.941374]';
NodesLatLon(355,1:2) = [40.831573, -73.943015]';
NodesLatLon(356,1:2) = [40.832786, -73.945855]';
NodesLatLon(357,1:2) = [40.834448, -73.949149]';


Coordinates(1,:)=NodesLatLon(:,1)';
Coordinates(2,:)=NodesLatLon(:,2)';

% link lengths
LinkLength = sparse(N,N);
NumLanes = 2; %sparse(N,N);
% define link capacities
%LinkNumVehicles = sparse(N,N);
%LinkCapacityLeft = sparse(N,N);
LinkTime = zeros(N,N);
LinkSpeed = sparse(N,N);

NominalCapacity = .2;
CAR_LENGTH = 6.5; % meters/car


% define link freeflow speed
LinkFreeFlow = sparse(N,N);

RefLocation = [40.689235 -74.044385];
tmpLocation = lla2flat(NodesLatLon, RefLocation, 28.8, 0);%28.8
NodesLoc(:,1) = tmpLocation(:,2); %x (longitude)
NodesLoc(:,2) = tmpLocation(:,1); %y (latitude)
%%% base of network

% make graph symmetric
if symmetricFlag == 1
    for i = 1:N
        for j = RoadGraph{i}
            if isempty(find(RoadGraph{j}==i,1))
                RoadGraph{j} = [RoadGraph{j} i];
            end
        end
    end
end

% find link lengths
for i = 1:N
    for j = RoadGraph{i}
        LinkLength(i,j) = norm(NodesLoc(j,:) - NodesLoc(i,:),1); %m
        LinkFreeFlow(i,j) = 8; %14; % m/s
        %NumLanes(i,j) = 1;
        LinkTime(i,j) = LinkLength(i,j)/LinkFreeFlow(i,j); %s
        %RoadCap(i,j) = NominalCapacity*NumLanes(i,j)*LinkFreeFlow(i,j)/CAR_LENGTH*3600;
        %RoadCap(i,j) = ceil(NominalCapacity*NumLanes(i,j)*LinkLength(i,j));
    end
end

% polyx = [1.8; 3.8; 3.6; 2.6; 2.75; 3; 3; 1; 0.9; -1; -1; 0.25; 1; 1.8];
% polyx = polyx*1000;
% polyy = [2.4; 5.15; 7; 8.1; 10; 12.5; 15; 17.5; 18.1; 18.1; 7; 4; 2.4; 2.4];
% polyy = polyy*1000;
%%
%RoadCap = RoadCap;

%%Graph Creation

G_road=digraph(LinkTime/(60)); %min
FreeFlowTime = G_road.Edges.Weight'; %min

load('NYC_Demands')
DemandS = zeros(357);
for iii=1:53000
    x = outdata(iii,2);
    y = outdata(iii,3);
    if x && y > 0
    DemandS(x,y) = DemandS(x,y) +1;
    end
end
ADJ=adjacency(G_road);
% figure(1)
% plot(graph(ADJ+ADJ'),'YData',NodesLoc(:,1)/1000,'XData',NodesLoc(:,2)/1000,'Marker','.','Linewidth',1,'EdgeColor','k')
% %pp=plot(G_road,'XData',NodesLoc(:,1),'YData',NodesLoc(:,2),'Linewidth',1)
% ylim([-2 6])
% axis equal
% xlabel('kilometers [km]')
% ylabel('kilometers [km]')
% grid on
% set(gca,'ticklabelinterpreter','Latex','fontsize',16)
%saveas(figure(1),strcat('Images/','NYC_networks.pdf'))
save('Graphs.mat', 'G_road','NodesLoc','LinkTime','FreeFlowTime','DemandS');

%% Pruned version
N_no = 80;
[idx,C] = kmeans(NodesLoc,N_no);
RoadCap = sparse(N_no,N_no);
NodesLoc2 = C;
% find link lengths
for i = 1:N_no
    oldNodes = find(idx==i);
    ConnectedNodes = ADJ(oldNodes,:);
    [a,b]=find(ConnectedNodes==1);
    IndNewNodes = idx(b);
    for k = 1:length(IndNewNodes)
        j = IndNewNodes(k);
        LinkLength1(i,j) = norm(NodesLoc2(j,:) - NodesLoc2(i,:)); %m
        LinkFreeFlow1(i,j) = 8; %14; % m/s
        LinkTime1(i,j) = LinkLength1(i,j)/LinkFreeFlow1(i,j); %s
        %RoadCap(i,j) = NominalCapacity*NumLanes(i,j)*LinkFreeFlow(i,j)/CAR_LENGTH*3600;
        RoadCap(i,j) = NominalCapacity*NumLanes*LinkLength(i,j);
    end
end

for i = 1:N_no
    LinkLength1(i,i) = 0; %m
    LinkFreeFlow1(i,i) = 0; %14; % m/s
    LinkTime1(i,i) = 0; %s
end
%%

G_road=digraph(LinkTime1/(60)); %min
FreeFlowTime = G_road.Edges.Weight'; %min
N_lanes = 2;
Capacity1 = FreeFlowTime*0 + NominalCapacity*N_lanes*8*3600/CAR_LENGTH;

Old=load('Graphs')
DemandS = zeros(N_no,N_no);
for ii = 1 : 357
    for jj = 1 : 357 % old nodes
        if Old.DemandS(ii,jj)~= 0
            ii1 = idx(ii);
            jj1 = idx(jj);
            if ii1 ~=jj1
            DemandS(ii1,jj1) = DemandS(ii1,jj1) + Old.DemandS(ii,jj);
            else
            NodAdj = find(LinkTime1(ii1,:)~=0);    
            DemandS(ii1,NodAdj(1)) = DemandS(ii1,NodAdj(1)) + Old.DemandS(ii,jj);    
            end
        end
    end
end
Capacity = Capacity1;

% DemandS = zeros(357);
% for iii=1:53000
%     x = outdata(iii,2);
%     y = outdata(iii,3);
%     if x && y > 0
%     DemandS(x,y) = DemandS(x,y) +1;
%     end
% end
% ADJ=adjacency(G_road);
% figure(1)
% plot(graph(ADJ+ADJ'),'YData',NodesLoc(:,1)/1000,'XData',NodesLoc(:,2)/1000,'Marker','.','Linewidth',1,'EdgeColor','k')
% %pp=plot(G_road,'XData',NodesLoc(:,1),'YData',NodesLoc(:,2),'Linewidth',1)
% ylim([-2 6])
% axis equal
% xlabel('kilometers [km]')
% ylabel('kilometers [km]')
% grid on
% set(gca,'ticklabelinterpreter','Latex','fontsize',16)
%saveas(figure(1),strcat('Images/','NYC_networks.pdf'))
NodesLoc = NodesLoc2;
save(strcat('Graphs',num2str(N_no),'.mat'), 'G_road','NodesLoc','FreeFlowTime','DemandS','Capacity');
