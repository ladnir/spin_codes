import SpinCodes.UniformMixer.Map22PacketProfileWitness
namespace Spin.UniformMixer
open Spin.Structured
set_option maxRecDepth 10000
set_option maxHeartbeats 0
def profileZeroWitness0 : ProfileZeroWitness :=
  ⟨0, 15, 18364758544493064720, 18364758544493064720, [1,2,4,8,16,32,64,128,256,512,1024]⟩
theorem profileZeroWitness0_checked : profileZeroWitness0.check = true := by decide
def profileZeroWitness1 : ProfileZeroWitness :=
  ⟨0, 51, 18355225778678027280, 18355225778678027280, [1,2,8,4,16,64,32,128,256,1024,512]⟩
theorem profileZeroWitness1_checked : profileZeroWitness1.check = true := by decide
def profileZeroWitness2 : ProfileZeroWitness :=
  ⟨0, 60, 18345655775496524850, 18355151303945097300, [1,2,13,4,16,98,32,128,256,1552,512]⟩
theorem profileZeroWitness2_checked : profileZeroWitness2.check = true := by decide
def profileZeroWitness3 : ProfileZeroWitness :=
  ⟨0, 85, 18282904370552071200, 18152560703022715200, [1,8,2,4,16,64,256,1024,32,128,512]⟩
theorem profileZeroWitness3_checked : profileZeroWitness3.check = true := by decide
def profileZeroWitness4 : ProfileZeroWitness :=
  ⟨0, 90, 18278119368961319985, 18152494474626995460, [1,13,2,4,16,98,256,1552,32,128,512]⟩
theorem profileZeroWitness4_checked : profileZeroWitness4.check = true := by decide
def profileZeroWitness5 : ProfileZeroWitness :=
  ⟨0, 102, 18210566405327185185, 18135671946718105860, [1,11,2,4,16,64,292,1168,32,128,512]⟩
theorem profileZeroWitness5_checked : profileZeroWitness5.check = true := by decide
def profileZeroWitness6 : ProfileZeroWitness :=
  ⟨0, 105, 18206344345100051760, 18135606233718461760, [1,14,2,4,16,98,292,1664,32,128,512]⟩
theorem profileZeroWitness6_checked : profileZeroWitness6.check = true := by decide
def profileZeroWitness7 : ProfileZeroWitness :=
  ⟨0, 150, 17129718987181028385, 13829104802042368260, [1,15,2,4,16,96,288,1680,32,128,512]⟩
theorem profileZeroWitness7_checked : profileZeroWitness7.check = true := by decide
def profileZeroWitness8 : ProfileZeroWitness :=
  ⟨0, 153, 17125496926953894960, 13829039089042724160, [1,10,2,4,16,66,288,1152,32,128,512]⟩
theorem profileZeroWitness8_checked : profileZeroWitness8.check = true := by decide
def profileZeroWitness9 : ProfileZeroWitness :=
  ⟨0, 165, 17057943963319760160, 13812216561133834560, [1,12,2,4,16,96,260,1536,32,128,512]⟩
theorem profileZeroWitness9_checked : profileZeroWitness9.check = true := by decide
def profileZeroWitness10 : ProfileZeroWitness :=
  ⟨0, 170, 17053158961729008945, 13812150332738114820, [1,9,2,4,16,66,260,1040,32,128,512]⟩
theorem profileZeroWitness10_checked : profileZeroWitness10.check = true := by decide
def profileZeroWitness11 : ProfileZeroWitness :=
  ⟨0, 195, 15905304964213405200, 13474449681378858000, [1,2,12,4,16,96,32,128,260,1536,512]⟩
theorem profileZeroWitness11_checked : profileZeroWitness11.check = true := by decide
def profileZeroWitness12 : ProfileZeroWitness :=
  ⟨0, 204, 15895734961031902770, 13474375206645928020, [1,2,9,4,16,66,32,128,260,1040,512]⟩
theorem profileZeroWitness12_checked : profileZeroWitness12.check = true := by decide
def profileZeroWitness13 : ProfileZeroWitness :=
  ⟨0, 240, 13445776909200815700, 13445776909200815700, [1,2,4,9,16,32,66,128,260,512,1040]⟩
theorem profileZeroWitness13_checked : profileZeroWitness13.check = true := by decide
def profileZeroWitness14 : ProfileZeroWitness :=
  ⟨1, 255, 18364758544493064720, 18364758544493064720, [1,2,4,8,16,32,64,128,256,512,1024]⟩
theorem profileZeroWitness14_checked : profileZeroWitness14.check = true := by decide
def profileZeroWitness15 : ProfileZeroWitness :=
  ⟨0, 771, 18336085189344794640, 18355113192123290640, [1,2,8,16,4,64,128,32,1024,256,512]⟩
theorem profileZeroWitness15_checked : profileZeroWitness15.check = true := by decide
def profileZeroWitness16 : ProfileZeroWitness :=
  ⟨0, 780, 18326515186163292210, 18355113190986879060, [1,2,13,16,4,98,128,32,1552,256,512]⟩
theorem profileZeroWitness16_checked : profileZeroWitness16.check = true := by decide
def profileZeroWitness17 : ProfileZeroWitness :=
  ⟨0, 816, 18316944598875150420, 18355113189832711320, [1,2,8,21,4,64,162,32,1288,256,512]⟩
theorem profileZeroWitness17_checked : profileZeroWitness17.check = true := by decide
def profileZeroWitness18 : ProfileZeroWitness :=
  ⟨0, 960, 18307374595693647990, 18355113188696299740, [1,2,13,21,4,98,162,32,1821,256,512]⟩
theorem profileZeroWitness18_checked : profileZeroWitness18.check = true := by decide
def profileZeroWitness19 : ProfileZeroWitness :=
  ⟨0, 1285, 18263763781218838560, 18152448116467978560, [1,8,2,16,4,64,1024,256,128,32,512]⟩
theorem profileZeroWitness19_checked : profileZeroWitness19.check = true := by decide
def profileZeroWitness20 : ProfileZeroWitness :=
  ⟨0, 1290, 18258978779628087345, 18152448115457398020, [1,13,2,16,4,98,1552,256,128,32,512]⟩
theorem profileZeroWitness20_checked : profileZeroWitness20.check = true := by decide
def profileZeroWitness21 : ProfileZeroWitness :=
  ⟨0, 1360, 18244623190749194340, 18152448114177399240, [1,8,2,21,4,64,1288,256,162,32,512]⟩
theorem profileZeroWitness21_checked : profileZeroWitness21.check = true := by decide
def profileZeroWitness22 : ProfileZeroWitness :=
  ⟨0, 1440, 18239838189158443125, 18152448113166818700, [1,13,2,21,4,98,1821,256,162,32,512]⟩
theorem profileZeroWitness22_checked : profileZeroWitness22.check = true := by decide
def profileZeroWitness23 : ProfileZeroWitness :=
  ⟨0, 1542, 18191425815993952545, 18135559360163369220, [1,11,2,16,4,64,1168,292,128,32,512]⟩
theorem profileZeroWitness23_checked : profileZeroWitness23.check = true := by decide
def profileZeroWitness24 : ProfileZeroWitness :=
  ⟨0, 1545, 18187203755766819120, 18135559359160653120, [1,14,2,16,4,98,1664,292,128,32,512]⟩
theorem profileZeroWitness24_checked : profileZeroWitness24.check = true := by decide
def profileZeroWitness25 : ProfileZeroWitness :=
  ⟨0, 1632, 18172285225524308325, 18135559357872789900, [1,11,2,21,4,64,1471,292,162,32,512]⟩
theorem profileZeroWitness25_checked : profileZeroWitness25.check = true := by decide
def profileZeroWitness26 : ProfileZeroWitness :=
  ⟨0, 1680, 18168063165297174900, 18135559356870073800, [1,14,2,21,4,98,1962,292,162,32,512]⟩
theorem profileZeroWitness26_checked : profileZeroWitness26.check = true := by decide
def profileZeroWitness27 : ProfileZeroWitness :=
  ⟨0, 2310, 17110578397847795745, 13828926503490703620, [1,15,2,16,4,96,1680,288,128,32,512]⟩
theorem profileZeroWitness27_checked : profileZeroWitness27.check = true := by decide
def profileZeroWitness28 : ProfileZeroWitness :=
  ⟨0, 2313, 17106356337620662320, 13828926502487987520, [1,10,2,16,4,66,1152,288,128,32,512]⟩
theorem profileZeroWitness28_checked : profileZeroWitness28.check = true := by decide
def profileZeroWitness29 : ProfileZeroWitness :=
  ⟨0, 2400, 17091437807378151525, 13828926501200124300, [1,15,2,21,4,96,1983,288,162,32,512]⟩
theorem profileZeroWitness29_checked : profileZeroWitness29.check = true := by decide
def profileZeroWitness30 : ProfileZeroWitness :=
  ⟨0, 2448, 17087215747151018100, 13828926500197408200, [1,10,2,21,4,66,1450,288,162,32,512]⟩
theorem profileZeroWitness30_checked : profileZeroWitness30.check = true := by decide
def profileZeroWitness31 : ProfileZeroWitness :=
  ⟨0, 2565, 17038803373986527520, 13812037747193958720, [1,12,2,16,4,96,1536,260,128,32,512]⟩
theorem profileZeroWitness31_checked : profileZeroWitness31.check = true := by decide
def profileZeroWitness32 : ProfileZeroWitness :=
  ⟨0, 2570, 17034018372395776305, 13812037746183378180, [1,9,2,16,4,66,1040,260,128,32,512]⟩
theorem profileZeroWitness32_checked : profileZeroWitness32.check = true := by decide
def profileZeroWitness33 : ProfileZeroWitness :=
  ⟨0, 2640, 17019662783516883300, 13812037744903379400, [1,12,2,21,4,96,1800,260,162,32,512]⟩
theorem profileZeroWitness33_checked : profileZeroWitness33.check = true := by decide
def profileZeroWitness34 : ProfileZeroWitness :=
  ⟨0, 2720, 17014877781926132085, 13812037743892798860, [1,9,2,21,4,66,1309,260,162,32,512]⟩
theorem profileZeroWitness34_checked : profileZeroWitness34.check = true := by decide
def profileZeroWitness35 : ProfileZeroWitness :=
  ⟨0, 3075, 15886164374880172560, 13474262621227602960, [1,2,12,16,4,96,128,32,1536,260,512]⟩
theorem profileZeroWitness35_checked : profileZeroWitness35.check = true := by decide
def profileZeroWitness36 : ProfileZeroWitness :=
  ⟨0, 3084, 15876594371698670130, 13474262620091191380, [1,2,9,16,4,66,128,32,1040,260,512]⟩
theorem profileZeroWitness36_checked : profileZeroWitness36.check = true := by decide
def profileZeroWitness37 : ProfileZeroWitness :=
  ⟨0, 3120, 15867023784410528340, 13474262618937023640, [1,2,12,21,4,96,162,32,1800,260,512]⟩
theorem profileZeroWitness37_checked : profileZeroWitness37.check = true := by decide
def profileZeroWitness38 : ProfileZeroWitness :=
  ⟨0, 3264, 15857453781229025910, 13474262617800612060, [1,2,9,21,4,66,162,32,1309,260,512]⟩
theorem profileZeroWitness38_checked : profileZeroWitness38.check = true := by decide
def profileZeroWitness39 : ProfileZeroWitness :=
  ⟨0, 3840, 8526495043095935640, 8526495043095935640, [1,2,4,8,17,32,64,130,256,516,1032]⟩
theorem profileZeroWitness39_checked : profileZeroWitness39.check = true := by decide
def profileZeroWitness40 : ProfileZeroWitness :=
  ⟨1, 3855, 18364683486789906960, 18364683486789906960, [1,2,4,16,8,32,128,64,512,256,1024]⟩
theorem profileZeroWitness40_checked : profileZeroWitness40.check = true := by decide
def profileZeroWitness41 : ProfileZeroWitness :=
  ⟨1, 4080, 18364608427941459540, 18364683484499327640, [1,2,4,25,8,32,194,64,772,256,1024]⟩
theorem profileZeroWitness41_checked : profileZeroWitness41.check = true := by decide
def profileZeroWitness42 : ProfileZeroWitness :=
  ⟨0, 4369, 18119083436115347520, 18119083436115347520, [1,8,16,2,4,1024,64,256,128,512,32]⟩
theorem profileZeroWitness42_checked : profileZeroWitness42.check = true := by decide
def profileZeroWitness43 : ProfileZeroWitness :=
  ⟨0, 4386, 18114298434524596305, 18119083435104766980, [1,13,16,2,4,1552,98,256,128,512,32]⟩
theorem profileZeroWitness43_checked : profileZeroWitness43.check = true := by decide
def profileZeroWitness44 : ProfileZeroWitness :=
  ⟨0, 4420, 18109513140880525410, 18119083433842524360, [1,8,21,2,4,1288,64,256,162,512,32]⟩
theorem profileZeroWitness44_checked : profileZeroWitness44.check = true := by decide
def profileZeroWitness45 : ProfileZeroWitness :=
  ⟨0, 4488, 18104728139289774195, 18119083432831943820, [1,13,21,2,4,1821,98,256,162,512,32]⟩
theorem profileZeroWitness45_checked : profileZeroWitness45.check = true := by decide
def profileZeroWitness46 : ProfileZeroWitness :=
  ⟨0, 4626, 18046745470890461505, 18119017207719627780, [1,11,16,2,4,1168,64,292,128,512,32]⟩
theorem profileZeroWitness46_checked : profileZeroWitness46.check = true := by decide
def profileZeroWitness47 : ProfileZeroWitness :=
  ⟨0, 4641, 18042523410663328080, 18119017206709078080, [1,14,16,2,4,1664,98,292,128,512,32]⟩
theorem profileZeroWitness47_checked : profileZeroWitness47.check = true := by decide
def profileZeroWitness48 : ProfileZeroWitness :=
  ⟨0, 4680, 18037175175655639395, 18119017205446804620, [1,11,21,2,4,1471,64,292,162,512,32]⟩
theorem profileZeroWitness48_checked : profileZeroWitness48.check = true := by decide
def profileZeroWitness49 : ProfileZeroWitness :=
  ⟨0, 4740, 18032953115428505970, 18119017204436254920, [1,14,21,2,4,1962,98,292,162,512,32]⟩
theorem profileZeroWitness49_checked : profileZeroWitness49.check = true := by decide
def profileZeroWitness50 : ProfileZeroWitness :=
  ⟨0, 5140, 17974403090506566210, 18118934486649487560, [1,8,19,2,4,1096,64,256,128,548,32]⟩
theorem profileZeroWitness50_checked : profileZeroWitness50.check = true := by decide
def profileZeroWitness51 : ProfileZeroWitness :=
  ⟨0, 5160, 17969618088915814995, 18118934485638907020, [1,13,19,2,4,1663,98,256,128,548,32]⟩
theorem profileZeroWitness51_checked : profileZeroWitness51.check = true := by decide
def profileZeroWitness52 : ProfileZeroWitness :=
  ⟨0, 5185, 17965958712358193760, 18118934484376733760, [1,8,22,2,4,1344,64,256,162,548,32]⟩
theorem profileZeroWitness52_checked : profileZeroWitness52.check = true := by decide
def profileZeroWitness53 : ProfileZeroWitness :=
  ⟨0, 5250, 17961173710767442545, 18118934483366153220, [1,13,22,2,4,1906,98,256,162,548,32]⟩
theorem profileZeroWitness53_checked : profileZeroWitness53.check = true := by decide
def profileZeroWitness54 : ProfileZeroWitness :=
  ⟨0, 6168, 17902065125281680195, 18118868258253767820, [1,11,19,2,4,1243,64,292,128,548,32]⟩
theorem profileZeroWitness54_checked : profileZeroWitness54.check = true := by decide
def profileZeroWitness55 : ProfileZeroWitness :=
  ⟨0, 6180, 17897843065054546770, 18118868257243218120, [1,14,19,2,4,1772,98,292,128,548,32]⟩
theorem profileZeroWitness55_checked : profileZeroWitness55.check = true := by decide
def profileZeroWitness56 : ProfileZeroWitness :=
  ⟨0, 6210, 17893620747133307745, 18118868255981014020, [1,11,22,2,4,1524,64,292,162,548,32]⟩
theorem profileZeroWitness56_checked : profileZeroWitness56.check = true := by decide
def profileZeroWitness57 : ProfileZeroWitness :=
  ⟨0, 6273, 17889398686906174320, 18118868254970464320, [1,14,22,2,4,1990,98,292,162,548,32]⟩
theorem profileZeroWitness57_checked : profileZeroWitness57.check = true := by decide
def profileZeroWitness58 : ProfileZeroWitness :=
  ⟨0, 8466, 16965898052744304705, 13778739295237016580, [1,15,16,2,4,1680,96,288,128,512,32]⟩
theorem profileZeroWitness58_checked : profileZeroWitness58.check = true := by decide
def profileZeroWitness59 : ProfileZeroWitness :=
  ⟨0, 8481, 16961675992517171280, 13778739294226466880, [1,10,16,2,4,1152,66,288,128,512,32]⟩
theorem profileZeroWitness59_checked : profileZeroWitness59.check = true := by decide
def profileZeroWitness60 : ProfileZeroWitness :=
  ⟨0, 8520, 16956327757509482595, 13778739292964193420, [1,15,21,2,4,1983,96,288,162,512,32]⟩
theorem profileZeroWitness60_checked : profileZeroWitness60.check = true := by decide
def profileZeroWitness61 : ProfileZeroWitness :=
  ⟨0, 8580, 16952105697282349170, 13778739291953643720, [1,10,21,2,4,1450,66,288,162,512,32]⟩
theorem profileZeroWitness61_checked : profileZeroWitness61.check = true := by decide
def profileZeroWitness62 : ProfileZeroWitness :=
  ⟨0, 8721, 16894123028883036480, 13778673066841327680, [1,12,16,2,4,1536,96,260,128,512,32]⟩
theorem profileZeroWitness62_checked : profileZeroWitness62.check = true := by decide
def profileZeroWitness63 : ProfileZeroWitness :=
  ⟨0, 8738, 16889338027292285265, 13778673065830747140, [1,9,16,2,4,1040,66,260,128,512,32]⟩
theorem profileZeroWitness63_checked : profileZeroWitness63.check = true := by decide
def profileZeroWitness64 : ProfileZeroWitness :=
  ⟨0, 8772, 16884552733648214370, 13778673064568504520, [1,12,21,2,4,1800,96,260,162,512,32]⟩
theorem profileZeroWitness64_checked : profileZeroWitness64.check = true := by decide
def profileZeroWitness65 : ProfileZeroWitness :=
  ⟨0, 8840, 16879767732057463155, 13778673063557923980, [1,9,21,2,4,1309,66,260,162,512,32]⟩
theorem profileZeroWitness65_checked : profileZeroWitness65.check = true := by decide
def profileZeroWitness66 : ProfileZeroWitness :=
  ⟨0, 9240, 16821217707135523395, 13778590345771156620, [1,15,19,2,4,1791,96,288,128,548,32]⟩
theorem profileZeroWitness66_checked : profileZeroWitness66.check = true := by decide
def profileZeroWitness67 : ProfileZeroWitness :=
  ⟨0, 9252, 16816995646908389970, 13778590344760606920, [1,10,19,2,4,1224,66,288,128,548,32]⟩
theorem profileZeroWitness67_checked : profileZeroWitness67.check = true := by decide
def profileZeroWitness68 : ProfileZeroWitness :=
  ⟨0, 9282, 16812773328987150945, 13778590343498402820, [1,15,22,2,4,2000,96,288,162,548,32]⟩
theorem profileZeroWitness68_checked : profileZeroWitness68.check = true := by decide
def profileZeroWitness69 : ProfileZeroWitness :=
  ⟨0, 9345, 16808551268760017520, 13778590342487853120, [1,10,22,2,4,1506,66,288,162,548,32]⟩
theorem profileZeroWitness69_checked : profileZeroWitness69.check = true := by decide
def profileZeroWitness70 : ProfileZeroWitness :=
  ⟨0, 10260, 16749442683274255170, 13778524117375467720, [1,12,19,2,4,1644,96,260,128,548,32]⟩
theorem profileZeroWitness70_checked : profileZeroWitness70.check = true := by decide
def profileZeroWitness71 : ProfileZeroWitness :=
  ⟨0, 10280, 16744657681683503955, 13778524116364887180, [1,9,19,2,4,1115,66,260,128,548,32]⟩
theorem profileZeroWitness71_checked : profileZeroWitness71.check = true := by decide
def profileZeroWitness72 : ProfileZeroWitness :=
  ⟨0, 10305, 16740998305125882720, 13778524115102713920, [1,12,22,2,4,1892,96,260,162,548,32]⟩
theorem profileZeroWitness72_checked : profileZeroWitness72.check = true := by decide
def profileZeroWitness73 : ProfileZeroWitness :=
  ⟨0, 10370, 16736213303535131505, 13778524114092133380, [1,9,22,2,4,1366,66,260,162,548,32]⟩
theorem profileZeroWitness73_checked : profileZeroWitness73.check = true := by decide
def profileZeroWitness74 : ProfileZeroWitness :=
  ⟨0, 12291, 13436094029115874320, 8517149923829371920, [1,2,8,20,4,64,160,32,1280,256,516]⟩
theorem profileZeroWitness74_checked : profileZeroWitness74.check = true := by decide
def profileZeroWitness75 : ProfileZeroWitness :=
  ⟨0, 12300, 13426524025934371890, 8517149922692960340, [1,2,13,20,4,98,160,32,1808,256,516]⟩
theorem profileZeroWitness75_checked : profileZeroWitness75.check = true := by decide
def profileZeroWitness76 : ProfileZeroWitness :=
  ⟨0, 12336, 13416953438646230100, 8517149921538792600, [1,2,8,17,4,64,130,32,1032,256,516]⟩
theorem profileZeroWitness76_checked : profileZeroWitness76.check = true := by decide
def profileZeroWitness77 : ProfileZeroWitness :=
  ⟨0, 12480, 13407383435464727670, 8517149920402381020, [1,2,13,17,4,98,130,32,1565,256,516]⟩
theorem profileZeroWitness77_checked : profileZeroWitness77.check = true := by decide
def profileZeroWitness78 : ProfileZeroWitness :=
  ⟨0, 13056, 8516962277280898200, 8516962277280898200, [1,2,8,4,17,64,32,130,256,1032,516]⟩
theorem profileZeroWitness78_checked : profileZeroWitness78.check = true := by decide
def profileZeroWitness79 : ProfileZeroWitness :=
  ⟨1, 13107, 18355113192123290640, 18336085189344794640, [1,2,16,4,8,128,32,64,512,1024,256]⟩
theorem profileZeroWitness79_checked : profileZeroWitness79.check = true := by decide
def profileZeroWitness80 : ProfileZeroWitness :=
  ⟨1, 13260, 18355075662699066930, 18336085187071971480, [1,2,25,4,8,194,32,64,772,1024,256]⟩
theorem profileZeroWitness80_checked : profileZeroWitness80.check = true := by decide
def profileZeroWitness81 : ProfileZeroWitness :=
  ⟨0, 15360, 8507392274099395770, 8516887802547968220, [1,2,13,4,17,98,32,130,256,1565,516]⟩
theorem profileZeroWitness81_checked : profileZeroWitness81.check = true := by decide
def profileZeroWitness82 : ProfileZeroWitness :=
  ⟨1, 15420, 18345542896888468530, 18335936239878934680, [1,2,21,4,8,162,32,64,512,1288,256]⟩
theorem profileZeroWitness82_checked : profileZeroWitness82.check = true := by decide
def profileZeroWitness83 : ProfileZeroWitness :=
  ⟨1, 15555, 18345505659522020880, 18335936237606180880, [1,2,28,4,8,224,32,64,772,1288,256]⟩
theorem profileZeroWitness83_checked : profileZeroWitness83.check = true := by decide
def profileZeroWitness84 : ProfileZeroWitness :=
  ⟨0, 16660, 15812642284523219010, 8357531243789762760, [1,8,23,2,4,1352,64,256,160,544,32]⟩
theorem profileZeroWitness84_checked : profileZeroWitness84.check = true := by decide
def profileZeroWitness85 : ProfileZeroWitness :=
  ⟨0, 16680, 15807857282932467795, 8357531242779182220, [1,13,23,2,4,1919,98,256,160,544,32]⟩
theorem profileZeroWitness85_checked : profileZeroWitness85.check = true := by decide
def profileZeroWitness86 : ProfileZeroWitness :=
  ⟨0, 16705, 15804197906374846560, 8357531241517008960, [1,8,18,2,4,1088,64,256,130,544,32]⟩
theorem profileZeroWitness86_checked : profileZeroWitness86.check = true := by decide
def profileZeroWitness87 : ProfileZeroWitness :=
  ⟨0, 16770, 15799412904784095345, 8357531240506428420, [1,13,18,2,4,1650,98,256,130,544,32]⟩
theorem profileZeroWitness87_checked : profileZeroWitness87.check = true := by decide
def profileZeroWitness88 : ProfileZeroWitness :=
  ⟨0, 16920, 15740304319298332995, 8357465015394043020, [1,11,23,2,4,1535,64,292,160,544,32]⟩
theorem profileZeroWitness88_checked : profileZeroWitness88.check = true := by decide
def profileZeroWitness89 : ProfileZeroWitness :=
  ⟨0, 16932, 15736082259071199570, 8357465014383493320, [1,14,23,2,4,1992,98,292,160,544,32]⟩
theorem profileZeroWitness89_checked : profileZeroWitness89.check = true := by decide
def profileZeroWitness90 : ProfileZeroWitness :=
  ⟨0, 16962, 15731859941149960545, 8357465013121289220, [1,11,18,2,4,1232,64,292,130,544,32]⟩
theorem profileZeroWitness90_checked : profileZeroWitness90.check = true := by decide
def profileZeroWitness91 : ProfileZeroWitness :=
  ⟨0, 17025, 15727637880922827120, 8357465012110739520, [1,14,18,2,4,1762,98,292,130,544,32]⟩
theorem profileZeroWitness91_checked : profileZeroWitness91.check = true := by decide
def profileZeroWitness92 : ProfileZeroWitness :=
  ⟨0, 17425, 15669087856000887360, 8357382294323972160, [1,8,20,2,4,1280,64,256,160,516,32]⟩
theorem profileZeroWitness92_checked : profileZeroWitness92.check = true := by decide
def profileZeroWitness93 : ProfileZeroWitness :=
  ⟨0, 17442, 15664302854410136145, 8357382293313391620, [1,13,20,2,4,1808,98,256,160,516,32]⟩
theorem profileZeroWitness93_checked : profileZeroWitness93.check = true := by decide
def profileZeroWitness94 : ProfileZeroWitness :=
  ⟨0, 17476, 15659517560766065250, 8357382292051149000, [1,8,17,2,4,1032,64,256,130,516,32]⟩
theorem profileZeroWitness94_checked : profileZeroWitness94.check = true := by decide
def profileZeroWitness95 : ProfileZeroWitness :=
  ⟨0, 17544, 15654732559175314035, 8357382291040568460, [1,13,17,2,4,1565,98,256,130,516,32]⟩
theorem profileZeroWitness95_checked : profileZeroWitness95.check = true := by decide
def profileZeroWitness96 : ProfileZeroWitness :=
  ⟨0, 18450, 15596749890776001345, 8357316065928252420, [1,11,20,2,4,1460,64,292,160,516,32]⟩
theorem profileZeroWitness96_checked : profileZeroWitness96.check = true := by decide
def profileZeroWitness97 : ProfileZeroWitness :=
  ⟨0, 18465, 15592527830548867920, 8357316064917702720, [1,14,20,2,4,1956,98,292,160,516,32]⟩
theorem profileZeroWitness97_checked : profileZeroWitness97.check = true := by decide
def profileZeroWitness98 : ProfileZeroWitness :=
  ⟨0, 18504, 15587179595541179235, 8357316063655429260, [1,11,17,2,4,1179,64,292,130,516,32]⟩
theorem profileZeroWitness98_checked : profileZeroWitness98.check = true := by decide
def profileZeroWitness99 : ProfileZeroWitness :=
  ⟨0, 18564, 15582957535314045810, 8357316062644879560, [1,14,17,2,4,1678,98,292,130,516,32]⟩
theorem profileZeroWitness99_checked : profileZeroWitness99.check = true := by decide
def profileZeroWitness100 : ProfileZeroWitness :=
  ⟨0, 20485, 13363772620989918240, 8314484848174059840, [1,8,2,20,4,64,1280,256,160,32,516]⟩
theorem profileZeroWitness100_checked : profileZeroWitness100.check = true := by decide
def profileZeroWitness101 : ProfileZeroWitness :=
  ⟨0, 20490, 13358987619399167025, 8314484847163479300, [1,13,2,20,4,98,1808,256,160,32,516]⟩
theorem profileZeroWitness101_checked : profileZeroWitness101.check = true := by decide
def profileZeroWitness102 : ProfileZeroWitness :=
  ⟨0, 20560, 13344632030520274020, 8314484845883480520, [1,8,2,17,4,64,1032,256,130,32,516]⟩
theorem profileZeroWitness102_checked : profileZeroWitness102.check = true := by decide
def profileZeroWitness103 : ProfileZeroWitness :=
  ⟨0, 20640, 13339847028929522805, 8314484844872899980, [1,13,2,17,4,98,1565,256,130,32,516]⟩
theorem profileZeroWitness103_checked : profileZeroWitness103.check = true := by decide
def profileZeroWitness104 : ProfileZeroWitness :=
  ⟨0, 21760, 8444640869154942120, 8314297201625586120, [1,8,2,4,17,64,256,1032,32,130,516]⟩
theorem profileZeroWitness104_checked : profileZeroWitness104.check = true := by decide
def profileZeroWitness105 : ProfileZeroWitness :=
  ⟨1, 21845, 18282773019571545120, 17863200012815733120, [1,16,2,4,8,128,512,1024,32,64,256]⟩
theorem profileZeroWitness105_checked : profileZeroWitness105.check = true := by decide
def profileZeroWitness106 : ProfileZeroWitness :=
  ⟨1, 21930, 18282754254859433265, 17863200010794572040, [1,25,2,4,8,194,772,1024,32,64,256]⟩
theorem profileZeroWitness106_checked : profileZeroWitness106.check = true := by decide
def profileZeroWitness107 : ProfileZeroWitness :=
  ⟨0, 23040, 8439855867564190905, 8314230973229866380, [1,13,2,4,17,98,256,1565,32,130,516]⟩
theorem profileZeroWitness107_checked : profileZeroWitness107.check = true := by decide
def profileZeroWitness108 : ProfileZeroWitness :=
  ⟨1, 23130, 18277987871954134065, 17863067556024293640, [1,21,2,4,8,162,512,1288,32,64,256]⟩
theorem profileZeroWitness108_checked : profileZeroWitness108.check = true := by decide
def profileZeroWitness109 : ProfileZeroWitness :=
  ⟨1, 23205, 18277969253270910240, 17863067554003194240, [1,28,2,4,8,224,772,1288,32,64,256]⟩
theorem profileZeroWitness109_checked : profileZeroWitness109.check = true := by decide
def profileZeroWitness110 : ProfileZeroWitness :=
  ⟨0, 24582, 13291434655765032225, 8297596091869450500, [1,11,2,20,4,64,1460,292,160,32,516]⟩
theorem profileZeroWitness110_checked : profileZeroWitness110.check = true := by decide
def profileZeroWitness111 : ProfileZeroWitness :=
  ⟨0, 24585, 13287212595537898800, 8297596090866734400, [1,14,2,20,4,98,1956,292,160,32,516]⟩
theorem profileZeroWitness111_checked : profileZeroWitness111.check = true := by decide
def profileZeroWitness112 : ProfileZeroWitness :=
  ⟨0, 24672, 13272294065295388005, 8297596089578871180, [1,11,2,17,4,64,1179,292,130,32,516]⟩
theorem profileZeroWitness112_checked : profileZeroWitness112.check = true := by decide
def profileZeroWitness113 : ProfileZeroWitness :=
  ⟨0, 24720, 13268072005068254580, 8297596088576155080, [1,14,2,17,4,98,1678,292,130,32,516]⟩
theorem profileZeroWitness113_checked : profileZeroWitness113.check = true := by decide
def profileZeroWitness114 : ProfileZeroWitness :=
  ⟨0, 26112, 8372302903930056105, 8297408445320976780, [1,11,2,4,17,64,292,1179,32,130,516]⟩
theorem profileZeroWitness114_checked : profileZeroWitness114.check = true := by decide
def profileZeroWitness115 : ProfileZeroWitness :=
  ⟨1, 26214, 18210432846767154465, 17829422500206514440, [1,19,2,4,8,128,548,1096,32,64,256]⟩
theorem profileZeroWitness115_checked : profileZeroWitness115.check = true := by decide
def profileZeroWitness116 : ProfileZeroWitness :=
  ⟨1, 26265, 18210416289668232240, 17829422498201082240, [1,26,2,4,8,194,800,1096,32,64,256]⟩
theorem profileZeroWitness116_checked : profileZeroWitness116.check = true := by decide
def profileZeroWitness117 : ProfileZeroWitness :=
  ⟨0, 26880, 8368080843702922680, 8297342732321332680, [1,14,2,4,17,98,292,1678,32,130,516]⟩
theorem profileZeroWitness117_checked : profileZeroWitness117.check = true := by decide
def profileZeroWitness118 : ProfileZeroWitness :=
  ⟨1, 26985, 18206210657692968240, 17829291074207226240, [1,22,2,4,8,162,548,1344,32,64,256]⟩
theorem profileZeroWitness118_checked : profileZeroWitness118.check = true := by decide
def profileZeroWitness119 : ProfileZeroWitness :=
  ⟨1, 27030, 18206194229443064865, 17829291072201855240, [1,31,2,4,8,224,800,1344,32,64,256]⟩
theorem profileZeroWitness119_checked : profileZeroWitness119.check = true := by decide
def profileZeroWitness120 : ProfileZeroWitness :=
  ⟨0, 33048, 14659456901152176195, 4017187102911431820, [1,15,23,2,4,2015,96,288,160,544,32]⟩
theorem profileZeroWitness120_checked : profileZeroWitness120.check = true := by decide
def profileZeroWitness121 : ProfileZeroWitness :=
  ⟨0, 33060, 14655234840925042770, 4017187101900882120, [1,10,23,2,4,1512,66,288,160,544,32]⟩
theorem profileZeroWitness121_checked : profileZeroWitness121.check = true := by decide
def profileZeroWitness122 : ProfileZeroWitness :=
  ⟨0, 33090, 14651012523003803745, 4017187100638678020, [1,15,18,2,4,1776,96,288,130,544,32]⟩
theorem profileZeroWitness122_checked : profileZeroWitness122.check = true := by decide
def profileZeroWitness123 : ProfileZeroWitness :=
  ⟨0, 33153, 14646790462776670320, 4017187099628128320, [1,10,18,2,4,1218,66,288,130,544,32]⟩
theorem profileZeroWitness123_checked : profileZeroWitness123.check = true := by decide
def profileZeroWitness124 : ProfileZeroWitness :=
  ⟨0, 33300, 14587681877290907970, 4017120874515742920, [1,12,23,2,4,1896,96,260,160,544,32]⟩
theorem profileZeroWitness124_checked : profileZeroWitness124.check = true := by decide
def profileZeroWitness125 : ProfileZeroWitness :=
  ⟨0, 33320, 14582896875700156755, 4017120873505162380, [1,9,23,2,4,1375,66,260,160,544,32]⟩
theorem profileZeroWitness125_checked : profileZeroWitness125.check = true := by decide
def profileZeroWitness126 : ProfileZeroWitness :=
  ⟨0, 33345, 14579237499142535520, 4017120872242989120, [1,12,18,2,4,1632,96,260,130,544,32]⟩
theorem profileZeroWitness126_checked : profileZeroWitness126.check = true := by decide
def profileZeroWitness127 : ProfileZeroWitness :=
  ⟨0, 33410, 14574452497551784305, 4017120871232408580, [1,9,18,2,4,1106,66,260,130,544,32]⟩
theorem profileZeroWitness127_checked : profileZeroWitness127.check = true := by decide
def profileZeroWitness128 : ProfileZeroWitness :=
  ⟨0, 33810, 14515902472629844545, 4017038153445641220, [1,15,20,2,4,1968,96,288,160,516,32]⟩
theorem profileZeroWitness128_checked : profileZeroWitness128.check = true := by decide
def profileZeroWitness129 : ProfileZeroWitness :=
  ⟨0, 33825, 14511680412402711120, 4017038152435091520, [1,10,20,2,4,1440,66,288,160,516,32]⟩
theorem profileZeroWitness129_checked : profileZeroWitness129.check = true := by decide
def profileZeroWitness130 : ProfileZeroWitness :=
  ⟨0, 33864, 14506332177395022435, 4017038151172818060, [1,15,17,2,4,1695,96,288,130,516,32]⟩
theorem profileZeroWitness130_checked : profileZeroWitness130.check = true := by decide
def profileZeroWitness131 : ProfileZeroWitness :=
  ⟨0, 33924, 14502110117167889010, 4017038150162268360, [1,10,17,2,4,1162,66,288,130,516,32]⟩
theorem profileZeroWitness131_checked : profileZeroWitness131.check = true := by decide
def profileZeroWitness132 : ProfileZeroWitness :=
  ⟨0, 34833, 14444127448768576320, 4016971925049952320, [1,12,20,2,4,1796,96,260,160,516,32]⟩
theorem profileZeroWitness132_checked : profileZeroWitness132.check = true := by decide
def profileZeroWitness133 : ProfileZeroWitness :=
  ⟨0, 34850, 14439342447177825105, 4016971924039371780, [1,9,20,2,4,1300,66,260,160,516,32]⟩
theorem profileZeroWitness133_checked : profileZeroWitness133.check = true := by decide
def profileZeroWitness134 : ProfileZeroWitness :=
  ⟨0, 34884, 14434557153533754210, 4016971922777129160, [1,12,17,2,4,1548,96,260,130,516,32]⟩
theorem profileZeroWitness134_checked : profileZeroWitness134.check = true := by decide
def profileZeroWitness135 : ProfileZeroWitness :=
  ⟨0, 34952, 14429772151943002995, 4016971921766548620, [1,9,17,2,4,1049,66,260,130,516,32]⟩
theorem profileZeroWitness135_checked : profileZeroWitness135.check = true := by decide
def profileZeroWitness136 : ProfileZeroWitness :=
  ⟨0, 36870, 12210587237618875425, 3990963235196784900, [1,15,2,20,4,96,1968,288,160,32,516]⟩
theorem profileZeroWitness136_checked : profileZeroWitness136.check = true := by decide
def profileZeroWitness137 : ProfileZeroWitness :=
  ⟨0, 36873, 12206365177391742000, 3990963234194068800, [1,10,2,20,4,66,1440,288,160,32,516]⟩
theorem profileZeroWitness137_checked : profileZeroWitness137.check = true := by decide
def profileZeroWitness138 : ProfileZeroWitness :=
  ⟨0, 36960, 12191446647149231205, 3990963232906205580, [1,15,2,17,4,96,1695,288,130,32,516]⟩
theorem profileZeroWitness138_checked : profileZeroWitness138.check = true := by decide
def profileZeroWitness139 : ProfileZeroWitness :=
  ⟨0, 37008, 12187224586922097780, 3990963231903489480, [1,10,2,17,4,66,1162,288,130,32,516]⟩
theorem profileZeroWitness139_checked : profileZeroWitness139.check = true := by decide
def profileZeroWitness140 : ProfileZeroWitness :=
  ⟨0, 38400, 7291455485783899305, 3990841300645239180, [1,15,2,4,17,96,288,1695,32,130,516]⟩
theorem profileZeroWitness140_checked : profileZeroWitness140.check = true := by decide
def profileZeroWitness141 : ProfileZeroWitness :=
  ⟨1, 38505, 17129568872025384240, 9216288212860410240, [1,30,2,4,8,226,804,1352,32,64,256]⟩
theorem profileZeroWitness141_checked : profileZeroWitness141.check = true := by decide
def profileZeroWitness142 : ProfileZeroWitness :=
  ⟨1, 38550, 17129552443775480865, 9216288210855039240, [1,23,2,4,8,160,544,1352,32,64,256]⟩
theorem profileZeroWitness142_checked : profileZeroWitness142.check = true := by decide
def profileZeroWitness143 : ProfileZeroWitness :=
  ⟨0, 39168, 7287233425556765880, 3990775587645595080, [1,10,2,4,17,66,288,1162,32,130,516]⟩
theorem profileZeroWitness143_checked : profileZeroWitness143.check = true := by decide
def profileZeroWitness144 : ProfileZeroWitness :=
  ⟨1, 39270, 17125346811800216865, 9216156786861183240, [1,27,2,4,8,192,804,1088,32,64,256]⟩
theorem profileZeroWitness144_checked : profileZeroWitness144.check = true := by decide
def profileZeroWitness145 : ProfileZeroWitness :=
  ⟨1, 39321, 17125330254701294640, 9216156784855751040, [1,18,2,4,8,130,544,1088,32,64,256]⟩
theorem profileZeroWitness145_checked : profileZeroWitness145.check = true := by decide
def profileZeroWitness146 : ProfileZeroWitness :=
  ⟨0, 40965, 12138812213757607200, 3974074478900040000, [1,12,2,20,4,96,1796,260,160,32,516]⟩
theorem profileZeroWitness146_checked : profileZeroWitness146.check = true := by decide
def profileZeroWitness147 : ProfileZeroWitness :=
  ⟨0, 40970, 12134027212166855985, 3974074477889459460, [1,9,2,20,4,66,1300,260,160,32,516]⟩
theorem profileZeroWitness147_checked : profileZeroWitness147.check = true := by decide
def profileZeroWitness148 : ProfileZeroWitness :=
  ⟨0, 41040, 12119671623287962980, 3974074476609460680, [1,12,2,17,4,96,1548,260,130,32,516]⟩
theorem profileZeroWitness148_checked : profileZeroWitness148.check = true := by decide
def profileZeroWitness149 : ProfileZeroWitness :=
  ⟨0, 41120, 12114886621697211765, 3974074475598880140, [1,9,2,17,4,66,1049,260,130,32,516]⟩
theorem profileZeroWitness149_checked : profileZeroWitness149.check = true := by decide
def profileZeroWitness150 : ProfileZeroWitness :=
  ⟨0, 42240, 7219680461922631080, 3973953059736705480, [1,12,2,4,17,96,260,1548,32,130,516]⟩
theorem profileZeroWitness150_checked : profileZeroWitness150.check = true := by decide
def profileZeroWitness151 : ProfileZeroWitness :=
  ⟨1, 42330, 17057793848197538865, 9182511731059071240, [1,29,2,4,8,226,768,1280,32,64,256]⟩
theorem profileZeroWitness151_checked : profileZeroWitness151.check = true := by decide
def profileZeroWitness152 : ProfileZeroWitness :=
  ⟨1, 42405, 17057775229514315040, 9182511729037971840, [1,20,2,4,8,160,516,1280,32,64,256]⟩
theorem profileZeroWitness152_checked : profileZeroWitness152.check = true := by decide
def profileZeroWitness153 : ProfileZeroWitness :=
  ⟨0, 43520, 7214895460331879865, 3973886831340985740, [1,9,2,4,17,66,260,1049,32,130,516]⟩
theorem profileZeroWitness153_checked : profileZeroWitness153.check = true := by decide
def profileZeroWitness154 : ProfileZeroWitness :=
  ⟨1, 43605, 17053008846609015840, 9182379274267693440, [1,24,2,4,8,192,768,1032,32,64,256]⟩
theorem profileZeroWitness154_checked : profileZeroWitness154.check = true := by decide
def profileZeroWitness155 : ProfileZeroWitness :=
  ⟨1, 43690, 17052990081896903985, 9182379272246532360, [1,17,2,4,8,130,516,1032,32,64,256]⟩
theorem profileZeroWitness155_checked : profileZeroWitness155.check = true := by decide
def profileZeroWitness156 : ProfileZeroWitness :=
  ⟨0, 49155, 10986173214651252240, 3636299352933684240, [1,2,12,20,4,96,160,32,1796,260,516]⟩
theorem profileZeroWitness156_checked : profileZeroWitness156.check = true := by decide
def profileZeroWitness157 : ProfileZeroWitness :=
  ⟨0, 49164, 10976603211469749810, 3636299351797272660, [1,2,9,20,4,66,160,32,1300,260,516]⟩
theorem profileZeroWitness157_checked : profileZeroWitness157.check = true := by decide
def profileZeroWitness158 : ProfileZeroWitness :=
  ⟨0, 49200, 10967032624181608020, 3636299350643104920, [1,2,12,17,4,96,130,32,1548,260,516]⟩
theorem profileZeroWitness158_checked : profileZeroWitness158.check = true := by decide
def profileZeroWitness159 : ProfileZeroWitness :=
  ⟨0, 49344, 10957462621000105590, 3636299349506693340, [1,2,9,17,4,66,130,32,1049,260,516]⟩
theorem profileZeroWitness159_checked : profileZeroWitness159.check = true := by decide
def profileZeroWitness160 : ProfileZeroWitness :=
  ⟨0, 49920, 6067041462816276120, 3636186179981728920, [1,2,12,4,17,96,32,130,260,1548,516]⟩
theorem profileZeroWitness160_checked : profileZeroWitness160.check = true := by decide
def profileZeroWitness161 : ProfileZeroWitness :=
  ⟨1, 49980, 15905154849375278130, 8574532997019209880, [1,2,29,4,8,226,32,64,768,1280,256]⟩
theorem profileZeroWitness161_checked : profileZeroWitness161.check = true := by decide
def profileZeroWitness162 : ProfileZeroWitness :=
  ⟨1, 50115, 15905117612008830480, 8574532994746456080, [1,2,20,4,8,160,32,64,516,1280,256]⟩
theorem profileZeroWitness162_checked : profileZeroWitness162.check = true := by decide
def profileZeroWitness163 : ProfileZeroWitness :=
  ⟨0, 52224, 6057471459634773690, 3636111705248798940, [1,2,9,4,17,66,32,130,260,1049,516]⟩
theorem profileZeroWitness163_checked : profileZeroWitness163.check = true := by decide
def profileZeroWitness164 : ProfileZeroWitness :=
  ⟨1, 52275, 15895584846198232080, 8574384047553419280, [1,2,24,4,8,192,32,64,768,1032,256]⟩
theorem profileZeroWitness164_checked : profileZeroWitness164.check = true := by decide
def profileZeroWitness165 : ProfileZeroWitness :=
  ⟨1, 52428, 15895547316774008370, 8574384045280596120, [1,2,17,4,8,130,32,64,516,1032,256]⟩
theorem profileZeroWitness165_checked : profileZeroWitness165.check = true := by decide
def profileZeroWitness166 : ProfileZeroWitness :=
  ⟨0, 61440, 3607513407803686620, 3607513407803686620, [1,2,4,9,17,32,66,130,260,516,1049]⟩
theorem profileZeroWitness166_checked : profileZeroWitness166.check = true := by decide
def profileZeroWitness167 : ProfileZeroWitness :=
  ⟨1, 61455, 13445626794939789840, 8526720218495988240, [1,2,4,24,8,32,192,64,768,256,1032]⟩
theorem profileZeroWitness167_checked : profileZeroWitness167.check = true := by decide
def profileZeroWitness168 : ProfileZeroWitness :=
  ⟨1, 61680, 13445551736091342420, 8526720216205408920, [1,2,4,17,8,32,130,64,516,256,1032]⟩
theorem profileZeroWitness168_checked : profileZeroWitness168.check = true := by decide
def profileZeroWitness169 : ProfileZeroWitness :=
  ⟨1, 65280, 8526495043095935640, 8526495043095935640, [1,2,4,8,17,32,64,130,256,516,1032]⟩
theorem profileZeroWitness169_checked : profileZeroWitness169.check = true := by decide
def profileZeroWitness170 : ProfileZeroWitness :=
  ⟨2, 65535, 18364758544493064720, 18364758544493064720, [1,2,4,8,16,32,64,128,256,512,1024]⟩
theorem profileZeroWitness170_checked : profileZeroWitness170.check = true := by decide
def profileZeroBank : List ProfileZeroWitness := [
  profileZeroWitness0,
  profileZeroWitness1,
  profileZeroWitness2,
  profileZeroWitness3,
  profileZeroWitness4,
  profileZeroWitness5,
  profileZeroWitness6,
  profileZeroWitness7,
  profileZeroWitness8,
  profileZeroWitness9,
  profileZeroWitness10,
  profileZeroWitness11,
  profileZeroWitness12,
  profileZeroWitness13,
  profileZeroWitness14,
  profileZeroWitness15,
  profileZeroWitness16,
  profileZeroWitness17,
  profileZeroWitness18,
  profileZeroWitness19,
  profileZeroWitness20,
  profileZeroWitness21,
  profileZeroWitness22,
  profileZeroWitness23,
  profileZeroWitness24,
  profileZeroWitness25,
  profileZeroWitness26,
  profileZeroWitness27,
  profileZeroWitness28,
  profileZeroWitness29,
  profileZeroWitness30,
  profileZeroWitness31,
  profileZeroWitness32,
  profileZeroWitness33,
  profileZeroWitness34,
  profileZeroWitness35,
  profileZeroWitness36,
  profileZeroWitness37,
  profileZeroWitness38,
  profileZeroWitness39,
  profileZeroWitness40,
  profileZeroWitness41,
  profileZeroWitness42,
  profileZeroWitness43,
  profileZeroWitness44,
  profileZeroWitness45,
  profileZeroWitness46,
  profileZeroWitness47,
  profileZeroWitness48,
  profileZeroWitness49,
  profileZeroWitness50,
  profileZeroWitness51,
  profileZeroWitness52,
  profileZeroWitness53,
  profileZeroWitness54,
  profileZeroWitness55,
  profileZeroWitness56,
  profileZeroWitness57,
  profileZeroWitness58,
  profileZeroWitness59,
  profileZeroWitness60,
  profileZeroWitness61,
  profileZeroWitness62,
  profileZeroWitness63,
  profileZeroWitness64,
  profileZeroWitness65,
  profileZeroWitness66,
  profileZeroWitness67,
  profileZeroWitness68,
  profileZeroWitness69,
  profileZeroWitness70,
  profileZeroWitness71,
  profileZeroWitness72,
  profileZeroWitness73,
  profileZeroWitness74,
  profileZeroWitness75,
  profileZeroWitness76,
  profileZeroWitness77,
  profileZeroWitness78,
  profileZeroWitness79,
  profileZeroWitness80,
  profileZeroWitness81,
  profileZeroWitness82,
  profileZeroWitness83,
  profileZeroWitness84,
  profileZeroWitness85,
  profileZeroWitness86,
  profileZeroWitness87,
  profileZeroWitness88,
  profileZeroWitness89,
  profileZeroWitness90,
  profileZeroWitness91,
  profileZeroWitness92,
  profileZeroWitness93,
  profileZeroWitness94,
  profileZeroWitness95,
  profileZeroWitness96,
  profileZeroWitness97,
  profileZeroWitness98,
  profileZeroWitness99,
  profileZeroWitness100,
  profileZeroWitness101,
  profileZeroWitness102,
  profileZeroWitness103,
  profileZeroWitness104,
  profileZeroWitness105,
  profileZeroWitness106,
  profileZeroWitness107,
  profileZeroWitness108,
  profileZeroWitness109,
  profileZeroWitness110,
  profileZeroWitness111,
  profileZeroWitness112,
  profileZeroWitness113,
  profileZeroWitness114,
  profileZeroWitness115,
  profileZeroWitness116,
  profileZeroWitness117,
  profileZeroWitness118,
  profileZeroWitness119,
  profileZeroWitness120,
  profileZeroWitness121,
  profileZeroWitness122,
  profileZeroWitness123,
  profileZeroWitness124,
  profileZeroWitness125,
  profileZeroWitness126,
  profileZeroWitness127,
  profileZeroWitness128,
  profileZeroWitness129,
  profileZeroWitness130,
  profileZeroWitness131,
  profileZeroWitness132,
  profileZeroWitness133,
  profileZeroWitness134,
  profileZeroWitness135,
  profileZeroWitness136,
  profileZeroWitness137,
  profileZeroWitness138,
  profileZeroWitness139,
  profileZeroWitness140,
  profileZeroWitness141,
  profileZeroWitness142,
  profileZeroWitness143,
  profileZeroWitness144,
  profileZeroWitness145,
  profileZeroWitness146,
  profileZeroWitness147,
  profileZeroWitness148,
  profileZeroWitness149,
  profileZeroWitness150,
  profileZeroWitness151,
  profileZeroWitness152,
  profileZeroWitness153,
  profileZeroWitness154,
  profileZeroWitness155,
  profileZeroWitness156,
  profileZeroWitness157,
  profileZeroWitness158,
  profileZeroWitness159,
  profileZeroWitness160,
  profileZeroWitness161,
  profileZeroWitness162,
  profileZeroWitness163,
  profileZeroWitness164,
  profileZeroWitness165,
  profileZeroWitness166,
  profileZeroWitness167,
  profileZeroWitness168,
  profileZeroWitness169,
  profileZeroWitness170
]
theorem profileZeroBank_checked : ∀ w ∈ profileZeroBank, w.check = true := by
  simp only [profileZeroBank, List.mem_cons, List.not_mem_nil, or_false]
  intro w h
  rcases h with h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h
  · subst w; exact profileZeroWitness0_checked
  · subst w; exact profileZeroWitness1_checked
  · subst w; exact profileZeroWitness2_checked
  · subst w; exact profileZeroWitness3_checked
  · subst w; exact profileZeroWitness4_checked
  · subst w; exact profileZeroWitness5_checked
  · subst w; exact profileZeroWitness6_checked
  · subst w; exact profileZeroWitness7_checked
  · subst w; exact profileZeroWitness8_checked
  · subst w; exact profileZeroWitness9_checked
  · subst w; exact profileZeroWitness10_checked
  · subst w; exact profileZeroWitness11_checked
  · subst w; exact profileZeroWitness12_checked
  · subst w; exact profileZeroWitness13_checked
  · subst w; exact profileZeroWitness14_checked
  · subst w; exact profileZeroWitness15_checked
  · subst w; exact profileZeroWitness16_checked
  · subst w; exact profileZeroWitness17_checked
  · subst w; exact profileZeroWitness18_checked
  · subst w; exact profileZeroWitness19_checked
  · subst w; exact profileZeroWitness20_checked
  · subst w; exact profileZeroWitness21_checked
  · subst w; exact profileZeroWitness22_checked
  · subst w; exact profileZeroWitness23_checked
  · subst w; exact profileZeroWitness24_checked
  · subst w; exact profileZeroWitness25_checked
  · subst w; exact profileZeroWitness26_checked
  · subst w; exact profileZeroWitness27_checked
  · subst w; exact profileZeroWitness28_checked
  · subst w; exact profileZeroWitness29_checked
  · subst w; exact profileZeroWitness30_checked
  · subst w; exact profileZeroWitness31_checked
  · subst w; exact profileZeroWitness32_checked
  · subst w; exact profileZeroWitness33_checked
  · subst w; exact profileZeroWitness34_checked
  · subst w; exact profileZeroWitness35_checked
  · subst w; exact profileZeroWitness36_checked
  · subst w; exact profileZeroWitness37_checked
  · subst w; exact profileZeroWitness38_checked
  · subst w; exact profileZeroWitness39_checked
  · subst w; exact profileZeroWitness40_checked
  · subst w; exact profileZeroWitness41_checked
  · subst w; exact profileZeroWitness42_checked
  · subst w; exact profileZeroWitness43_checked
  · subst w; exact profileZeroWitness44_checked
  · subst w; exact profileZeroWitness45_checked
  · subst w; exact profileZeroWitness46_checked
  · subst w; exact profileZeroWitness47_checked
  · subst w; exact profileZeroWitness48_checked
  · subst w; exact profileZeroWitness49_checked
  · subst w; exact profileZeroWitness50_checked
  · subst w; exact profileZeroWitness51_checked
  · subst w; exact profileZeroWitness52_checked
  · subst w; exact profileZeroWitness53_checked
  · subst w; exact profileZeroWitness54_checked
  · subst w; exact profileZeroWitness55_checked
  · subst w; exact profileZeroWitness56_checked
  · subst w; exact profileZeroWitness57_checked
  · subst w; exact profileZeroWitness58_checked
  · subst w; exact profileZeroWitness59_checked
  · subst w; exact profileZeroWitness60_checked
  · subst w; exact profileZeroWitness61_checked
  · subst w; exact profileZeroWitness62_checked
  · subst w; exact profileZeroWitness63_checked
  · subst w; exact profileZeroWitness64_checked
  · subst w; exact profileZeroWitness65_checked
  · subst w; exact profileZeroWitness66_checked
  · subst w; exact profileZeroWitness67_checked
  · subst w; exact profileZeroWitness68_checked
  · subst w; exact profileZeroWitness69_checked
  · subst w; exact profileZeroWitness70_checked
  · subst w; exact profileZeroWitness71_checked
  · subst w; exact profileZeroWitness72_checked
  · subst w; exact profileZeroWitness73_checked
  · subst w; exact profileZeroWitness74_checked
  · subst w; exact profileZeroWitness75_checked
  · subst w; exact profileZeroWitness76_checked
  · subst w; exact profileZeroWitness77_checked
  · subst w; exact profileZeroWitness78_checked
  · subst w; exact profileZeroWitness79_checked
  · subst w; exact profileZeroWitness80_checked
  · subst w; exact profileZeroWitness81_checked
  · subst w; exact profileZeroWitness82_checked
  · subst w; exact profileZeroWitness83_checked
  · subst w; exact profileZeroWitness84_checked
  · subst w; exact profileZeroWitness85_checked
  · subst w; exact profileZeroWitness86_checked
  · subst w; exact profileZeroWitness87_checked
  · subst w; exact profileZeroWitness88_checked
  · subst w; exact profileZeroWitness89_checked
  · subst w; exact profileZeroWitness90_checked
  · subst w; exact profileZeroWitness91_checked
  · subst w; exact profileZeroWitness92_checked
  · subst w; exact profileZeroWitness93_checked
  · subst w; exact profileZeroWitness94_checked
  · subst w; exact profileZeroWitness95_checked
  · subst w; exact profileZeroWitness96_checked
  · subst w; exact profileZeroWitness97_checked
  · subst w; exact profileZeroWitness98_checked
  · subst w; exact profileZeroWitness99_checked
  · subst w; exact profileZeroWitness100_checked
  · subst w; exact profileZeroWitness101_checked
  · subst w; exact profileZeroWitness102_checked
  · subst w; exact profileZeroWitness103_checked
  · subst w; exact profileZeroWitness104_checked
  · subst w; exact profileZeroWitness105_checked
  · subst w; exact profileZeroWitness106_checked
  · subst w; exact profileZeroWitness107_checked
  · subst w; exact profileZeroWitness108_checked
  · subst w; exact profileZeroWitness109_checked
  · subst w; exact profileZeroWitness110_checked
  · subst w; exact profileZeroWitness111_checked
  · subst w; exact profileZeroWitness112_checked
  · subst w; exact profileZeroWitness113_checked
  · subst w; exact profileZeroWitness114_checked
  · subst w; exact profileZeroWitness115_checked
  · subst w; exact profileZeroWitness116_checked
  · subst w; exact profileZeroWitness117_checked
  · subst w; exact profileZeroWitness118_checked
  · subst w; exact profileZeroWitness119_checked
  · subst w; exact profileZeroWitness120_checked
  · subst w; exact profileZeroWitness121_checked
  · subst w; exact profileZeroWitness122_checked
  · subst w; exact profileZeroWitness123_checked
  · subst w; exact profileZeroWitness124_checked
  · subst w; exact profileZeroWitness125_checked
  · subst w; exact profileZeroWitness126_checked
  · subst w; exact profileZeroWitness127_checked
  · subst w; exact profileZeroWitness128_checked
  · subst w; exact profileZeroWitness129_checked
  · subst w; exact profileZeroWitness130_checked
  · subst w; exact profileZeroWitness131_checked
  · subst w; exact profileZeroWitness132_checked
  · subst w; exact profileZeroWitness133_checked
  · subst w; exact profileZeroWitness134_checked
  · subst w; exact profileZeroWitness135_checked
  · subst w; exact profileZeroWitness136_checked
  · subst w; exact profileZeroWitness137_checked
  · subst w; exact profileZeroWitness138_checked
  · subst w; exact profileZeroWitness139_checked
  · subst w; exact profileZeroWitness140_checked
  · subst w; exact profileZeroWitness141_checked
  · subst w; exact profileZeroWitness142_checked
  · subst w; exact profileZeroWitness143_checked
  · subst w; exact profileZeroWitness144_checked
  · subst w; exact profileZeroWitness145_checked
  · subst w; exact profileZeroWitness146_checked
  · subst w; exact profileZeroWitness147_checked
  · subst w; exact profileZeroWitness148_checked
  · subst w; exact profileZeroWitness149_checked
  · subst w; exact profileZeroWitness150_checked
  · subst w; exact profileZeroWitness151_checked
  · subst w; exact profileZeroWitness152_checked
  · subst w; exact profileZeroWitness153_checked
  · subst w; exact profileZeroWitness154_checked
  · subst w; exact profileZeroWitness155_checked
  · subst w; exact profileZeroWitness156_checked
  · subst w; exact profileZeroWitness157_checked
  · subst w; exact profileZeroWitness158_checked
  · subst w; exact profileZeroWitness159_checked
  · subst w; exact profileZeroWitness160_checked
  · subst w; exact profileZeroWitness161_checked
  · subst w; exact profileZeroWitness162_checked
  · subst w; exact profileZeroWitness163_checked
  · subst w; exact profileZeroWitness164_checked
  · subst w; exact profileZeroWitness165_checked
  · subst w; exact profileZeroWitness166_checked
  · subst w; exact profileZeroWitness167_checked
  · subst w; exact profileZeroWitness168_checked
  · subst w; exact profileZeroWitness169_checked
  · subst w; exact profileZeroWitness170_checked
end Spin.UniformMixer
