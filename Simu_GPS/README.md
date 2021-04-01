# attenuation_in_block

#运行说明

#在本目录外建一个空目录，任意取名（例如BG_simu-build)
#进入BG_simu-build
#执行cmake ../BG_simu
#执行make
#得到attenuationInBlock文件


#模拟gamma本底
#执行./attenuationInBlock "gamma"
#打开图形化界面（visualization）同时在目录内得到result.root文件
#图形化界面输入：/vis/open OGL 可以看到探测器几何
#输入run/BeamOn 1 可以显示一个事件

#gamma谱实验数据：BG_simu/include/spectra.h
#入射粒子抽样（BG_simu/src/primaryGeneratorAction.cc）： 随机产生半球状分布，能量分布按照spectra.h中的实验数据抽样
#result.root文件中originals是粒子最终通过屏蔽体后在探测器中的能量分布

#需要修改的文件
#入射粒子：BG_simu/src/primaryGeneratorAction.cc
#探测器屏蔽体几何：BG_simu/src/detectorConstruction.cc

#获得能量沉积：src/steppingAction.cc：step->GetTotalEnergyDeposit();
#能量沉积填入originals: src/eventAction.cc: FillOriginalHistogram(fEdep);

