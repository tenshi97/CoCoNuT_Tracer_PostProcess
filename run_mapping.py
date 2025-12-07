import h5reader_3d
import map3d

y = h5reader_3d.hydrof(model='s3.5_envel_new_cut3', yinyang=True)
map3d.map3d(y, y.xzl(), y.xzr(), y.yzl(), y.yzr(), y.zzl(), y.zzr(), -1.2e13, 1.2e13, -1.2e13, 1.2e13, 50, 50, 50, 1e-3)
