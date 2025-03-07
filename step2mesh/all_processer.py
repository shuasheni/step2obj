from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.Bnd import Bnd_Box

from step2mesh.stepfile_processer import STEPFileProcessor


class AllShapeProcessor(STEPFileProcessor):
    """
    这个类的pocess方法会输出全部形状的mesh
    在is_remove中保留全部形状
    """

    def __init__(self, step_file, lin_deflection=0.2, ang_deflection=0.36):
        super().__init__(step_file, lin_deflection=lin_deflection, ang_deflection=ang_deflection)

    def is_remove(self, shape, vertices):
        # 具体条件：保留全部形状
        return False
