from collections import defaultdict

from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.Bnd import Bnd_Box
import numpy as np
from scipy.spatial import KDTree, QhullError
from step2mesh.stepfile_processer import STEPFileProcessor
from scipy.spatial import ConvexHull

def point_to_hull_surface_distance(point, hull):
    """
    计算点到三维凸包表面的最小距离（近似值）

    参数：
        point : numpy数组, 形状(3,) 三维点坐标
        hull : scipy.spatial.ConvexHull对象

    返回：
        float: 点到凸包表面的最小距离（正数）
    """
    min_distance = np.inf
    for eq in hull.equations:
        # 提取平面方程参数 Ax + By + Cz + D = 0
        A, B, C, D = eq

        # 计算有符号距离（无需绝对值）
        signed_distance = A * point[0] + B * point[1] + C * point[2] + D

        # 计算法向量模长
        norm = np.sqrt(A ** 2 + B ** 2 + C ** 2)

        # 计算绝对距离
        distance = abs(signed_distance) / norm

        # 更新最小距离
        if distance < min_distance:
            min_distance = distance

    return min_distance

class BounderBoxProcessor(STEPFileProcessor):

    def __init__(self, step_file,dist_tor = 5, lin_deflection=0.2, ang_deflection=0.36, show_detail=False, scale=1.0, p=0):
        super().__init__(step_file, lin_deflection=lin_deflection, ang_deflection=ang_deflection, show_detail=show_detail)
        self.scale = scale
        self.p = p
        self.cut_box = [9999.0, 9999.0, 9999.0, -9999.0, -9999.0, -9999.0]
        self.bounder_box = [9999.0, 9999.0, 9999.0, -9999.0, -9999.0, -9999.0]
        self.hull_exist = defaultdict(bool)
        self.convex_hulls = []
        self.dist_tor = dist_tor


    def is_point_in_3d_convex_hull(self, point, hull, tolerance=1e-10):
        """
        判断点是否在三维凸包内
        :param point: 待测点，形状为 (3,) 的 numpy 数组
        :param hull: scipy.spatial.ConvexHull 对象
        :param tolerance: 浮点误差容忍度
        :return: True（在内部）或 False
        """

        for eq in hull.equations:
            # 计算点代入平面方程的值：Ax + By + Cz + D
            val = np.dot(eq[:-1], point) + eq[-1]
            if val > -tolerance:  # 点在平面外侧
                return False
        return True

    def compute_convex_hull_all(self):
        convex_hulls = []
        for (c_name, transformed_shape, vertices, face_meshes, result) in self.final_shapes:
            try:
                # 计算三维凸包
                np_vertices = np.array(vertices)
                hull = ConvexHull(np_vertices)
                convex_hulls.append((c_name, hull))
                if self.show_detail:
                    print(f"形状 '{c_name}' 生成凸包")
            except QhullError as e:
                print(f"警告：形状 '{c_name}' 无法生成凸包（原因：{str(e)}）")

        # 按凸包体积从大到小排序
        sorted_convex_hulls = sorted(
            convex_hulls,
            key=lambda item: item[1].volume,
            reverse=True
        )

        self.filter_non_contained_convex_hulls(sorted_convex_hulls)

    def filter_non_contained_convex_hulls(self, sorted_convex_hulls):
        """
        过滤被其他凸包完全包含的凸包
        :param sorted_convex_hulls: 已按体积降序排序的凸包列表
        :return: 过滤后的凸包列表 accepted_convex_hulls
        """

        for shape_name, hull in sorted_convex_hulls:
            # 获取当前凸包的所有顶点
            is_covered = False
            for accepted_hull_name, accepted_hull in self.convex_hulls:
                point_not_covered = False
                for vertex in hull.points[hull.vertices]:
                    if not self.is_point_in_3d_convex_hull(vertex, accepted_hull):
                        point_not_covered = True
                        break
                if not point_not_covered:
                    is_covered = True
                    break

            if self.show_detail:
                print(f"形状 '{shape_name}' 凸包保留? {not is_covered}")

            if not is_covered:
                self.convex_hulls.append((shape_name, hull))
                self.hull_exist[shape_name] = True
            else:
                self.hull_exist[shape_name] = False

    def get_bounding_box(self, shape):
        bounding_box = Bnd_Box()
        brepbndlib.Add(shape, bounding_box)
        return bounding_box

    def check_shape(self, shape, vertices):
        bounding_box = self.get_bounding_box(shape)
        min_point = bounding_box.CornerMin()
        max_point = bounding_box.CornerMax()
        box_num = [min_point.X(), min_point.Y(), min_point.Z(),
                   max_point.X(), max_point.Y(), max_point.Z()]
        if self.show_detail:
            print(f"    Min point: ({min_point.X()}, {min_point.Y()}, {min_point.Z()})")
            print(f"    Max point: ({max_point.X()}, {max_point.Y()}, {max_point.Z()})")
        if box_num[0] < self.bounder_box[0]:
            self.bounder_box[0] = box_num[0]
        if box_num[1] < self.bounder_box[1]:
            self.bounder_box[1] = box_num[1]
        if box_num[2] < self.bounder_box[2]:
            self.bounder_box[2] = box_num[2]
        if box_num[3] > self.bounder_box[3]:
            self.bounder_box[3] = box_num[3]
        if box_num[4] > self.bounder_box[4]:
            self.bounder_box[4] = box_num[4]
        if box_num[5] > self.bounder_box[5]:
            self.bounder_box[5] = box_num[5]
        return box_num

    def generate_face_samples(self, face):
        xl = self.bounder_box[3] - self.bounder_box[0]
        yl = self.bounder_box[4] - self.bounder_box[1]
        zl = self.bounder_box[5] - self.bounder_box[2]
        xs = xl * (1 - self.scale) * 0.5
        ys = yl * (1 - self.scale) * 0.5
        zs = zl * (1 - self.scale) * 0.5

        if face == 'front':
            x = np.linspace(self.bounder_box[0] + xs, self.bounder_box[3] - xs, 10)
            y = np.linspace(self.bounder_box[1], self.bounder_box[4], 10)
            xx, yy = np.meshgrid(x, y)
            zz = np.full_like(xx, self.bounder_box[5])
        elif face == 'back':
            x = np.linspace(self.bounder_box[0] + xs, self.bounder_box[3] - xs, 10)
            y = np.linspace(self.bounder_box[1] + ys, self.bounder_box[4] - ys, 10)
            xx, yy = np.meshgrid(x, y)
            zz = np.full_like(xx, self.bounder_box[2])
        elif face == 'left':
            y = np.linspace(self.bounder_box[1] + ys, self.bounder_box[4] - ys, 10)
            z = np.linspace(self.bounder_box[2] + zs, self.bounder_box[5] - zs, 10)
            yy, zz = np.meshgrid(y, z)
            xx = np.full_like(yy, self.bounder_box[0])
        elif face == 'right':
            y = np.linspace(self.bounder_box[1] + ys, self.bounder_box[4] - ys, 10)
            z = np.linspace(self.bounder_box[2] + zs, self.bounder_box[5] - zs, 10)
            yy, zz = np.meshgrid(y, z)
            xx = np.full_like(yy, self.bounder_box[3])
        elif face == 'top':
            x = np.linspace(self.bounder_box[0] + xs, self.bounder_box[3] - xs, 10)
            z = np.linspace(self.bounder_box[2] + zs, self.bounder_box[5] - zs, 10)
            xx, zz = np.meshgrid(x, z)
            yy = np.full_like(xx, self.bounder_box[4])
        elif face == 'bottom':
            x = np.linspace(self.bounder_box[0] + xs, self.bounder_box[3] - xs, 10)
            z = np.linspace(self.bounder_box[2] + zs, self.bounder_box[5] - zs, 10)
            xx, zz = np.meshgrid(x, z)
            yy = np.full_like(xx, self.bounder_box[1])
        else:
            raise ValueError("Invalid face name")

        samples = np.stack([xx, yy, zz], axis=-1).reshape(-1, 3)
        # print(f"{face} samples: {samples}")
        return samples

    def compute_max_distance_per_face(self):
        npv = np.array(self.vertices)
        # print(npv)
        tree = KDTree(npv)
        faces = ['front', 'back', 'left', 'right', 'top', 'bottom']
        v_point = {}

        for face in faces:
            samples = self.generate_face_samples(face)
            distances, indexes = tree.query(samples)
            sorted_pairs = sorted(zip(distances, indexes), key=lambda x: x[0], reverse=True)
            d,i = sorted_pairs[int(self.p * 100)]
            v_point[face] = self.vertices[i]

        return v_point

    def cal_cutbox(self):
        # v_farest = self.compute_max_distance_per_face()
        # self.cut_box[0] = v_farest['left'][0]
        # self.cut_box[1] = v_farest['bottom'][1]
        # self.cut_box[2] = v_farest['back'][2]
        # self.cut_box[3] = v_farest['right'][0]
        # self.cut_box[4] = v_farest['top'][1]
        # self.cut_box[5] = v_farest['front'][2]
        self.cut_box[0] = 500
        self.cut_box[1] = 151
        self.cut_box[2] = -2187.5
        self.cut_box[3] = 6900
        self.cut_box[4] = 2165
        self.cut_box[5] = 37.5

    def is_remove(self, cname, shape, vertices, result):
        # 具体条件：检查是否存在至少一个顶点不在任何已接受凸包内
        if self.hull_exist[cname]:
            # 直接跳过凸包未被接受的形状
            for point in vertices:
                # 检测点是否被任何一个已接受的凸包包含
                is_covered = False

                for accepted_name, accepted_hull in self.convex_hulls:
                    if accepted_name == cname:
                        continue
                    if self.is_point_in_3d_convex_hull(point, accepted_hull):
                        is_covered = True
                        break  # 发现包含立即跳出循环
                # 只要有一个点未被覆盖，则保留当前形状
                if not is_covered:
                    if self.show_detail:
                        print(f"{cname}点{point}未被其他凸包覆盖，保留形状")
                    return False

        # 这里对初步确定要删除的形状再加条件判断是否保留
            if (result[3] < self.cut_box[0] or result[4] < self.cut_box[1] or result[5] < self.cut_box[2] or
                    result[0] > self.cut_box[3] or result[1] > self.cut_box[4] or result[2] > self.cut_box[5]):
                print("out")
                # 保留判断条件写在这

        # 删除当前形状
        if self.show_detail:
            print(f"{cname}点全部被其他凸包覆盖，删除形状")
        return True
