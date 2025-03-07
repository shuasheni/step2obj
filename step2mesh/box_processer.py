from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.Bnd import Bnd_Box
import numpy as np
from scipy.spatial import KDTree, QhullError
from step2mesh.stepfile_processer import STEPFileProcessor
from scipy.spatial import ConvexHull

class BounderBoxProcessor(STEPFileProcessor):

    def __init__(self, step_file, lin_deflection=0.2, ang_deflection=0.36, show_detail=False, scale=1.0):
        super().__init__(step_file, lin_deflection=lin_deflection, ang_deflection=ang_deflection, show_detail=show_detail)
        self.scale = scale
        self.cut_box = [9999.0, 9999.0, 9999.0, -9999.0, -9999.0, -9999.0]
        self.bounder_box = [9999.0, 9999.0, 9999.0, -9999.0, -9999.0, -9999.0]
        self.convex_hulls = []

    def is_point_in_3d_convex_hull(self, point, hull, tolerance=1e-10):
        """
        判断点是否在三维凸包内（包括面上）
        :param point: 待测点，形状为 (3,) 的 numpy 数组
        :param hull: scipy.spatial.ConvexHull 对象
        :param tolerance: 浮点误差容忍度
        :return: True（在内部或面上）或 False
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
        accepted_convex_hulls = []

        for shape_name, hull in sorted_convex_hulls:
            # 获取当前凸包的所有顶点
            vertices = hull.points[hull.vertices]
            print(f"{shape_name} hulls points: {vertices}")

            # 检查是否存在至少一个顶点不在任何已接受凸包内
            has_exterior_point = False

            for point in vertices:
                # 检测点是否被任何一个已接受的凸包包含
                is_covered = False
                for accepted_name, accepted_hull in accepted_convex_hulls:
                    if self.is_point_in_3d_convex_hull(point, accepted_hull):
                        is_covered = True
                        break  # 发现包含立即跳出循环
                # 只要有一个点未被覆盖，则保留当前凸包
                if not is_covered:
                    has_exterior_point = True
                    break  # 发现外部点立即停止检测

            if self.show_detail:
                print(f"形状 '{shape_name}' 凸包保留? {has_exterior_point}")

            if has_exterior_point:
                self.convex_hulls.append((shape_name, hull))

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
        return 0

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
        v_farest = {}

        for face in faces:
            samples = self.generate_face_samples(face)
            distances, indexes = tree.query(samples)
            max_distances = 0
            v_index = -1
            for d,i in zip(distances, indexes):
                # print(f"d: {d}, i: {i}, p: {self.vertices[i]}")
                if d > max_distances:
                    max_distances = d
                    v_index = i
            v_farest[face] = self.vertices[v_index]

        return v_farest

    def cal_cutbox(self):
        v_farest = self.compute_max_distance_per_face()
        self.cut_box[0] = v_farest['left'][0]
        self.cut_box[1] = v_farest['bottom'][1]
        self.cut_box[2] = v_farest['back'][2]
        self.cut_box[3] = v_farest['right'][0]
        self.cut_box[4] = v_farest['top'][1]
        self.cut_box[5] = v_farest['front'][2]

    def is_remove(self, cname, shape, vertices, result):
        # 具体条件：检查是否存在至少一个顶点不在任何已接受凸包内

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
        # 所有点被覆盖，则删除当前形状
        if self.show_detail:
            print(f"{cname}点全部被其他凸包覆盖，删除形状")
        return True
