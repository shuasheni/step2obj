from abc import ABC, abstractmethod
import os

from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.STEPCAFControl import STEPCAFControl_Reader
from OCC.Core.TDF import TDF_Label, TDF_LabelSequence
from OCC.Core.TDocStd import TDocStd_Document
from OCC.Core.XCAFDoc import XCAFDoc_DocumentTool
from OCC.Core.TopAbs import TopAbs_REVERSED
from OCC.Core.TopoDS import topods
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.BRepTools import breptools
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh


class STEPFileProcessor(ABC):

    def __init__(self, step_filename="in.step2mesh", lin_deflection=0.2, ang_deflection=0.36, show_detail=False):
        self.step_filename = step_filename
        self.lin_deflection = lin_deflection
        self.ang_deflection = ang_deflection
        self.show_detail = show_detail

        self.doc = None
        self.shape_tool = None
        self.obj = None
        self.vnum = 0
        self.fnum = 0
        self.step_doc_loaded = False
        self.obj_file_opened = False
        self.final_shapes = []
        self.vertices = []

        if not os.path.exists(step_filename):
            print(f"文件 {step_filename} 不存在")
        elif not self.read_file(step_filename):
            print(f"STEP文件 {step_filename} 加载失败")
        else:
            self.step_doc_loaded = True

    def read_file(self, fname):
        self.doc = TDocStd_Document("pythonocc-doc")
        self.shape_tool = XCAFDoc_DocumentTool.ShapeTool(self.doc.Main())

        step_reader = STEPCAFControl_Reader()
        step_reader.SetColorMode(True)
        step_reader.SetLayerMode(True)
        step_reader.SetNameMode(True)
        step_reader.SetMatMode(True)

        status = step_reader.ReadFile(fname)
        if status == IFSelect_RetDone:
            step_reader.Transfer(self.doc)
            return True
        return False

    # @abstractmethod
    def check_shape(self, shape, vertices):
        """计算"""
        pass

    @abstractmethod
    def is_remove(self, cname, shape, vertices, result):
        """由子类实现的移除判断条件"""
        pass

    def shape_to_obj(self, shape):
        topo_exp = TopologyExplorer(shape)
        all_vertices = []
        all_face_meshes = []
        vertex_offset = 0
        for i, face in enumerate(topo_exp.faces()):
            breptools.Clean(face)
            BRepMesh_IncrementalMesh(face, self.lin_deflection, True, self.ang_deflection).Perform()
            rev_flag = False
            if face.Orientation() == TopAbs_REVERSED:
                rev_flag = True

            loc = TopLoc_Location()
            face_triangulation = BRep_Tool.Triangulation(topods.Face(face), loc)
            if face_triangulation is None:
                continue

            vertices = []
            face_mesh = []
            for j in range(1, face_triangulation.NbNodes() + 1):
                vertex = face_triangulation.Node(j)
                vertices.append([vertex.X(), vertex.Y(), vertex.Z()])


            for j in range(1, face_triangulation.NbTriangles() + 1):
                triangle = face_triangulation.Triangle(j)
                if rev_flag:
                    v1_idx = triangle.Value(3) - 1
                    v2_idx = triangle.Value(2) - 1
                    v3_idx = triangle.Value(1) - 1
                else:
                    v1_idx = triangle.Value(1) - 1
                    v2_idx = triangle.Value(2) - 1
                    v3_idx = triangle.Value(3) - 1
                face_vertex_indices = [v1_idx + vertex_offset, v2_idx + vertex_offset, v3_idx + vertex_offset]
                face_mesh.append(face_vertex_indices)

            all_vertices.extend(vertices)
            all_face_meshes.append(face_mesh)
            vertex_offset += len(vertices)

        self.vertices.extend(all_vertices)
        return all_vertices, all_face_meshes

    def add_obj_to_file(self, shape_name, all_vertices, all_face_meshes):
        assert self.obj_file_opened, "未打开OBJ文件"

        self.obj.write("#\n")
        self.obj.write(f"# {shape_name}\n")
        self.obj.write("#\n")
        for v in all_vertices:
            self.obj.write(f"v {v[0]} {v[1]} {v[2]}\n")
        self.obj.write("\n")

        self.obj.write(f"g {shape_name}\n")
        for face_mesh in all_face_meshes:
            for f_verts in face_mesh:
                self.obj.write(f"f {f_verts[0] + 1 + self.vnum}// {f_verts[1] + 1 + self.vnum}// {f_verts[2] + 1 + self.vnum}//\n")
                self.fnum += 1
            self.obj.write("\n")
        self.vnum += len(all_vertices)
        if self.show_detail:
            print(f"Conversion complete: {shape_name}")

    def travel_components(self, comps, parent_loc=None):
        for j in range(comps.Length()):
            c_label = comps.Value(j + 1)

            # 获取当前节点的坐标系变换
            child_loc = self.shape_tool.GetLocation(c_label)
            child_trsf = child_loc.Transformation()
            current_trsf = parent_loc.Transformation()
            current_trsf.Multiply(child_trsf)

            # 处理引用（引用实例或直接形状）
            ref_label = TDF_Label()
            if self.shape_tool.GetReferredShape(c_label, ref_label):
                target_label = ref_label
            else:
                target_label = c_label

            # 判断是否为装配体并递归处理
            if self.shape_tool.IsAssembly(target_label):
                sub_comps = TDF_LabelSequence()
                self.shape_tool.GetComponents(target_label, sub_comps)  # 仅取直接子组件
                self.travel_components(sub_comps, TopLoc_Location(current_trsf))  # 传递累积的坐标系变换
            else:
                # 处理最终形状
                c_name = c_label.GetLabelName()
                shape = self.shape_tool.GetShape(target_label)

                analyzer = BRepCheck_Analyzer(shape)
                if not analyzer.IsValid():
                    if self.show_detail:
                        print(f"not valid: {c_name}")
                    continue
                if self.show_detail:
                    print(f"travel: {c_name}:")
                transformed_shape = BRepBuilderAPI_Transform(shape, current_trsf, True).Shape()
                all_vertices, all_face_meshes = self.shape_to_obj(transformed_shape)
                result = self.check_shape(transformed_shape, all_vertices)
                shape_info = (c_name, transformed_shape, all_vertices, all_face_meshes, result)
                self.final_shapes.append(shape_info)

    def travel(self):
        """遍历全部形状进行计算并转化为mesh"""
        assert self.step_doc_loaded, "未加载STEP文件"
        self.final_shapes = []

        labels = TDF_LabelSequence()
        self.shape_tool.GetShapes(labels)
        rootlabel = labels.Value(1)  # 1号形状标签代表整个装配体
        if self.shape_tool.IsAssembly(rootlabel):
            top_comps = TDF_LabelSequence()
            subchilds = False
            is_assy = self.shape_tool.GetComponents(rootlabel, top_comps, subchilds)
            top_loc = self.shape_tool.GetLocation(rootlabel)
            if is_assy and top_comps.Length():
                self.travel_components(top_comps, top_loc)
                return True
        return False

    def output(self, obj_filename="out.obj"):
        with open(obj_filename, 'w') as obj:
            self.obj = obj
            self.obj_file_opened = True
            self.vnum = 0
            self.fnum = 0
            for (c_name, transformed_shape, all_vertices, all_face_meshes, result) in self.final_shapes:
                if not self.is_remove(c_name, transformed_shape, all_vertices, result):
                    self.add_obj_to_file(c_name, all_vertices, all_face_meshes)

        self.obj_file_opened = False
        print(f"简化MESH输出完成, 面数量:{self.fnum}")
        return self.fnum

    def process_components(self, comps, parent_loc=None):
        for j in range(comps.Length()):
            c_label = comps.Value(j + 1)

            # 获取当前节点的坐标系变换
            child_loc = self.shape_tool.GetLocation(c_label)
            child_trsf = child_loc.Transformation()
            current_trsf = parent_loc.Transformation()
            current_trsf.Multiply(child_trsf)

            # 处理引用（引用实例或直接形状）
            ref_label = TDF_Label()
            if self.shape_tool.GetReferredShape(c_label, ref_label):
                target_label = ref_label
            else:
                target_label = c_label

            # 判断是否为装配体并递归处理
            if self.shape_tool.IsAssembly(target_label):
                sub_comps = TDF_LabelSequence()
                self.shape_tool.GetComponents(target_label, sub_comps)  # 仅取直接子组件
                self.process_components(sub_comps, TopLoc_Location(current_trsf))  # 传递累积的坐标系变换
            else:
                # 处理最终形状
                c_name = c_label.GetLabelName()
                shape = self.shape_tool.GetShape(target_label)

                analyzer = BRepCheck_Analyzer(shape)
                if not analyzer.IsValid():
                    if self.show_detail:
                        print(f"not valid: {c_name}")
                    continue
                if self.show_detail:
                    print(f"process: {c_name}:")
                transformed_shape = BRepBuilderAPI_Transform(shape, current_trsf, True).Shape()
                all_vertices, all_face_meshes = self.shape_to_obj(transformed_shape)

                if not self.is_remove(transformed_shape, all_vertices):
                    self.add_obj_to_file(c_name, all_vertices, all_face_meshes)

    def pocess(self, obj_filename="out.obj"):
        """移除部分形状并将剩余形状转化为OBJ输出"""
        assert self.step_doc_loaded, "未加载STEP文件"

        labels = TDF_LabelSequence()
        self.shape_tool.GetShapes(labels)
        rootlabel = labels.Value(1)  # 1号形状标签代表整个装配体
        if self.shape_tool.IsAssembly(rootlabel):
            top_comps = TDF_LabelSequence()
            subchilds = False
            is_assy = self.shape_tool.GetComponents(rootlabel, top_comps, subchilds)
            top_loc = self.shape_tool.GetLocation(rootlabel)
            if is_assy and top_comps.Length():
                with open(obj_filename, 'w') as obj:
                    self.obj = obj
                    self.obj_file_opened = True
                    self.vnum = 0
                    self.fnum = 0
                    self.process_components(top_comps, top_loc)
                self.obj_file_opened = False
                print(f"简化MESH输出完成, 面数量:{self.fnum}")
                return self.fnum
        return 0